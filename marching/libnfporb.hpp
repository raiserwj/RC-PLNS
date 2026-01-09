#ifndef NFP_HPP_
#define NFP_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <exception>
#include <limits>

#include "boost/geometry/algorithms/intersects.hpp"

#include "geometry.hpp"
#include "svg.hpp"
#include "wkt.hpp"
#include "translation_vector.hpp"
#include "history.hpp"
#include "marching/algo/touching_point.hpp"
#include "marching/algo/trim_vector.hpp"
#include "marching/algo/select_next.hpp"
#include "marching/algo/find_feasible.hpp"
#include "marching/algo/search_start.hpp"
#include "marching/algo/slide.hpp"
#include "vector"
#include <iomanip>

using namespace std;
namespace libnfporb {

/**
 * Delete oscillations and loops from the given ring.
 * @param ring A reference to the ring to be shortened.
 * @return true if the ring has changed.
 */
    inline bool delete_consecutive_repeating_point_patterns(polygon_t::ring_type &ring) {
        size_t startLen = ring.size();
        off_t len = ring.size();
        int i, j, counter;
        for (i = 1; i <= len / 2; ++i) {
            for (j = i, counter = 0; j < len; ++j) {
                if (equals(ring[j], ring[j - i]))
                    counter++;
                else
                    counter = 0;
                if (counter > 2 && counter == i) {
                    counter = 0;
                    std::copy(ring.begin() + j, ring.begin() + len, ring.begin() + (j - i));
                    j -= i;
                    len -= i;
                }
            }
            ring.resize(j);
        }

        size_t start = 0, cnt = 0;
        point_t c, l = ring[0];

        for (size_t i = 1; i < ring.size(); ++i) {
            c = ring[i];
            if (equals(c, l)) {
                if (cnt == 0)
                    start = i - 1;

                ++cnt;
            } else {
                if (cnt > 1) {
                    ring.erase(ring.begin() + start + 1, ring.begin() + start + cnt);
                    if (start + cnt >= ring.size())
                        break;
                }
            }
            l = c;
        }
        return ring.size() != startLen;
    }

/**
 * Remove co-linear points from a ring.
 * @param r A reference to a ring.
 */
    inline void remove_co_linear(polygon_t::ring_type &r) {
        assert(r.size() > 2);
        psize_t nextI;
        psize_t prevI = 0;
        segment_t segment(r[r.size() - 2], r[0]);
        polygon_t::ring_type newR;

        for (psize_t i = 1; i < r.size() + 1; ++i) {
            if (i >= r.size())
                nextI = i % r.size() + 1;
            else
                nextI = i;

            if (get_alignment(segment, r[nextI]) != ON) {
                newR.push_back(r[prevI]);
            }
            segment = {segment.second, r[nextI]};
            prevI = nextI;
        }

        r = newR;
    }

/**
 * Remove · points from a polygon.
 * @param p A reference to a polygon.
 */
    inline void remove_co_linear(polygon_t &p) {
        remove_co_linear(p.outer());
        for (auto &r: p.inners())
            remove_co_linear(r);

        bg::correct(p);
    }

/**
 * Generate the NFP for the given polygons pA and pB. Optionally check the input polygons for validity.
 * @param pA polygon A (the stationary polygon).
 * @param pB polygon B (the orbiting polygon).
 *  * @param checkValidity Check the input polygons for validity if true. Defaults to true.
 * @return The generated NFP.
 */
    inline nfp_t generate_nfp(polygon_t &pA, polygon_t &pB, const bool checkValidity = true) {
        remove_co_linear(pA);
        remove_co_linear(pB);

        if (checkValidity) {
            std::string reason;
            if (!bg::is_valid(pA, reason))
                throw std::runtime_error("Polygon A is invalid: " + reason);

            if (!bg::is_valid(pB, reason))
                throw std::runtime_error("Polygon B is invalid: " + reason);
        }

        nfp_t nfp;

#ifdef NFP_DEBUG
        write_svg("start.svg", {pA, pB});
#endif

        DEBUG_VAL(bg::wkt(pA));
        DEBUG_VAL(bg::wkt(pB));

        //prevent double vertex connections at start because we might come back the same way we go which would end the nfp prematurely
        std::vector<psize_t> yAminI = find_minimum_y(pA);
        std::vector<psize_t> yBmaxI = find_maximum_y(pB);

        point_t pAstart;
        point_t pBstart;

        if (yAminI.size() > 1 || yBmaxI.size() > 1) {
            //find right-most of A and left-most of B to prevent double connection at start
            coord_t maxX = MIN_COORD;
            psize_t iRightMost = 0;
            for (psize_t &ia: yAminI) {
                const point_t &candidateA = pA.outer()[ia];
                if (larger(candidateA.x_, maxX)) {
                    maxX = candidateA.x_;
                    iRightMost = ia;
                }
            }

            coord_t minX = MAX_COORD;
            psize_t iLeftMost = 0;
            for (psize_t &ib: yBmaxI) {
                const point_t &candidateB = pB.outer()[ib];
                if (smaller(candidateB.x_, minX)) {
                    minX = candidateB.x_;
                    iLeftMost = ib;
                }
            }
            pAstart = pA.outer()[iRightMost];
            pBstart = pB.outer()[iLeftMost];
        } else {
            pAstart = pA.outer()[yAminI.front()];
            pBstart = pB.outer()[yBmaxI.front()];
        }

        nfp.push_back({});
        point_t transB = {pAstart - pBstart};


        SlideResult res;
        if ((res = slide(pA, pA.outer(), pB.outer(), nfp, transB, false)) != LOOP) {
            throw std::runtime_error("Unable to complete outer nfp loop: " + std::to_string(res));
        }

        DEBUG_VAL("##### outer #####");
        point_t startTrans;
        while (true) {
            SearchStartResult res = search_start_translation(pA.outer(), pB.outer(), nfp, false, startTrans);
            if (res == FOUND) {
                nfp.push_back({});
                DEBUG_VAL("##### interlock start #####");
                polygon_t::ring_type rifsB;
                boost::geometry::transform(pB.outer(), rifsB,
                                           trans::translate_transformer<coord_t, 2, 2>(startTrans.x_, startTrans.y_));
                if (in_nfp(rifsB.front(), nfp)) {
                    continue;
                }
                SlideResult sres = slide(pA, pA.outer(), pB.outer(), nfp, startTrans, true);
                if (sres != LOOP) {
                    if (sres == NO_TRANSLATION) {
                        //no initial slide found -> jigsaw
                        if (!in_nfp(pB.outer().front(), nfp)) {
                            nfp.push_back({});
                            nfp.back().push_back(pB.outer().front());
                        }
                    }
                }
                DEBUG_VAL("##### interlock end #####");
            } else if (res == FIT) {
                DEBUG_VAL("##### perfect fit #####");
                point_t reference = pB.outer().front();
                point_t translated;
                trans::translate_transformer<coord_t, 2, 2> translate(startTrans.x_, startTrans.y_);
                boost::geometry::transform(reference, translated, translate);
                if (!in_nfp(translated, nfp)) {
                    nfp.push_back({});
                    nfp.back().push_back(translated);
                }
                break;
            } else
                break;
        }

        for (auto &rA: pA.inners()) {
            while (true) {
                SearchStartResult res = search_start_translation(rA, pB.outer(), nfp, true, startTrans);
                if (res == FOUND) {
                    nfp.push_back({});
                    DEBUG_VAL("##### hole start #####");
                    slide(pA, rA, pB.outer(), nfp, startTrans, true);
                    DEBUG_VAL("##### hole end #####");
                } else if (res == FIT) {
                    point_t reference = pB.outer().front();
                    point_t translated;
                    trans::translate_transformer<coord_t, 2, 2> translate(startTrans.x_, startTrans.y_);
                    boost::geometry::transform(reference, translated, translate);
                    if (!in_nfp(translated, nfp)) {
                        nfp.push_back({});
                        nfp.back().push_back(translated);
                    }
                    break;
                } else
                    break;
            }
        }

#ifdef NFP_DEBUG
        write_svg("nfp.svg", pA, pB, nfp);
#endif

        for (auto &r: nfp) {
            while (delete_consecutive_repeating_point_patterns(r));
            bg::correct(r);
        }
        return nfp;
    }
}

class nfpClass {
public:
    vector<float> nfpCalculate(float *A, int nA, float *B, int nB, float width, float height) {
        libnfporb::polygon_t a, b;
        float minxA = 100000, maxxA = -100000, minyA = 100000, maxyA = -100000;
        float minxB = 100000, maxxB = -100000, minyB = 100000, maxyB = -100000;
        for (int i = 0; i < nA; i++) {
            bg::append(a.outer(), libnfporb::point_t(A[i * 2] + 10000, A[i * 2 + 1] + 10000));
            if (A[i * 2] < minxA) minxA = A[i * 2];
            if (A[i * 2] > maxxA) maxxA = A[i * 2];
            if (A[i * 2 + 1] < minyA) minyA = A[i * 2 + 1];
            if (A[i * 2 + 1] > maxyA) maxyA = A[i * 2 + 1];
        }
        for (int i = 0; i < nB; i++) {
            bg::append(b.outer(), libnfporb::point_t(B[i * 2] + 10000, B[i * 2 + 1] + 10000));
//            std::cout<<B[i * 2]<<B[i * 2 + 1]<<std::endl;
            if (B[i * 2] < minxB) minxB = B[i * 2];
            if (B[i * 2] > maxxB) maxxB = B[i * 2];
            if (B[i * 2 + 1] < minyB) minyB = B[i * 2 + 1];
            if (B[i * 2 + 1] > maxyB) maxyB = B[i * 2 + 1];
        }
        float origArea = (maxxA - minxA) * (maxyA - minyA) + (maxxB - minxB) * (maxyB - minyB);
        std::cout << fixed << setprecision(2) << maxxA - minxA << "  " << maxyA - minyA << "  " << maxxB - minxB << "  " << maxyB - minyB << std::endl;
        libnfporb::polygon_t b_simple, a_simple;
        boost::geometry::simplify(b, b_simple, 10);
        boost::geometry::simplify(a, a_simple, 10);
        float bSimpleX = b_simple.outer()[0].x_.val() - 10000;
        float bSimpleY = b_simple.outer()[0].y_.val() - 10000;
        minxB = minxB - bSimpleX;
        maxxB = maxxB - bSimpleX;
        minyB = minyB - bSimpleY;
        maxyB = maxyB - bSimpleY;

        libnfporb::nfp_t nfpv = generate_nfp(a_simple, b_simple, true);


        vector<float> bestLocation = {0, 0, 0, 0, 0, 0, 0};
        float bestReward = 10;
        float bestDistance = 1000000;
        for (int i = 0; i < nfpv[0].size(); i++) {
            float reward, tempDistance;
            float tempX = nfpv[0][i].x_.val() - 10000;
            float tempY = nfpv[0][i].y_.val() - 10000;

            if (max(maxxA, tempX + maxxB) - min(minxA, tempX + minxB) > width || max(maxyA, tempY + maxyB) - min(minyA, tempY + minyB) > height) continue;
            reward = origArea - (max(maxxA, tempX + maxxB) - min(minxA, tempX + minxB)) *
                                (max(maxyA, tempY + maxyB) - min(minyA, tempY + minyB));

            if (reward > bestReward || (reward == bestReward && tempDistance < bestDistance)) {
                bestReward = reward;
                bestDistance = tempDistance;
                bestLocation[0] = reward;
                bestLocation[1] = max(maxxA, tempX + maxxB) - min(minxA, tempX + minxB);
                bestLocation[2] = max(maxyA, tempY + maxyB) - min(minyA, tempY + minyB);
                float xMiddle = (max(maxxA, tempX + maxxB) + min(minxA, tempX + minxB)) / 2;
                float yMiddle = (max(maxyA, tempY + maxyB) + min(minyA, tempY + minyB)) / 2;
                bestLocation[3] = (minxA + maxxA) / 2 - xMiddle;
                bestLocation[4] = (minyA + maxyA) / 2 - yMiddle;
                bestLocation[5] = (tempX + minxB + tempX + maxxB) / 2 - xMiddle;
                bestLocation[6] = (tempY + minyB + tempY + maxyB) / 2 - yMiddle;
            }
        }
        return bestLocation;
    }

    libnfporb::polygon_t MakeBox(float xmin, float ymin, float xmax, float ymax) {
        libnfporb::polygon_t poly;
        poly.outer().assign({
                                    {xmin, ymin},
                                    {xmax, ymin},
                                    {xmax, ymax},
                                    {xmin, ymax}
                            });
        return poly;
    }

    vector<vector<float>> nfpCalculate_(vector<vector<vector<float>>> A, vector<vector<float>> B, float width, float height) {
        vector<vector<float>> ret;
        vector<libnfporb::polygon_t> a;
        vector<vector<float>> boundingboxA;
        libnfporb::polygon_t b;
        float minxB = 100000, maxxB = -100000, minyB = 100000, maxyB = -100000;
        float minxA = 100000, maxxA = -100000, minyA = 100000, maxyA = -100000;
        for (int i = 0; i < A.size(); i++) {
            libnfporb::polygon_t tempA, tempASimple;
            float minxATemp = 100000, maxxATemp = -100000, minyATemp = 100000, maxyATemp = -100000;
            for (int j = 0; j < A[i].size(); j++) {
                bg::append(tempA.outer(), libnfporb::point_t(A[i][j][0] + 10000, A[i][j][1] + 10000));
                if (A[i][j][0] < minxATemp) minxATemp = A[i][j][0];
                if (A[i][j][0] > maxxATemp) maxxATemp = A[i][j][0];
                if (A[i][j][1] < minyATemp) minyATemp = A[i][j][1];
                if (A[i][j][1] > maxyATemp) maxyATemp = A[i][j][1];
            }
            bg::append(tempA.outer(), libnfporb::point_t(A[i][0][0] + 10000, A[i][0][1] + 10000));
            boundingboxA.push_back({minxATemp, maxxATemp, minyATemp, maxyATemp});
            if (minxATemp < minxA) minxA = minxATemp;
            if (maxxATemp > maxxA) maxxA = maxxATemp;
            if (minyATemp < minyA) minyA = minyATemp;
            if (maxyATemp > maxyA) maxyA = maxyATemp;
            boost::geometry::simplify(tempA, tempASimple, 10);
            a.push_back(tempASimple);
        }
        for (int i = 0; i < B.size(); i++) {
            bg::append(b.outer(), libnfporb::point_t(B[i][0] + 10000, B[i][1] + 10000));
            if (B[i][0] < minxB) minxB = B[i][0];
            if (B[i][0] > maxxB) maxxB = B[i][0];
            if (B[i][1] < minyB) minyB = B[i][1];
            if (B[i][1] > maxyB) maxyB = B[i][1];
        }
        bg::append(b.outer(), libnfporb::point_t(B[0][0] + 10000, B[0][1] + 10000));

        float origArea = (maxxA - minxA) * (maxyA - minyA) + (maxxB - minxB) * (maxyB - minyB);
        libnfporb::polygon_t b_simple;
        boost::geometry::simplify(b, b_simple, 10);
        float bSimpleX = b_simple.outer()[0].x_.val() - 10000;
        float bSimpleY = b_simple.outer()[0].y_.val() - 10000;
        minxB = minxB - bSimpleX;
        maxxB = maxxB - bSimpleX;
        minyB = minyB - bSimpleY;
        maxyB = maxyB - bSimpleY;

        vector<libnfporb::polygon_t> collection;
        collection.push_back(MakeBox(-100000, -100000, 100000, 100000));
        for (int i = 0; i < A.size(); i++) {
            libnfporb::nfp_t nfp = generate_nfp(a[i], b_simple, true);

            vector<libnfporb::polygon_t> tempCollection;
            for (auto polygon: collection) {
                vector<libnfporb::polygon_t> tempCollection_;
                boost::geometry::difference(polygon, nfp[0], tempCollection_);
                tempCollection.insert(tempCollection.end(), tempCollection_.begin(), tempCollection_.end());
            }
            collection = tempCollection;
        }
        vector<float> bestLocation;
        float bestReward = 10, bestDistance = 1000000;
        for (auto polygon: collection) {
            vector<libnfporb::point_t> points;
            for (auto point: polygon.outer()) {
                if (point.x_.val() == -100000 || point.x_.val() == 100000 || point.y_.val() == -100000 || point.y_.val() == 100000) continue;
                points.push_back(point);
            }
            for (auto ring: polygon.inners()) {
                for (auto point: ring) {
                    if (point.x_.val() == -100000 || point.x_.val() == 100000 || point.y_.val() == -100000 || point.y_.val() == 100000) continue;
                    points.push_back(point);
                }
            }
            for (auto point: points) {
                float reward;
                float distance;
                float tempX = point.x_.val() - 10000;
                float tempY = point.y_.val() - 10000;
                if (max(maxxA, tempX + maxxB) - min(minxA, tempX + minxB) > width || max(maxyA, tempY + maxyB) - min(minyA, tempY + minyB) > height) continue;
                reward = origArea - (max(maxxA, tempX + maxxB) - min(minxA, tempX + minxB)) *
                                    (max(maxyA, tempY + maxyB) - min(minyA, tempY + minyB));
                float xMiddle = (max(maxxA, tempX + maxxB) + min(minxA, tempX + minxB)) / 2;
                float yMiddle = (max(maxyA, tempY + maxyB) + min(minyA, tempY + minyB)) / 2;
                float xBMiddle = (tempX + maxxB + tempX + minxB) / 2, yBMiddle = (tempY + maxyB + tempY + minyB) / 2;
                distance = pow(pow(xBMiddle - xMiddle, 2) + pow(yBMiddle - yMiddle, 2), 0.5);
                if (reward > bestReward || (reward == bestReward && distance < bestDistance)) { ;
                    bestReward = reward;
                    bestDistance = distance;

                    ret.clear();
                    ret.push_back({reward, max(maxxA, tempX + maxxB) - min(minxA, tempX + minxB), max(maxyA, tempY + maxyB) - min(minyA, tempY + minyB), tempX - bSimpleX, tempY - bSimpleY});
                    for (int i = 0; i < A.size(); i++) {
                        vector<float> location;
                        location.push_back((boundingboxA[i][0] + boundingboxA[i][1]) / 2 - xMiddle);
                        location.push_back((boundingboxA[i][2] + boundingboxA[i][3]) / 2 - yMiddle);
                        ret.push_back(location);
                    }
                    vector<float> location;
                    location.push_back(xBMiddle - xMiddle);
                    location.push_back(yBMiddle - yMiddle);
                    ret.push_back(location);
                }
            }
        }
        return ret;
    }
    vector<vector<float>> nfpbetweentwo(vector<vector<vector<float>>> A, vector<vector<float>> B, float width, float height) {
        vector<vector<float>> ret={};
        vector<libnfporb::polygon_t> a;
        vector<vector<float>> boundingboxA;
        libnfporb::polygon_t b;
        float minxB = 100000, maxxB = -100000, minyB = 100000, maxyB = -100000;
        float minxA = 100000, maxxA = -100000, minyA = 100000, maxyA = -100000;
        for (int i = 0; i < A.size(); i++) {
            libnfporb::polygon_t tempA, tempASimple;
            float minxATemp = 100000, maxxATemp = -100000, minyATemp = 100000, maxyATemp = -100000;
            for (int j = 0; j < A[i].size(); j++) {
                bg::append(tempA.outer(), libnfporb::point_t(A[i][j][0] +10000, A[i][j][1] +10000));
                if (A[i][j][0] < minxATemp) minxATemp = A[i][j][0];
                if (A[i][j][0] > maxxATemp) maxxATemp = A[i][j][0];
                if (A[i][j][1] < minyATemp) minyATemp = A[i][j][1];
                if (A[i][j][1] > maxyATemp) maxyATemp = A[i][j][1];
            }
            bg::append(tempA.outer(), libnfporb::point_t(A[i][0][0]+10000, A[i][0][1]+10000));
            boundingboxA.push_back({minxATemp, maxxATemp, minyATemp, maxyATemp});
            if (minxATemp < minxA) minxA = minxATemp;
            if (maxxATemp > maxxA) maxxA = maxxATemp;
            if (minyATemp < minyA) minyA = minyATemp;
            if (maxyATemp > maxyA) maxyA = maxyATemp;
            boost::geometry::simplify(tempA, tempASimple, 0);
            a.push_back(tempASimple);
        }
        for (int i = 0; i < B.size(); i++) {
            bg::append(b.outer(), libnfporb::point_t(B[i][0]+10000, B[i][1]+10000));
            if (B[i][0] < minxB) minxB = B[i][0];
            if (B[i][0] > maxxB) maxxB = B[i][0];
            if (B[i][1] < minyB) minyB = B[i][1];
            if (B[i][1] > maxyB) maxyB = B[i][1];
        }
        bg::append(b.outer(), libnfporb::point_t(B[0][0]+10000 , B[0][1]+10000 ));

        float origArea = (maxxA - minxA) * (maxyA - minyA) + (maxxB - minxB) * (maxyB - minyB);
        libnfporb::polygon_t b_simple;
        boost::geometry::simplify(b, b_simple, 0);
        float bSimpleX = b_simple.outer()[0].x_.val() ;
        float bSimpleY = b_simple.outer()[0].y_.val() ;
        minxB = minxB - bSimpleX;
        maxxB = maxxB - bSimpleX;
        minyB = minyB - bSimpleY;
        maxyB = maxyB - bSimpleY;

        vector<libnfporb::polygon_t> collection;
        collection.push_back(MakeBox(-100000, -100000, 100000, 100000));
        for (int i = 0; i < A.size(); i++) {
            libnfporb::nfp_t nfp = generate_nfp(a[i], b_simple, true);

            vector<libnfporb::polygon_t> tempCollection;
            for (auto polygon: collection) {
                vector<libnfporb::polygon_t> tempCollection_;
                boost::geometry::difference(polygon, nfp[0], tempCollection_);
                tempCollection.insert(tempCollection.end(), tempCollection_.begin(), tempCollection_.end());
            }
            collection = tempCollection;
            for(int i=0;i<nfp[0].size();i++)
            {
                float x=nfp[0][i].x_.val()-10000.0-B[0][0];
                float y=nfp[0][i].y_.val()-10000.0-B[0][1];
                ret.push_back({x,y});
            }
        }

//        vector<float> bestLocation;
//        float bestReward = 10, bestDistance = 1000000;
//        for (auto polygon: collection) {
//            vector<libnfporb::point_t> points;
//            for (auto point: polygon.outer()) {
//                if (point.x_.val() == -100000 || point.x_.val() == 100000 || point.y_.val() == -100000 || point.y_.val() == 100000) continue;
//                points.push_back(point);
//            }
//            for (auto ring: polygon.inners()) {
//                for (auto point: ring) {
//                    if (point.x_.val() == -100000 || point.x_.val() == 100000 || point.y_.val() == -100000 || point.y_.val() == 100000) continue;
//                    points.push_back(point);
//                }
//            }
//            for (auto point: points) {
//                float reward;
//                float distance;
//                float tempX = point.x_.val() - 10000;
//                float tempY = point.y_.val() - 10000;
//                if (max(maxxA, tempX + maxxB) - min(minxA, tempX + minxB) > width || max(maxyA, tempY + maxyB) - min(minyA, tempY + minyB) > height) continue;
//                reward = origArea - (max(maxxA, tempX + maxxB) - min(minxA, tempX + minxB)) *
//                                    (max(maxyA, tempY + maxyB) - min(minyA, tempY + minyB));
//                float xMiddle = (max(maxxA, tempX + maxxB) + min(minxA, tempX + minxB)) / 2;
//                float yMiddle = (max(maxyA, tempY + maxyB) + min(minyA, tempY + minyB)) / 2;
//                float xBMiddle = (tempX + maxxB + tempX + minxB) / 2, yBMiddle = (tempY + maxyB + tempY + minyB) / 2;
//                distance = pow(pow(xBMiddle - xMiddle, 2) + pow(yBMiddle - yMiddle, 2), 0.5);
//                if (reward > bestReward || (reward == bestReward && distance < bestDistance)) { ;
//                    bestReward = reward;
//                    bestDistance = distance;
//
//                    ret.clear();
//                    ret.push_back({reward, max(maxxA, tempX + maxxB) - min(minxA, tempX + minxB), max(maxyA, tempY + maxyB) - min(minyA, tempY + minyB), tempX - bSimpleX, tempY - bSimpleY});
//                    for (int i = 0; i < A.size(); i++) {
//                        vector<float> location;
//                        location.push_back((boundingboxA[i][0] + boundingboxA[i][1]) / 2 - xMiddle);
//                        location.push_back((boundingboxA[i][2] + boundingboxA[i][3]) / 2 - yMiddle);
//                        ret.push_back(location);
//                    }
//                    vector<float> location;
//                    location.push_back(xBMiddle - xMiddle);
//                    location.push_back(yBMiddle - yMiddle);
//                    ret.push_back(location);
//                }
//            }
//        }
        return ret;
    }
};

#endif
