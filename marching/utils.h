//
// Created by raiser on 2022/4/26.
//

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include<vector>
#include<math.h>
#include "libnfporb.hpp"
#include "Matching.h"

using std::vector;

inline float dis(vector<float> p1, vector<float> p2) {
    return (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]);
}

inline bool isrect(vector<vector<float>> polygon) {
    float x_c = (polygon[0][0] + polygon[1][0] + polygon[2][0] + polygon[3][0]) / 4;
    float y_c = (polygon[0][1] + polygon[1][1] + polygon[2][1] + polygon[3][1]) / 4;
    float d1 = dis(polygon[0], {x_c, y_c});
    float d2 = dis(polygon[1], {x_c, y_c});
    float d3 = dis(polygon[2], {x_c, y_c});
    float d4 = dis(polygon[3], {x_c, y_c});
    if (d1 == d2 and d2 == d3 and d3 == d4) {
        return true;
    } else { return false; }
}

inline vector<vector<float>> rotate_polygon(vector<vector<float>> polygon, int angle, float x_, float y_) {
    vector<vector<float>> rotated = {};
    float cos_, sin_;
    switch (angle) {
        case 0:
            cos_ = 1;
            sin_ = 0;
            break;
        case 90:
            cos_ = 0;
            sin_ = 1;
            break;
        case 180:
            cos_ = -1;
            sin_ = 0;
            break;
        case 270:
            cos_ = 0;
            sin_ = -1;
            break;
    }
    float x;
    float y;
    for (int i = 0; i < polygon.size(); i++) {
        x = polygon[i][0];
        y = polygon[i][1];
        rotated.push_back({x * cos_ - y * sin_ + x_, x * sin_ + y * cos_ + y_});
    }
    return rotated;
}

inline void rotate_polygon_(vector<vector<float>> &polygon, int angle) {
    float cos_, sin_;
    switch (angle) {
        case 0:
            cos_ = 1;
            sin_ = 0;
            break;
        case 90:
            cos_ = 0;
            sin_ = 1;
            break;
        case 180:
            cos_ = -1;
            sin_ = 0;
            break;
        case 270:
            cos_ = 0;
            sin_ = -1;
            break;
    }
    float x;
    float y;
    for (int i = 0; i < polygon.size(); i++) {
        x = polygon[i][0];
        y = polygon[i][1];
        polygon[i] = {x * cos_ - y * sin_, x * sin_ + y * cos_};
    }
}

inline void shift_polygon_(vector<vector<float>> &polygon, float x, float y) {
    for (int i = 0; i < polygon.size(); i++) {
        polygon[i] = {polygon[i][0] + x, polygon[i][1] + y};
    }
}

inline void rotate_polygons_(vector<vector<vector<float>>> &polygons, int angle) {
    float cos_, sin_;
    switch (angle) {
        case 0:
            cos_ = 1;
            sin_ = 0;
            break;
        case 90:
            cos_ = 0;
            sin_ = 1;
            break;
        case 180:
            cos_ = -1;
            sin_ = 0;
            break;
        case 270:
            cos_ = 0;
            sin_ = -1;
            break;
    }
    float x;
    float y;
    for (int i = 0; i < polygons.size(); i++) {
        for (int j = 0; j < polygons[i].size(); j++) {
            x = polygons[i][j][0];
            y = polygons[i][j][1];
            polygons[i][j] = {x * cos_ - y * sin_, x * sin_ + y * cos_};
        }
    }
}


inline vector<vector<float>> group_two_(vector<vector<vector<float>>> A, vector<vector<float>> B, vector<int> rotateA, vector<int> rotateB, float width, float height, long int &totaltime) {
    vector<vector<float>> bestret;
    float bestReward = 0;
    nfpClass nfpC;
//    for (int i = 0; i < rotateA.size(); i++) {
//        for (int j = 0; j < rotateB.size(); j++) {
//
//        }
//    }

    for (int i = 0; i < rotateA.size(); i++) {
        for (int j = 0; j < rotateB.size(); j++) {
            vector<vector<vector<float>>> a;
            for (int k = 0; k < A.size(); k++) {
                a.push_back(rotate_polygon(A[k], rotateA[i], 0, 0));
            }
            vector<vector<float>> b = rotate_polygon(B, rotateB[j], 0, 0);
            vector<vector<float>> ret = nfpC.nfpCalculate_(a, b, width, height,totaltime);
            if (ret.size() != 0 && ret[0][0] > bestReward) {
                bestReward = ret[0][0];
                bestret.clear();
                bestret.insert(bestret.end(), ret.begin(), ret.end());
                bestret[0].push_back(rotateA[i]);
                bestret[0].push_back(rotateB[j]);
            }
        }
    }
    return bestret;
}

inline float squarearea(vector<vector<float>> polygon) {
    float xmax = -99999;
    float ymax = -99999;
    float xmin = 99999;
    float ymin = 99999;
    for (int i = 0; i < polygon.size(); i++) {
        if (polygon[i][0] > xmax) {
            xmax = polygon[i][0];
        }
        if (polygon[i][0] < xmin) {
            xmin = polygon[i][0];
        }
        if (polygon[i][1] > ymax) {
            ymax = polygon[i][1];
        }
        if (polygon[i][1] < ymin) {
            ymin = polygon[i][1];
        }
    }
    return (xmax - xmin) * (ymax - ymin);
}

inline float polygonarea(vector<vector<float>> polygon) {
    float area = 0;
    vector<float> q = polygon[polygon.size() - 1];
    for (int i = 0; i < polygon.size(); i++) {
        area += polygon[i][0] * q[1] - q[0] * polygon[i][1];
        q = polygon[i];
    }
    return std::abs(area / 2);
}

//vector<vector<vector<float>>> matchfun(vector<vector<vector<float>>> parts, vector<vector<int>> rotate) {
//    vector<vector<vector<float>>> square;
//    vector<vector<vector<float>>> notsquare;
//    vector<vector<int>> nsrotate;
//    vector<vector<int>> srotate;
//    vector<vector<vector<int>>> resultrotate;
//    vector<vector<vector<float>>> ret;
//    vector<vector<float>> combine;
//    vector<vector<vector<float>>> resultsquare;
//    vector<float> bestret = {0, 0, 0};
//    vector<float> ret_ = {0, 0, 0};
//    vector<int> bestrotate = {0, 0};
//    vector<int> rotate_ = {0, 0};
//    int squarenum = 0;
//    int notsquarenum = 0;
//    map<int, int> squareedge;
//    map<int, int> notsquareedge;
//    for (int i = 0; i < parts.size(); i++) {
////        if(parts[i].size()==4){
////            if(isrect(parts[i])){
////                square.push_back(parts[i]);
////            }
////            else if(polygonarea(parts[i])>squarearea(parts[i])*0.9){
////                square.push_back(parts[i]);
////            }
////        }
////        else
//        if (polygonarea(parts[i]) > squarearea(parts[i]) * 0.9) {
//            square.push_back(parts[i]);
//            srotate.push_back((rotate[i]));
//            squareedge[squarenum] = i;
//            squarenum++;
//        } else {
//            notsquare.push_back(parts[i]);
//            nsrotate.push_back((rotate[i]));
//            notsquareedge[notsquarenum] = i;
//            notsquarenum++;
//        }
//    }
//
//    for (int i = 0; i < notsquarenum; i++) {
//        resultrotate.emplace_back();
//        for (int j = 0; j < notsquarenum; j++) {
//            resultrotate[i].push_back({0, 0});
//        }
//    }
//    for (int i = 0; i < notsquarenum; i++) {
//        ret.emplace_back();
//        for (int j = 0; j < notsquarenum; j++) {
//            ret[i].push_back({0, 0});
//        }
//    }
//    Graph G(2 * notsquarenum);
//    vector<double> cost;
//    for (int i = 0; i < notsquarenum; i++) {
//        for (int j = i + 1; j < notsquarenum; j++) {
//            ret[i][j] = group_two(notsquare[i], notsquare[j], nsrotate[i], nsrotate[j], resultrotate[i][j]);
//            if (ret[i][j][2] > 0) {
//                G.AddEdge(i, j);
//                cost.push_back(-ret[i][j][2]);
//                G.AddEdge(i + notsquarenum, j + notsquarenum);
//                cost.push_back(-ret[i][j][2]);
//            }
//        }
//        G.AddEdge(i, i + notsquarenum);
//        cost.push_back(0);
//    }
//    bool flag;
//    Matching M(G);
//    pair<list<int>, double> solution = M.SolveMinimumCostPerfectMatching(cost);
//    list<int> matching(solution.first);
//    map<int, int> edge;
//    for (list<int>::iterator it = matching.begin(); it != matching.end(); it++) {
//        pair<int, int> e = G.GetEdge(*it);
//        if (e.first < notsquarenum and e.second < notsquarenum) {
//            int min_index = min(e.first, e.second);
//            int max_index = max(e.first, e.second);
//            edge[min_index] = max_index;
//            edge[max_index] = min_index;
//        }
//    }
//    for (map<int, int>::iterator it = edge.begin(); it != edge.end(); it++) {
//        if (it->second > it->first) {
//            combine.push_back({float(notsquareedge[it->first]), 1, 0, 0, float(resultrotate[it->first][it->second][0])});
//            combine.push_back({float(notsquareedge[it->second]), 1, ret[it->first][it->second][0], ret[it->first][it->second][1], float(resultrotate[it->first][it->second][1])});
//            resultsquare.push_back(combine);
//            combine.clear();
//        }
//    }
//    int index;
//    vector<int> usedsquare;
//    for (int i = 0; i < notsquare.size(); i++) {
//        index = -1;
//        flag = false;
//        for (map<int, int>::iterator it = edge.begin(); it != edge.end(); it++) {
//            if (i == it->first) {
//                flag = true;
//                continue;
//            }
//        }
//        if (flag == 1) {
//            continue;
//        }
//        for (int j = 0; j < square.size(); j++) {
//            flag = false;
//            for (int k = 0; k < usedsquare.size(); k++) {
//                if (j == usedsquare[k]) {
//                    flag = true;
//                    continue;
//                }
//
//            }
//            if (flag == true) {
//                continue;
//            }
//            bestret = {0, 0, 0};
//            if (squarearea(square[j]) <= squarearea(notsquare[i]) - polygonarea(notsquare[i])) {
//                ret_ = group_two(notsquare[i], square[j], nsrotate[i], srotate[j], rotate_);
//                if (bestret[2] < ret_[2]) {
//                    index = j;
//                    bestret = ret_;
//                    bestrotate = rotate_;
//                }
//            }
//        }
//        if (index != -1) {
//            combine.push_back({float(notsquareedge[i]), 1, 0, 0, float(bestrotate[0])});
//            combine.push_back({float(squareedge[index]), 0, bestret[0], bestret[1], float(bestrotate[1])});
//            resultsquare.push_back(combine);
//            usedsquare.push_back(index);
//            combine.clear();
//        } else {
//            combine.push_back({float(notsquareedge[i]), 1, 0, 0, 0});
//            resultsquare.push_back(combine);
//            combine.clear();
//        }
//    }
//    for (int j = 0; j < square.size(); j++) {
//        flag = false;
//        for (int i = 0; i < usedsquare.size(); i++) {
//            if (j == usedsquare[i]) {
//                flag = true;
//                continue;
//            }
//        }
//        if (not flag) {
//            combine.push_back({float(squareedge[j]), 0, 0, 0, 0});
//            resultsquare.push_back(combine);
//            combine.clear();
//        }
//    }
//    return resultsquare;
//}

#endif //SRC_UTILS_H