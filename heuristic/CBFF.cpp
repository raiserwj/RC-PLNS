#include "CBFF.h"
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <chrono>
using std::numeric_limits;
using std::max;
using std::min;
typedef boost::geometry::model::d2::point_xy<int> geopoint;
typedef boost::geometry::model::multi_point<geopoint> multiPoint;

vector<Part_Ptr_V> maxRects(Part_Ptr_V parts, float width, float height, int binNum, std::mt19937 &mt, float *score1, float *score2, bool junyun) {
    std::shuffle(parts.begin(), parts.end(), mt);

    vector<Part_Ptr_V> ret;
    vector<rbp::MaxRectsBinPack> detailpack = {};
    int bestRectIndex = -1;
    double bestScore1 = numeric_limits<double>::max();
    bool big_flag = true;
    int big_length = 1;

    for (int i = 0; i < max(binNum - 1, 1); i++) {
        detailpack.push_back(MaxRectsBinPack());
        detailpack.back().Init(width, height);
        ret.push_back(Part_Ptr_V());
        if (!junyun) continue;
        if (big_flag) {
            bestRectIndex = -1;
            bestScore1 = numeric_limits<double>::max();
            for (size_t i1 = 0; i1 < parts.size(); ++i1) {
                Part_Ptr current_part = parts[i1];
                if (parts[i1]->smallitem) {
                    continue;
                }
                double score1;
                score1 = detailpack[i].Insertscore(current_part->width, current_part->height, current_part->rotate,
                                                   detailpack[i].RectBestShortSideFit);
                if (score1 < bestScore1) {
                    bestScore1 = score1;
                    bestRectIndex = i1;
                }
            }

            if (bestRectIndex != -1) {
                Part_Ptr current_part = parts[bestRectIndex];
                detailpack[i].Insert(current_part->width, current_part->height, current_part->rotate,
                                     detailpack[i].RectBestShortSideFit);
                parts.erase(parts.begin() + bestRectIndex);
                ret.back().push_back(current_part);
            } else {
                big_flag = false;
            }
        }
    }
    while (parts.size() > 0) {
        bestScore1 = numeric_limits<double>::max();
        vector<vector<int>> indexs;
        for (int i1 = 0; i1 < parts.size(); ++i1) {
            double score1;
            Part_Ptr current_part = parts[i1];
            for (int j = 0; j < big_length; j++) {
                score1 = detailpack[j].Insertscore(current_part->width, current_part->height, current_part->rotate,
                                                   detailpack[0].RectBestShortSideFit);
                if (score1 < bestScore1) {
                    bestScore1 = score1;
                    indexs.clear();
                    indexs.push_back({j, i1});
                } else if (score1 == bestScore1 && score1 != numeric_limits<double>::max())indexs.push_back({j, i1});
            }
        }
        if (bestScore1 == numeric_limits<double>::max()) {
            if (big_length < binNum - 1) {
                big_length += 1;
            } else if (big_length == binNum - 1) {
                big_length += 1;

                detailpack.push_back(MaxRectsBinPack());
                detailpack.back().Init(width, height);
                ret.push_back(Part_Ptr_V());
                bestScore1 = numeric_limits<double>::max();
                bestRectIndex = -1;
                for (size_t i1 = 0; i1 < parts.size(); ++i1) {
                    Part_Ptr current_part = parts[i1];
                    if (big_flag && current_part->smallitem) {
                        continue;
                    }
                    double score1;
                    score1 = detailpack.back().Insertscore(current_part->width, current_part->height,
                                                           current_part->rotate,
                                                           detailpack.back().RectBestShortSideFit);
                    if (score1 < bestScore1) {
                        bestScore1 = score1;
                        bestRectIndex = i1;
                    }
                }
                if (bestRectIndex != -1) {
                    detailpack.back().Insert(parts[bestRectIndex]->width, parts[bestRectIndex]->height,
                                             parts[bestRectIndex]->rotate,
                                             detailpack.back().RectBestShortSideFit);
                    ret.back().push_back(parts[bestRectIndex]);
                    parts.erase(parts.begin() + bestRectIndex);
                }
            } else {
                break;
            }
        } else {
            std::uniform_int_distribution<> dis(0, indexs.size() - 1);
            int bestIndex = dis(mt);
            int bestRectIndex_ = indexs[bestIndex][1];
            int bestBinIndex_ = indexs[bestIndex][0];

            Part_Ptr current_part = parts[bestRectIndex_];
            detailpack[bestBinIndex_].Insert(current_part->width, current_part->height, current_part->rotate,
                                             rbp::MaxRectsBinPack::RectBestShortSideFit);
            ret[bestBinIndex_].push_back(current_part);
            parts.erase(parts.begin() + bestRectIndex_);
        }
    }
    *score1 = 0;
    for (auto part: parts) *score1 += part->size;
    *score2 = 0;
//    ret.erase(ret.begin() + detailpack.size(), ret.end());
    for (auto part: ret.back()) *score2 += part->size;
    ret.push_back(parts);
    return ret;
}
vector<Part_Ptr_V> maxRectBRKGA(Part_Ptr_V parts, float width, float height,  std::mt19937 &mt) {

    vector<Part_Ptr_V> ret;
    vector<rbp::MaxRectsBinPack> detailpack = {};
    int bestRectIndex = -1;
    double bestScore1 = numeric_limits<double>::max();
    bool big_flag = true;
    int big_length = 0;
    while (parts.size() > 0) {
        bestScore1 = numeric_limits<double>::max();
        vector<int> indexs;
        double score1;
        Part_Ptr current_part = parts[0];
        for (int j = 0; j < big_length; j++) {
            score1 = detailpack[j].Insertscore(current_part->width, current_part->height, current_part->rotate,
                                               detailpack[0].RectBestShortSideFit);
            if (score1 < bestScore1) {
                bestScore1 = score1;
                indexs.clear();
                indexs.push_back({j});
            } else if (score1 == bestScore1 && score1 != numeric_limits<double>::max()){
                indexs.clear();
                indexs.push_back({j});}
        }
        if (bestScore1 == numeric_limits<double>::max()) {
                big_length += 1;
                detailpack.push_back(MaxRectsBinPack());
                detailpack.back().Init(width, height);
                ret.push_back(Part_Ptr_V());
                bestScore1 = numeric_limits<double>::max();
                bestRectIndex = -1;
                    Part_Ptr current_part = parts[0];
                    if (big_flag && current_part->smallitem) {
                        continue;
                    }
                    double score1;
                    score1 = detailpack.back().Insertscore(current_part->width, current_part->height,
                                                           current_part->rotate,
                                                           detailpack.back().RectBestShortSideFit);
                    if (score1 < bestScore1) {
                        bestScore1 = score1;
                    }
                if (bestRectIndex != -1) {
                    detailpack.back().Insert(parts[bestRectIndex]->width, parts[bestRectIndex]->height,
                                             parts[bestRectIndex]->rotate,
                                             detailpack.back().RectBestShortSideFit);
                    ret.back().push_back(parts[bestRectIndex]);
                    parts.erase(parts.begin() + bestRectIndex);
                }
        } else {
            std::uniform_int_distribution<> dis(0, indexs.size() - 1);
            int bestIndex = dis(mt);
            int bestBinIndex_ = indexs[bestIndex];

            Part_Ptr current_part = parts[0];
            detailpack[bestBinIndex_].Insert(current_part->width, current_part->height, current_part->rotate,
                                             rbp::MaxRectsBinPack::RectBestShortSideFit);
            ret[bestBinIndex_].push_back(current_part);
            parts.erase(parts.begin());
        }
    }
    return ret;
}

bool cmpArea(Part_Ptr a, Part_Ptr b) {
    return a->size > b->size;
}

bool cmpPerimeter(Part_Ptr a, Part_Ptr b) {
    return (a->width + a->height) > (b->width + b->height);
}

bool cmpLongSide(Part_Ptr a, Part_Ptr b) {
    return max(a->width, a->height) > max(b->width, b->height);
}

bool cmpShortSide(Part_Ptr a, Part_Ptr b) {
    return min(a->width, a->height) > min(b->width, b->height);
}

vector<Part_Ptr_V> maxRectsSorts(Part_Ptr_V parts, float width, float height, int binNum, std::mt19937 &mt, float *score1, float *score2) {
    vector<Part_Ptr_V> bestBinRet;
    float bestScore1 = 10000000000, bestScore2 = 10000000000;

    for (int i = 0; i < 6; i++) {
        using namespace std::chrono;

        auto t0 = steady_clock::now();
        //面积，周长，长边，短边，随机1，随机2
        switch (i) {
            case 0:
                std::sort(parts.begin(), parts.end(), cmpArea);
                break;
            case 1:
                std::sort(parts.begin(), parts.end(), cmpPerimeter);
                break;
            case 2:
                std::sort(parts.begin(), parts.end(), cmpLongSide);
                break;
            case 3:
                std::sort(parts.begin(), parts.end(), cmpShortSide);
                break;
            default:
                std::shuffle(parts.begin(), parts.end(), mt);
                break;
        }
        auto t1 = steady_clock::now();

        // 毫秒
        double ms = duration<double, std::milli>(t1 - t0).count();
        // 或者用双精度毫秒

//        std::cout << ms << " ms\n";
        vector<rbp::MaxRectsBinPack> detailpack = {};
        vector<Part_Ptr_V> tempRet;
        detailpack.push_back(MaxRectsBinPack());
        detailpack.back().Init(width, height);
        tempRet.push_back(Part_Ptr_V());

        int j;
        for (j = 0; j < parts.size(); j++) {
            bool success = false;
            for (int k = 0; k < detailpack.size(); k++) {
                Rect rect = detailpack[k].Insert(parts[j]->width, parts[j]->height, parts[j]->rotate,
                                                 MaxRectsBinPack::RectBestShortSideFit);
                if (rect.height != 0) {
                    tempRet[k].push_back(parts[j]);
                    success = true;
                    break;
                }
            }
            if (!success && detailpack.size() >= binNum) {
                break;
            }
            if (!success && detailpack.size() < binNum) {
                detailpack.push_back(MaxRectsBinPack());
                detailpack.back().Init(width, height);
                tempRet.push_back(Part_Ptr_V());
                Rect rect = detailpack.back().Insert(parts[j]->width, parts[j]->height, parts[j]->rotate,
                                                     MaxRectsBinPack::RectBestShortSideFit);
                if (rect.height == 0) std::cout << "Error:CBFF" << std::endl;
                else tempRet.back().push_back(parts[j]);
            }
        }
        float tempScore1 = 0;
        for (int k = j; k < parts.size(); k++) tempScore1 += parts[k]->size;
        float tempScore2 = 0;
        for (auto part: tempRet.back()) tempScore2 += part->size;
        if (tempScore1 < bestScore1 || (tempScore1 == bestScore1 && tempScore2 < bestScore2)) {
            bestScore1 = tempScore1;
            bestScore2 = tempScore2;

            bestBinRet = tempRet;
            Part_Ptr_V tempParts(parts.begin() + j, parts.end());
            bestBinRet.push_back(tempParts);
        }
        auto t2 = steady_clock::now();

        // 毫秒
        double ms2 = duration<double, std::milli>(t2 - t0).count();
        // 或者用双精度毫秒

//        std::cout << ms2 << " ms\n";
    }
    *score1 = bestScore1;
    *score2 = bestScore2;
    return bestBinRet;
}

Part_Ptr_V insert(Bin_Ptr bin, Part_Ptr part, std::mt19937 &mt, bool initial_flag, bool return_flag, bool sort) {
    Part_Ptr_V bestBinParts;
    Part_Ptr_V bestRemainParts;
    Part_Ptr_V parts(bin->parts);
    parts.push_back(part);

    float score1MaxRects, score1Sort, score2MaxRects, score2Sort;
    vector<Part_Ptr_V> retMaxRects = maxRects(parts, bin->width, bin->height, 1, mt, &score1MaxRects, &score2MaxRects, false);
    vector<Part_Ptr_V> retSort;

    if (sort) {
        retSort = maxRectsSorts(parts, bin->width, bin->height, 1, mt, &score1Sort, &score2Sort);
    } else {
        score1Sort = numeric_limits<float>::max();
    }

    //成功的话，修改掉bin的parts，没成功的话还是用原来的parts
    if (initial_flag) {
        if (score1MaxRects > 0) {
            return {part};
        } else if (score1MaxRects == 0) {
            bin->parts = retMaxRects[0];
            return retMaxRects[1];
        } else {
            std::cout << "error" << std::endl;
            exit(-2);
        }
    } else {
        bool success = (score1MaxRects < part->size) | (score1Sort < part->size);
        if (!success) {
            return {part};
        } else if (score1MaxRects <= score1Sort) {
            bin->parts = retMaxRects[0];
            return retMaxRects[1];
        } else if (score1MaxRects > score1Sort) {
            bin->parts = retSort[0];
            return retSort[1];
        }
    }
}


Bin_Ptr_V regroup(Bin_Ptr_V bins, bool *success_flag, std::mt19937 &mt, bool return_flag, bool junyun, bool sort) {
    Bin_Ptr_V bins_ret = {};
    Part_Ptr_V parts;
    Bin_Ptr_V::iterator it;
    for (it = bins.begin(); it != bins.end(); it++) {
        parts.insert(parts.end(), (*it)->parts.begin(), (*it)->parts.end());
    }

    float areaSum = 0;
    for (int i = 0; i < parts.size(); i++) {
        areaSum += parts[i]->size;
    }

    float score1MaxRects, score1Sort, score2MaxRects, score2Sort;
    vector<Part_Ptr_V> retMaxRects = maxRects(parts, bins[0]->width, bins[0]->height, bins.size(), mt, &score1MaxRects, &score2MaxRects, junyun);
    bool successMaxRects = (score1MaxRects == 0);
    vector<Part_Ptr_V> retSort;
    bool successSort;
    if (sort) {
        retSort = maxRectsSorts(parts, bins[0]->width, bins[0]->height, bins.size(), mt, &score1Sort, &score2Sort);
        successSort = (score1Sort == 0);
    } else {
        successSort = false;
    }
    if (successMaxRects | successSort) {
        bool chooseMaxRects;
        if (successMaxRects && !successSort)chooseMaxRects = true;
        else if (!successMaxRects && successSort)chooseMaxRects = false;
        else {
            if (retMaxRects.size() < retSort.size()) chooseMaxRects = true;
            else if (retMaxRects.size() > retSort.size()) chooseMaxRects = false;
            else if (score2MaxRects < score2Sort) chooseMaxRects = true;
            else chooseMaxRects = false;
        }

        if (chooseMaxRects) {
            for (int i = 0; i < retMaxRects.size() - 1; i++) {
                if (retMaxRects[i].size() > 0) {
                    shared_ptr<amopt_pack::Bin> p_bin(
                            new amopt_pack::Bin(bins[0]->width, bins[0]->height, bins[i]->backornot, retMaxRects[i]));
                    bins_ret.push_back(p_bin);
                }
            }
        } else {
            for (int i = 0; i < retSort.size() - 1; i++) {
                if (retSort[i].size() > 0) {
                    shared_ptr<amopt_pack::Bin> p_bin(
                            new amopt_pack::Bin(bins[0]->width, bins[0]->height, bins[i]->backornot, retSort[i]));
                    bins_ret.push_back(p_bin);
                }
            }
        }
        *success_flag = true;
        return bins_ret;
    } else {
        *success_flag = false;
        return {};
    }
}

bool output_pack(Bin_Ptr bin) {
    rbp::MaxRectsBinPack maxrects = MaxRectsBinPack();
    maxrects.Init(bin->width, bin->height);
    for (int i = 0; i < bin->parts.size(); i++) {
        Rect rect = maxrects.Insert(bin->parts[i]->width, bin->parts[i]->height,
                                    bin->parts[i]->rotate,
                                    maxrects.RectBestShortSideFit);
        if (rect.height == 0) {
            std::cout << "false" << std::endl;
            exit(-1);
        }
        bin->parts[i]->rect = rect;
    }
    maxrects.free();
    return true;
}

bool output_pack_last(Bin_Ptr bin, int index, bool return_flag) {
    std::cout << "jinru" << std::endl;
    bool flag = false;
    if (index == 1) {
        int lb = 0;
        int ub = bin->height;
        while (ub - lb > 1) {
            int current = int((lb + ub) / 2);
            std::cout << lb << " " << current << "  " << ub << std::endl;
            vector<Part_Ptr> parts = bin->parts;
            rbp::MaxRectsBinPack maxrects = MaxRectsBinPack();
            maxrects.Init(bin->width, current);
            bool flag = true;
            while (parts.size() > 0 && flag) {
                double bestScore = numeric_limits<double>::max();
                double bestIndex = -1;
                for (int i = 0; i < parts.size(); i++) {
                    Part_Ptr current_part = parts[i];
                    double score = maxrects.Insertscore(current_part->width, current_part->height,
                                                        current_part->rotate,
                                                        maxrects.RectBestShortSideFit);
                    if (score < bestScore) {
                        bestScore = score;
                        bestIndex = i;
                    }
                }
                if (bestIndex != -1) {
                    parts[bestIndex]->temp_rect = maxrects.Insert(parts[bestIndex]->width, parts[bestIndex]->height,
                                                                  parts[bestIndex]->rotate,
                                                                  maxrects.RectBestShortSideFit);
                    parts.erase(parts.begin() + bestIndex);
                } else {
                    flag = false;
                }
            }
            if (parts.size() == 0) {
                ub = current;
                for (int i = 0; i < bin->parts.size(); i++) {
                    bin->parts[i]->rect = bin->parts[i]->temp_rect;
                    flag = true;
                }
            } else {
                lb = current + 1;
            }
            maxrects.free();
        }
        return flag;
    }
    if (index == 2) {
        int lb = 0;
        int ub = bin->width;
        while (ub - lb > 1) {
            int current = int((lb + ub) / 2);
            std::cout << lb << " " << current << "  " << ub << std::endl;
            vector<Part_Ptr> parts = bin->parts;
            rbp::MaxRectsBinPack maxrects = MaxRectsBinPack();
            maxrects.Init(current, bin->height);
            bool flag = true;
            while (parts.size() > 0 && flag) {
                double bestScore = numeric_limits<double>::max();
                double bestIndex = -1;
                for (int i = 0; i < parts.size(); i++) {
                    Part_Ptr current_part = parts[i];
                    double score = maxrects.Insertscore(current_part->width, current_part->height,
                                                        current_part->rotate,
                                                        maxrects.RectBestShortSideFit);
                    if (score < bestScore) {
                        bestScore = score;
                        bestIndex = i;
                    }
                }
                if (bestIndex != -1) {
                    parts[bestIndex]->temp_rect = maxrects.Insert(parts[bestIndex]->width, parts[bestIndex]->height,
                                                                  parts[bestIndex]->rotate,
                                                                  maxrects.RectBestShortSideFit);
                    parts.erase(parts.begin() + bestIndex);
                } else {
                    flag = false;
                }
            }
            if (parts.size() == 0) {
                ub = current;
                flag = true;
                for (int i = 0; i < bin->parts.size(); i++) {
                    bin->parts[i]->rect = bin->parts[i]->temp_rect;
                }
            } else {
                lb = current + 1;
            }
            maxrects.free();
        }
        return flag;
    }

}

bool output_pack_last_new(Bin_Ptr bin, int index, bool return_flag) {
    std::cout << "jinru" << std::endl;
    bool flag = false;
    if (index == 1) {
        int lb = 0;
        int ub = bin->height;
        double score = 9999;
        double maxhw;
        bool flag_ = false;
        multiPoint points;
        while (true) {
            int current = int((lb + ub) / 2);
            if (ub - lb <= 1) {
                current = ub;
                for (int n = 0; n < 500; n++) {
                    vector<Part_Ptr> parts = bin->parts;
                    rbp::MaxRectsBinPack maxrects = MaxRectsBinPack();
                    std::vector<int> sort;
                    for (int i = 0; i < parts.size(); i++) {
                        sort.push_back(i);
                    };
                    std::shuffle(sort.begin(), sort.end(), std::mt19937(std::random_device()()));
                    int sortflag = 0;
                    maxrects.Init(bin->width, current);
                    while (parts.size() > sortflag) {
                        int RectIndex = sort[sortflag];
                        double score1;

                        score1 = maxrects.Insertscore(parts[RectIndex]->width, parts[RectIndex]->height, parts[RectIndex]->rotate,
                                                      maxrects.RectBestShortSideFit);

                        if (score1 == numeric_limits<double>::max()) {
                            break;
                        } else {
                            parts[RectIndex]->temp_rect = maxrects.Insert(parts[RectIndex]->width, parts[RectIndex]->height,
                                                                          parts[RectIndex]->rotate,
                                                                          maxrects.RectBestShortSideFit);
                            sortflag += 1;
                        }
                    }
                    if (parts.size() == sortflag) {
                        maxhw = 0;
                        for (int i = 0; i < bin->parts.size(); i++) {
                            if (bin->parts[i]->temp_rect.y + bin->parts[i]->temp_rect.height +
                                bin->parts[i]->temp_rect.x + bin->parts[i]->temp_rect.width > maxhw) {
                                maxhw = bin->parts[i]->temp_rect.y + bin->parts[i]->temp_rect.height +
                                        bin->parts[i]->temp_rect.x + bin->parts[i]->temp_rect.width;
                            }
                        }
                        if (score > maxhw) {
                            score = maxhw;
                            for (int i = 0; i < bin->parts.size(); i++) {
                                bin->parts[i]->rect = bin->parts[i]->temp_rect;
                                flag = true;
                                flag_ = true;
                            }
                        }
                    }
                }
                return flag_;
            }

            std::cout << lb << " " << current << "  " << ub << std::endl;
            vector<Part_Ptr> parts = bin->parts;
            rbp::MaxRectsBinPack maxrects = MaxRectsBinPack();
            maxrects.Init(bin->width, current);
            bool flag = true;
            while (parts.size() > 0 && flag) {
                double bestScore = numeric_limits<double>::max();
                double bestIndex = -1;
                for (int i = 0; i < parts.size(); i++) {
                    Part_Ptr current_part = parts[i];
                    double score = maxrects.Insertscore(current_part->width, current_part->height,
                                                        current_part->rotate,
                                                        maxrects.RectBestShortSideFit);
                    if (score < bestScore) {
                        bestScore = score;
                        bestIndex = i;
                    }
                }
                if (bestIndex != -1) {
                    parts[bestIndex]->temp_rect = maxrects.Insert(parts[bestIndex]->width, parts[bestIndex]->height,
                                                                  parts[bestIndex]->rotate,
                                                                  maxrects.RectBestShortSideFit);
                    parts.erase(parts.begin() + bestIndex);
                } else {
                    flag = false;
                }
            }
            if (parts.size() == 0) {
                ub = current;
                for (int i = 0; i < bin->parts.size(); i++) {
                    bin->parts[i]->rect = bin->parts[i]->temp_rect;
                    flag = true;
                    flag_ = true;
                }
            } else {
                lb = current + 1;
            }
            maxrects.free();
        }
        return flag_;
    }
    if (index == 2) {
        int lb = 0;
        int ub = bin->width;
        double score = 9999;
        double maxhw;
        bool flag_ = false;
        multiPoint points;
        while (true) {
            int current = int((lb + ub) / 2);
            if (ub - lb <= 1) {
                for (int n = 0; n < 500; n++) {
                    vector<Part_Ptr> parts = bin->parts;
                    rbp::MaxRectsBinPack maxrects = MaxRectsBinPack();
                    std::vector<int> sort;
                    for (int i = 0; i < parts.size(); i++) {
                        sort.push_back(i);
                    };
                    std::shuffle(sort.begin(), sort.end(), std::mt19937(std::random_device()()));
                    int sortflag = 0;
                    maxrects.Init(current, bin->height);
                    while (parts.size() > sortflag) {
                        int RectIndex = sort[sortflag];
                        double score1;

                        score1 = maxrects.Insertscore(parts[RectIndex]->width, parts[RectIndex]->height, parts[RectIndex]->rotate,
                                                      maxrects.RectBestShortSideFit);

                        if (score1 == numeric_limits<double>::max()) {
                            break;
                        } else {
                            parts[RectIndex]->temp_rect = maxrects.Insert(parts[RectIndex]->width, parts[RectIndex]->height,
                                                                          parts[RectIndex]->rotate,
                                                                          maxrects.RectBestShortSideFit);
                            sortflag += 1;
                        }
                    }
                    if (parts.size() == sortflag) {

                        maxhw = 0;
                        for (int i = 0; i < bin->parts.size(); i++) {
                            if (bin->parts[i]->temp_rect.y + bin->parts[i]->temp_rect.height +
                                bin->parts[i]->temp_rect.x + bin->parts[i]->temp_rect.width > maxhw) {
                                maxhw = bin->parts[i]->temp_rect.y + bin->parts[i]->temp_rect.height +
                                        bin->parts[i]->temp_rect.x + bin->parts[i]->temp_rect.width;
                            }
                        }
                        if (score > maxhw) {
                            score = maxhw;
                            for (int i = 0; i < bin->parts.size(); i++) {
                                bin->parts[i]->rect = bin->parts[i]->temp_rect;
                                flag = true;
                                flag_ = true;
                            }
                        }
                    }
                }
                return flag_;
            }

            std::cout << lb << " " << current << "  " << ub << std::endl;
            vector<Part_Ptr> parts = bin->parts;
            rbp::MaxRectsBinPack maxrects = MaxRectsBinPack();
            maxrects.Init(current, bin->height);
            bool flag = true;
            while (parts.size() > 0 && flag) {
                double bestScore = numeric_limits<double>::max();
                double bestIndex = -1;
                for (int i = 0; i < parts.size(); i++) {
                    Part_Ptr current_part = parts[i];
                    double score = maxrects.Insertscore(current_part->width, current_part->height,
                                                        current_part->rotate,
                                                        maxrects.RectBestShortSideFit);
                    if (score < bestScore) {
                        bestScore = score;
                        bestIndex = i;
                    }
                }
                if (bestIndex != -1) {
                    parts[bestIndex]->temp_rect = maxrects.Insert(parts[bestIndex]->width, parts[bestIndex]->height,
                                                                  parts[bestIndex]->rotate,
                                                                  maxrects.RectBestShortSideFit);
                    parts.erase(parts.begin() + bestIndex);
                } else {
                    flag = false;
                }
            }
            if (parts.size() == 0) {
                ub = current;
                flag = true;
                flag_ = true;
                for (int i = 0; i < bin->parts.size(); i++) {
                    bin->parts[i]->rect = bin->parts[i]->temp_rect;
                }
            } else {
                lb = current + 1;
            }
            maxrects.free();
        }
        return flag_;
    }
}

void last_insert(vector<Bin_Ptr> bins, int index, std::mt19937 &mt) {
    int i = 0;
    Part_Ptr_V tempParts(bins[index]->parts);
    while (i < tempParts.size()) {
        bool flag = false;
        for (int j = 0; j < bins.size(); j++) {
            if (j == index) continue;
            Part_Ptr_V allParts(bins[j]->parts);
            allParts.push_back(tempParts[i]);

            float score1MaxRects, score1Sort, score2MaxRects, score2Sort;
            vector<Part_Ptr_V> retMaxRects = maxRects(allParts, bins[j]->width, bins[j]->height, 1, mt, &score1MaxRects, &score2MaxRects, false);
            if (score1MaxRects == 0) {
                bins[j]->parts = retMaxRects[0];
                flag = true;
                break;
            }

            vector<Part_Ptr_V> retSort;
            retSort = maxRectsSorts(allParts, bins[j]->width, bins[j]->height, 1, mt, &score1Sort, &score2Sort);
            if (score1Sort == 0) {
                bins[j]->parts = retSort[0];
                flag = true;
                break;
            }
        }
        if (!flag) i += 1;
        else {
            tempParts.erase(tempParts.begin() + i);
        }
    }
    bins[index]->parts = tempParts;

//    for (int i = 0; i < bins.size(); i++) {
//        if (i == index) {
//            continue;
//        }
//        vector<Part_Ptr> parts = bins[i]->parts;
//        int j = 0;
//        while (bins[index]->parts.size() > j) {
//            Part_Ptr part_ = bins[index]->parts[j];
//            rbp::MaxRectsBinPack maxrects = MaxRectsBinPack();
//            maxrects.Init(bins[i]->width, bins[i]->height);
//            bool flag = true;
//            while (parts.size() > 0 && flag) {
//                double bestScore = numeric_limits<double>::max();
//                double bestIndex = -1;
//                for (int i = 0; i < parts.size(); i++) {
//                    Part_Ptr current_part = parts[i];
//                    double score = maxrects.Insertscore(current_part->width, current_part->height,
//                                                        current_part->rotate,
//                                                        maxrects.RectBestShortSideFit);
//                    if (score < bestScore) {
//                        bestScore = score;
//                        bestIndex = i;
//                    }
//                }
//                if (bestIndex != -1) {
//                    parts[bestIndex]->temp_rect = maxrects.Insert(parts[bestIndex]->width, parts[bestIndex]->height,
//                                                                  parts[bestIndex]->rotate,
//                                                                  maxrects.RectBestShortSideFit);
//                    parts.erase(parts.begin() + bestIndex);
//                } else {
//                    flag = false;
//                    j++;
//                    break;
//                }
//            }
//            if (flag) {
//                double score = maxrects.Insertscore(bins[index]->parts[j]->width, bins[index]->parts[j]->height,
//                                                    bins[index]->parts[j]->rotate,
//                                                    maxrects.RectBestShortSideFit);
//                if (score < std::numeric_limits<double>::max()) {
//                    part_->temp_rect = maxrects.Insert(part_->width, part_->height,
//                                                       part_->rotate,
//                                                       maxrects.RectBestShortSideFit);
//                    bins[i]->parts.push_back(part_);
//                    bins[i]->normal_insert_flag = true;
//                    bins[index]->parts.erase(bins[index]->parts.begin() + j);
//
//                } else { j++; }
//            }
//        }
//        for (j = 0; j < bins[i]->parts.size(); j++) {
//            bins[i]->parts[j]->rect = bins[i]->parts[j]->temp_rect;
//        }
//    }
}