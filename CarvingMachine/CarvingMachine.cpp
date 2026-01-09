//
// Created by raiser on 2021/12/13.
//

#include "CarvingMachine.h"
#include <chrono>
#include <sys/time.h>
#include<math.h>
#include <algorithm>
#include <limits>
#include<cassert>
using namespace std;
using Clock = std::chrono::steady_clock;
using microseconds = std::chrono::microseconds;
vector<vector<double>> Input(string str) {
    Json::CharReaderBuilder rbuilder;
    Json::CharReader *reader(rbuilder.newCharReader());
//    rbuilder["collectComments"] = false;
    Json::Value root_group;
    JSONCPP_STRING errs;
    int partnum=0;
    vector<vector<double>> parts={};
    reader->parse(str.data(), str.data() + str.size(), &root_group, &errs);
    delete reader;
    for (int i = 0; i < root_group["items"].size(); i++) {
        float x_min = 100000, x_max = -100000, y_min = 100000, y_max = -100000;
        std::vector<std::vector<float>> points;
        for (int j = 0; j < root_group["items"][i]["points"].size(); j++) {
            std::vector<double> point;
            point.push_back(root_group["items"][i]["points"][j][0].asFloat());
            point.push_back(root_group["items"][i]["points"][j][1].asFloat());
            if (point[0] < x_min) x_min = point[0];
            if (point[0] > x_max) x_max = point[0];
            if (point[1] < y_min) y_min = point[1];
            if (point[1] > y_max) y_max = point[1];
        }
        parts.push_back({0,0,x_max-x_min,y_max-y_min});
    }
    return parts;
}
vector<double> Inputbinsize(string str){
    Json::CharReaderBuilder rbuilder;
    Json::CharReader *reader(rbuilder.newCharReader());
//    rbuilder["collectComments"] = false;
    Json::Value root_group;
    JSONCPP_STRING errs;
    int partnum=0;
    vector<vector<double>> parts={};
    reader->parse(str.data(), str.data() + str.size(), &root_group, &errs);
    delete reader;
    return {root_group["plates"][0]["width"].asDouble(),root_group["plates"][0]["height"].asDouble()};
}
static inline double binFillArea(const Bin& b) {
    double s = 0.0;
    for (const auto& r : b.rects) {
        if (r.w == 0.0 || r.h == 0.0) continue; // 若无“无效矩形”，可删掉
        s += r.w * r.h;
    }
    return s;
}
std::vector<Rect> destroy_solution_1(std::vector<Bin> &sol,
                                     int binnum,
                                     std::mt19937 &mt)
{
    std::vector<Rect> destroy_set;

    const int n = (int)sol.size();
    if (n == 0 || binnum <= 0) return destroy_set;

    const double binArea = 1;
    if (binArea <= 0.0) return destroy_set;

    // 1) 计算每个 bin 的利用率
    std::vector<double> util(n, 0.0);
    for (int i = 0; i < n; ++i) {
        util[i] = binFillArea(sol[i]) / binArea;
    }

    // 2) 按利用率从高到低排序索引（高利用率在前，低利用率在后）
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](int a, int b) { return util[a] > util[b]; });

    // 3) 取“后 75%”：去掉前 25%（最高利用率的 25%）
    int cut =0;
    if(binnum>3){
         cut = (int)std::ceil(0.25 * n);
    }     // 前 25% 的数量
    else{
         cut = (int)std::ceil(0 * n);
    }
    std::vector<int> candidates;
    candidates.reserve(n - cut);
    for (int k = cut; k < n; ++k) candidates.push_back(idx[k]);

    // 候选不足时退化：用全部 bin
    if ((int)candidates.size() < 1) candidates = idx;

    // 4) 从 candidates 中随机选 num 个要删除的 bin
    const int num = std::min(binnum, (int)candidates.size());
    std::shuffle(candidates.begin(), candidates.end(), mt);
    candidates.resize(num);

    // 5) 标记删除 + 预留空间
    std::vector<char> selected(n, 0);
    size_t total_rects = 0;
    for (int id : candidates) {
        selected[id] = 1;
        total_rects += sol[id].rects.size();
    }
    destroy_set.reserve(total_rects);

    // 6) 重建 remaining
    std::vector<Bin> remaining;
    remaining.reserve(n - num);

    for (int i = 0; i < n; ++i) {
        if (selected[i]) {
            auto& rs = sol[i].rects;
            destroy_set.insert(destroy_set.end(),
                               std::make_move_iterator(rs.begin()),
                               std::make_move_iterator(rs.end()));
        } else {
            remaining.push_back(std::move(sol[i]));
        }
    }

    sol.swap(remaining);
    return destroy_set;
}
std::vector<Rect> destroy_solution_3(std::vector<Bin>& sol, int binnum) {
    std::vector<Rect> destroy_set;
    if (sol.empty() || binnum <= 0) return destroy_set;

    const size_t n = sol.size();
    const size_t k = std::min<size_t>((size_t)binnum, n);

    // 1) 计算每个 bin 的填充面积（sum(w*h)）
    std::vector<double> fillarea(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (const auto& r : sol[i].rects) {
            if (r.w == 0.0 || r.h == 0.0) continue; // 若你不需要无效判断可去掉
            s += r.w * r.h;
        }
        fillarea[i] = s;
    }

    // 2) 生成索引并按 fillarea 从小到大排序
    std::vector<int> idx((int)n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](int a, int b) { return fillarea[(size_t)a] < fillarea[(size_t)b]; });

    // 3) 取最小的 k 个 bin，下标从大到小删除，避免 erase 造成下标偏移
    std::vector<int> to_remove(idx.begin(), idx.begin() + (int)k);
    std::sort(to_remove.begin(), to_remove.end(), std::greater<int>());

    // （可选）预留 destroy_set 大概容量，减少扩容
    size_t approx = 0;
    for (int i : to_remove) approx += sol[(size_t)i].rects.size();
    destroy_set.reserve(approx);

    for (int i : to_remove) {
        // 把该 bin 的矩形加入 destroy_set（这里用 move，避免拷贝）
        auto& rects = sol[(size_t)i].rects;
        destroy_set.insert(destroy_set.end(),
                           std::make_move_iterator(rects.begin()),
                           std::make_move_iterator(rects.end()));

        // 删除该 bin
        sol.erase(sol.begin() + i);
    }

    return destroy_set;
}

std::vector<Rect> destroy_solution_2(std::vector<Bin>& sol) {
    if (sol.empty()) return {};

    int bestPos = -1;
    double bestAvg = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < sol.size(); ++i) {
        const auto& bin = sol[i];
        if (bin.rects.empty()) continue;  // 避免除零；也避免选到空 bin

        double sumArea = 0.0;
        for (const auto& r : bin.rects) {
            // 若你用 h==0 表示无效，可在这里跳过无效矩形
            if (r.w == 0.0 || r.h == 0.0) continue;
            sumArea += r.w * r.h;
        }

        // 注意：你原代码用 sol[i].size() 作为分母（包含可能的“空条目”）。
        // 如果你希望严格等价，就用 bin.rects.size()；如果你希望更合理，就用有效计数。
//        const double avg = sumArea / (double)bin.rects.size();
        const double avg = sumArea ;

        if (avg < bestAvg) {
            bestAvg = avg;
            bestPos = (int)i;
        }
    }

    if (bestPos == -1) {
        // 全是空 bin 的情况：按原逻辑会出错；这里选择不破坏 sol
        return {};
    }

    std::vector<Rect> destroy_set = std::move(sol[bestPos].rects);
    sol.erase(sol.begin() + bestPos);
    return destroy_set;
}
std::vector<Rect> destroy_solution_4(std::vector<Bin>& sol) {
    if (sol.empty()) return {};

    int bestPos1 = -1; // 最小 avg
    int bestPos2 = -1; // 第二小 avg
    double bestAvg1 = std::numeric_limits<double>::infinity();
    double bestAvg2 = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < sol.size(); ++i) {
        const auto& bin = sol[i];
        if (bin.rects.empty()) continue; // 避免除零/空 bin

        double sumArea = 0.0;
        for (const auto& r : bin.rects) {
            if (r.w == 0.0 || r.h == 0.0) continue; // 若你没有“无效矩形”，可删掉
            sumArea += r.w * r.h;
        }

//        const double avg = sumArea / (double)bin.rects.size();
        const double avg = sumArea;

        if (avg < bestAvg1) {
            bestAvg2 = bestAvg1; bestPos2 = bestPos1;
            bestAvg1 = avg;      bestPos1 = (int)i;
        } else if (avg < bestAvg2) {
            bestAvg2 = avg;      bestPos2 = (int)i;
        }
    }

    // 没有任何非空 bin
    if (bestPos1 == -1) return {};

    // 只有一个非空 bin：退化为最小
    int target = (bestPos2 != -1) ? bestPos2 : bestPos1;

    std::vector<Rect> destroy_set = std::move(sol[target].rects);
    sol.erase(sol.begin() + target);
    return destroy_set;
}
std::vector<Rect> destroy_solution_5(std::vector<Bin>& sol, std::mt19937& mt) {
    if (sol.empty()) return {};

    std::vector<int> nonEmpty;
    nonEmpty.reserve(sol.size());
    for (int i = 0; i < (int)sol.size(); ++i) {
        if (!sol[(size_t)i].rects.empty()) nonEmpty.push_back(i);
    }
    if (nonEmpty.empty()) return {};

    std::uniform_int_distribution<int> dis(0, (int)nonEmpty.size() - 1);
    const int pos = nonEmpty[(size_t)dis(mt)];

    std::vector<Rect> destroy_set = std::move(sol[(size_t)pos].rects);
    sol.erase(sol.begin() + pos);
    return destroy_set;
}
std::vector<Rect> destroy_solution_6(std::vector<Bin>& sol) {
    if (sol.empty()) return {};

    int bestPos1 = -1, bestPos2 = -1; // 最小、第二小
    double bestSum1 = std::numeric_limits<double>::infinity();
    double bestSum2 = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < sol.size(); ++i) {
        const auto& bin = sol[i];
        if (bin.rects.empty()) continue;  // 避免选到空 bin

        double sumArea = 0.0;
        for (const auto& r : bin.rects) {
            if (r.w == 0.0 || r.h == 0.0) continue; // 跳过无效矩形（如你有此约定）
            sumArea += r.w * r.h;
        }

        // 用总面积做比较：维护最小与第二小
        if (sumArea < bestSum1) {
            bestSum2 = bestSum1; bestPos2 = bestPos1;
            bestSum1 = sumArea;  bestPos1 = (int)i;
        } else if (sumArea < bestSum2) {
            bestSum2 = sumArea;  bestPos2 = (int)i;
        }
    }

    // 不足两个非空 bin：不破坏 sol
    if (bestPos2 == -1) return {};

    std::vector<Rect> destroy_set = std::move(sol[bestPos2].rects);
    sol.erase(sol.begin() + bestPos2);
    return destroy_set;
}
void repair_solution1(
        std::vector<Bin>& solution,
        const std::vector<Rect>& destroy_set,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
) {
    if (destroy_set.empty()) return;

    bool flag = false;
    std::vector<Bin> bins = regroup_1(destroy_set, &flag, binSize, destroy_num, mt);

    solution.insert(solution.end(),
                    std::make_move_iterator(bins.begin()),
                    std::make_move_iterator(bins.end()));
}
void repair_solution2(
        std::vector<Bin>& solution,
        const std::vector<Rect>& destroy_set,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
) {
    if (destroy_set.empty()) return;

    bool flag = false; // 原代码没用 flag，这里保留，便于你后续扩展
    std::vector<Bin> bins = regroup_2(destroy_set, &flag, binSize, destroy_num, mt);

    solution.insert(solution.end(),
                    std::make_move_iterator(bins.begin()),
                    std::make_move_iterator(bins.end()));
}
bool repair_solution3(
        std::vector<Bin>& solution,
        const std::vector<Rect>& destroy_set,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
) {
    if (destroy_set.empty()) return false;

    bool flag = false;
    std::vector<Bin> temp_solution = std::move(regroup_3(solution, destroy_set, &flag, binSize, destroy_num, mt));
    if(flag==true){
        solution=temp_solution;
    }
    return flag;
}
bool repair_solution5(
        std::vector<Bin>& solution,
        const std::vector<Rect>& destroy_set,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
) {
    if (destroy_set.empty()) return false;

    bool flag = false;
    std::vector<Bin> temp_solution = std::move(regroup_3_(solution, destroy_set, &flag, binSize, destroy_num, mt));
    if(flag==true){
        solution=temp_solution;
    }
    return flag;
}
bool repair_solution4(
        std::vector<Bin>& solution,
        const std::vector<Rect>& destroy_set,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
) {
    if (destroy_set.empty()) return false;

    bool flag = false;
    std::vector<Bin> newBins = regroup_4(destroy_set, &flag, binSize, destroy_num, mt);

    solution.insert(solution.end(),
                    std::make_move_iterator(newBins.begin()),
                    std::make_move_iterator(newBins.end()));

    // 你原来这里的 count_in_bins 未完成；如果需要统计，可在此补充
    return flag;
}

    // 统计 bins 中所有 vector<double> 的总数
//    for (const auto& bin : bins) {
//        count_in_bins += bin.size();
//    }

    // destroy_set 中的数量就是 destroy_set.size()
//    std::size_t count_in_destroy = destroy_set.size();
//    if (flag && count_in_bins == count_in_destroy){
//        if (sumarea(bins)!=sumarea({destroy_set})){
//            flag=false;
//        }
//        int binnum=0;
//        for (std::size_t i = 0; i < temp_bins->size(); ++i) {
//            if (!(*temp_bins)[i].empty()) {
//                ++binnum;
//            }
//        }
//        std::cout<<binnum-27;
//    }
//    else{
//        flag=false;
//    }

//bool repair_solution5(vector<vector<vector<double>>> &solution, vector<vector<double>> destroy_set,vector<double> binsize,int destroy_num,std::mt19937 mt) {
//    if (destroy_set.size() == 0) return false;
//    bool flag;
//    vector<vector<vector<double>>> *temp_bins;
//    temp_bins = &solution;
//    vector<vector<vector<double>>> bins = regroup_5(destroy_set, &flag, binsize, destroy_num,mt);
//
//    (*temp_bins).insert((*temp_bins).end(), bins.begin(),
//                        bins.end());
//    std::size_t count_in_bins = 0;
//
//    // 统计 bins 中所有 vector<double> 的总数
//    for (const auto& bin : bins) {
//        count_in_bins += bin.size();
//    }
//
//    // destroy_set 中的数量就是 destroy_set.size()
//    std::size_t count_in_destroy = destroy_set.size();
//    if (flag && count_in_bins == count_in_destroy){
//        if (sumarea(bins)!=sumarea({destroy_set})){
//            flag=false;
//        }
////        int binnum=0;
////        for (std::size_t i = 0; i < temp_bins->size(); ++i) {
////            if (!(*temp_bins)[i].empty()) {
////                ++binnum;
////            }
////        }
////        std::cout<<binnum-27;
//    }
//    else{
//        flag=false;
//    }
//
//    return flag;
//}
bool rectsOverlap(const vector<double>& a, const vector<double>& b) {
    double ax1 = a[0];
    double ay1 = a[1];
    double ax2 = a[0] + a[2];
    double ay2 = a[1] + a[3];

    double bx1 = b[0];
    double by1 = b[1];
    double bx2 = b[0] + b[2];
    double by2 = b[1] + b[3];

    // 不重叠的四种情况（<= 允许边界刚好贴住）
    if (ax2 <= bx1 + 1e-9) return false; // a 在 b 左侧
    if (bx2 <= ax1 + 1e-9) return false; // b 在 a 左侧
    if (ay2 <= by1 + 1e-9) return false; // a 在 b 下方
    if (by2 <= ay1 + 1e-9) return false; // b 在 a 下方

    // 其他情况都视为有重叠
    return true;
}

// 检查所有 bin 是否存在重叠，存在则返回 false，并打印出问题位置
bool check_bins_no_overlap(const vector<vector<vector<double>>>& bins) {
    bool ok = true;
    for (size_t b = 0; b < bins.size(); ++b) {
        const auto& bin = bins[b];
        for (size_t i = 0; i < bin.size(); ++i) {
            for (size_t j = i + 1; j < bin.size(); ++j) {
                if (rectsOverlap(bin[i], bin[j])) {
                    ok = false;
                }
            }
        }
    }
    return ok;
}
std::vector<Bin> regroup_1(
        std::vector<Rect> destroy_set,   // 按值：内部可直接传给 maxRects
        bool* success_flag,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
) {
    std::vector<Bin> bins_ret;
    if (destroy_set.empty()) {
        if (success_flag) *success_flag = false;
        return bins_ret;
    }

    float score1MaxRects = 0.0f, score2MaxRects = 0.0f;

    // 直接把 destroy_set 作为 parts 传入（maxRects 内部会 shuffle/erase）
    std::vector<Bin> retMaxRects = maxRects(
            std::move(destroy_set),
            binSize.w,
            binSize.h,
            destroy_num,
            &score1MaxRects,
            &score2MaxRects,
            mt
    );

    const bool successMaxRects = (score1MaxRects == 0.0f);

    if (!successMaxRects) {
        if (success_flag) *success_flag = false;
        return {}; // 与原逻辑一致
    }

    // 成功：只保留非空 bin
    bins_ret.reserve(retMaxRects.size());
    for (auto& b : retMaxRects) {
        if (!b.rects.empty()) {
            bins_ret.push_back(std::move(b)); // move：避免拷贝
        }
    }

    if (success_flag) *success_flag = true;
    return bins_ret;
}

std::vector<Bin> regroup_2(
        std::vector<Rect> destroy_set,   // 按值：内部要 sort
        bool* success_flag,
        const Rect& binSize,             // 只用 binSize.w / binSize.h
        int destroy_num,
        std::mt19937& mt                 // 关键：按引用
) {
    std::sort(destroy_set.begin(), destroy_set.end(), cmpArea);

    std::vector<Bin> bins;
    bins.reserve(std::min<size_t>(destroy_set.size(), (destroy_num > 0 ? (size_t)destroy_num : destroy_set.size())));
    bins.push_back(Bin{});

    // 每个 bin 一份 free-rectangles
    std::vector<std::vector<Rect>> freeRects;
    freeRects.reserve(bins.capacity());
    freeRects.push_back({ Rect{0, 0, binSize.w, binSize.h} });

    for (const auto& item : destroy_set) {
        bool placed = false;

        // 尝试放入已有 bin
        for (size_t k = 0; k < bins.size(); ++k) {
            Rect pos = Insert(item.w, item.h, freeRects[k], mt);
            if (pos.h != 0) {
                bins[k].rects.push_back(pos);
                placed = true;
                break;
            }
        }

        // 放不进去则新开一个 bin
        if (!placed) {
            bins.push_back(Bin{});
            freeRects.push_back({ Rect{0, 0, binSize.w, binSize.h} });

            Rect pos = Insert(item.w, item.h, freeRects.back(), mt);

            // 你原代码这里是 push destroy_set[j]（忽略了 pos）。更合理的做法是：
            if (pos.h != 0) bins.back().rects.push_back(pos);
            else            bins.back().rects.push_back(item); // 极端情况下 Insert 仍失败则保底
        }
    }

    if (success_flag) {
        // 你原函数没设置 success_flag；这里给一个常见判定（可按你 repair 的逻辑自行调整）
        *success_flag = (destroy_num <= 0) ? true : ((int)bins.size() <= destroy_num);
    }
    return bins;
}
std::vector<Bin> regroup_3(
        std::vector<Bin> solution,                 // 按值：返回新解（与原函数一致）
        const std::vector<Rect>& destroy_set,       // const&：避免拷贝
        bool* success_flag,
        const Rect& binSize,
        int /*destroy_num*/,                        // 原函数里没用到，保留占位
        std::mt19937& mt                            // 关键：按引用
) {
    std::vector<Bin> temp_bins = std::move(solution);

    if (success_flag) *success_flag = false;

    // 与原逻辑一致：空解则直接把 destroy_set 作为一个新 bin 塞进去
    if (temp_bins.empty()) {
        Bin b;
        b.rects = destroy_set;
        temp_bins.push_back(std::move(b));
        if (success_flag) *success_flag = true;
        return temp_bins;
    }

    // parts = destroy_set 并 shuffle（保留 destroy_set 原顺序用于失败分支）
    std::vector<Rect> parts = destroy_set;
    std::shuffle(parts.begin(), parts.end(), mt);

    std::vector<Rect> remain_parts;
    remain_parts.reserve(parts.size());

    // 用“数组 + 标记”替代 map：索引范围固定为 [0, temp_bins.size())
    const size_t nBins = temp_bins.size();
    std::vector<Bin> modified(nBins);       // 存放被修改过的 bin
    std::vector<char> touched(nBins, 0);    // 是否已拷贝到 modified
    std::vector<int> touched_idx;
    touched_idx.reserve(std::min(nBins, parts.size()));

    // 随机选择要尝试插入的 bin
    int minIdx = -1;
    double minA = std::numeric_limits<double>::infinity();
    for (int i = 0; i < (int)nBins; ++i) {
        const double a = sumArea(temp_bins[(size_t)i].rects); // 你已有的 sumArea(rects)
        if (a < minA) { minA = a; minIdx = i; }
    }

// 随机选 bin：排除 minIdx
    if (nBins <= 1 || minIdx < 0) {
        return {}; // 或者按你的逻辑处理：无法排除/无法选择
    }
    std::uniform_int_distribution<int> pickBin(0, (int)nBins - 2);
    int r = pickBin(mt);
    int idx = (r >= minIdx) ? (r + 1) : r;

    for (const auto& p : parts) {
        const int idx = pickBin(mt);
        bool temp_flag=false;
        modified[(size_t)idx] = temp_bins[(size_t)idx]; // 第一次命中才复制

        // 尝试插入；失败则返回 {p}
        std::vector<Rect> leftover = insertparts_1(modified[(size_t)idx].rects, p, binSize, mt,&temp_flag);
        if (temp_flag==false){
            remain_parts.insert(remain_parts.end(),
                                std::make_move_iterator(leftover.begin()),
                                std::make_move_iterator(leftover.end()));
        }
        else if (!leftover.empty()) {
            remain_parts.insert(remain_parts.end(),
                                std::make_move_iterator(leftover.begin()),
                                std::make_move_iterator(leftover.end()));
            temp_bins[(size_t)idx]=modified[(size_t)idx];
        }
        else{
            temp_bins[(size_t)idx]=modified[(size_t)idx];
        }
    }

    // 如果存在放不进去的 parts：尝试用 regroup_1 再打一个新 bin
    if (!remain_parts.empty()) {
        bool flag_pack = false;
        std::vector<Bin> new_bins = regroup_4(remain_parts, &flag_pack, binSize, 1, mt);

        if (flag_pack && !new_bins.empty() && !new_bins[0].rects.empty()) {

            // 写回被修改的 bins
            // 加入新 bin
            temp_bins.push_back(std::move(new_bins[0]));
            if (success_flag) *success_flag = true;
        } else {
            // 与原逻辑一致：失败则把原 destroy_set 整体作为一个新 bin 追加
            Bin b;
            b.rects = destroy_set;
            temp_bins.push_back(std::move(b));
            *success_flag = false;
        }
//        int sumnum=0;
//        for (int i=0;i<temp_bins.size();i++){
//            sumnum+=temp_bins[i].rects.size();
//        }
//        if (sumnum>100 && *success_flag==true){
//            std::cout<<1;
//        }
        return temp_bins;
    }

    // remain_parts 为空：仅写回被修改的 bins

    *success_flag = true;
    return temp_bins;
}
std::vector<Bin> regroup_3_(
        std::vector<Bin> solution,                 // 按值：返回新解（与原函数一致）
        const std::vector<Rect>& destroy_set,       // const&：避免拷贝
        bool* success_flag,
        const Rect& binSize,
        int /*destroy_num*/,                        // 原函数里没用到，保留占位
        std::mt19937& mt                            // 关键：按引用
) {
    std::vector<Bin> temp_bins = std::move(solution);

    if (success_flag) *success_flag = false;

    // 与原逻辑一致：空解则直接把 destroy_set 作为一个新 bin 塞进去
    if (temp_bins.empty()) {
        Bin b;
        b.rects = destroy_set;
        temp_bins.push_back(std::move(b));
        if (success_flag) *success_flag = true;
        return temp_bins;
    }

    // parts = destroy_set 并 shuffle（保留 destroy_set 原顺序用于失败分支）
    std::vector<Rect> parts = destroy_set;
    std::shuffle(parts.begin(), parts.end(), mt);

    std::vector<Rect> remain_parts;
    remain_parts.reserve(parts.size());

    // 用“数组 + 标记”替代 map：索引范围固定为 [0, temp_bins.size())
    const size_t nBins = temp_bins.size();
    std::vector<Bin> modified(nBins);       // 存放被修改过的 bin
    std::vector<char> touched(nBins, 0);    // 是否已拷贝到 modified
    std::vector<int> touched_idx;
    touched_idx.reserve(std::min(nBins, parts.size()));

    // 随机选择要尝试插入的 bin
    int minIdx = -1;
    double minA = std::numeric_limits<double>::infinity();
    for (int i = 0; i < (int)nBins; ++i) {
        const double a = sumArea(temp_bins[(size_t)i].rects); // 你已有的 sumArea(rects)
        if (a < minA) { minA = a; minIdx = i; }
    }

// 随机选 bin：排除 minIdx
    if (nBins <= 1 || minIdx < 0) {
        return {}; // 或者按你的逻辑处理：无法排除/无法选择
    }
    std::uniform_int_distribution<int> pickBin(0, (int)nBins - 1);
    int r = pickBin(mt);
//    int idx = (r >= minIdx) ? (r + 1) : r;

    for (const auto& p : parts) {
        const int idx = pickBin(mt);
        bool temp_flag=false;
        modified[(size_t)idx] = temp_bins[(size_t)idx]; // 第一次命中才复制

        // 尝试插入；失败则返回 {p}
        std::vector<Rect> leftover = insertparts_1(modified[(size_t)idx].rects, p, binSize, mt,&temp_flag);
        if (temp_flag==false){
            remain_parts.insert(remain_parts.end(),
                                std::make_move_iterator(leftover.begin()),
                                std::make_move_iterator(leftover.end()));
        }
        else if (!leftover.empty()) {
            remain_parts.insert(remain_parts.end(),
                                std::make_move_iterator(leftover.begin()),
                                std::make_move_iterator(leftover.end()));
            temp_bins[(size_t)idx]=modified[(size_t)idx];
        }
        else{
            temp_bins[(size_t)idx]=modified[(size_t)idx];
        }
    }

    // 如果存在放不进去的 parts：尝试用 regroup_1 再打一个新 bin
    if (!remain_parts.empty()) {
        bool flag_pack = false;
        std::vector<Bin> new_bins = regroup_4(remain_parts, &flag_pack, binSize, 1, mt);

        if (flag_pack && !new_bins.empty() && !new_bins[0].rects.empty()) {

            // 写回被修改的 bins
            // 加入新 bin
            temp_bins.push_back(std::move(new_bins[0]));
            if (success_flag) *success_flag = true;
        } else {
            // 与原逻辑一致：失败则把原 destroy_set 整体作为一个新 bin 追加
            Bin b;
            b.rects = destroy_set;
            temp_bins.push_back(std::move(b));
            *success_flag = false;
        }
//        int sumnum=0;
//        for (int i=0;i<temp_bins.size();i++){
//            sumnum+=temp_bins[i].rects.size();
//        }
//        if (sumnum>100){
//            std::cout<<1;
//        }
        return temp_bins;
    }

    // remain_parts 为空：仅写回被修改的 bins

    *success_flag = true;
    return temp_bins;
}
static inline double sumArea(const std::vector<Rect>& rs) {
    double s = 0.0;
    for (const auto& r : rs) {
        if (r.w == 0.0 || r.h == 0.0) continue;
        s += r.w * r.h;
    }
    return s;
}

struct PackResult {
    std::vector<Bin> bins;      // 已装入的 bins（最多 destroy_num 个）
    double score1{0.0};         // 未装入面积
    double score2{0.0};         // 最后一个 bin 的面积
};

static PackResult pack_with_order(
        std::vector<Rect> parts,        // 已按规则排序/打乱
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
) {
    PackResult out;

    if (destroy_num <= 0) destroy_num = 1;

    std::vector<Bin> bins;
    bins.reserve(std::min<int>(destroy_num, (int)parts.size()));
    bins.push_back(Bin{});

    std::vector<std::vector<Rect>> freeRects;
    freeRects.reserve(bins.capacity());
    freeRects.push_back({ Rect{0,0,binSize.w,binSize.h} });

    int j = 0;
    for (; j < (int)parts.size(); ++j) {
        const Rect& item = parts[j];
        bool placed = false;

        for (size_t k = 0; k < bins.size(); ++k) {
            Rect pos = Insert(item.w, item.h, freeRects[k], mt);
            if (pos.h != 0) {
                bins[k].rects.push_back(pos);
                placed = true;
                break;
            }
        }

        if (!placed) {
            if ((int)bins.size() >= destroy_num) {
                break; // 达到 bin 上限且放不下：停止
            }
            // 新开一个 bin 再放
            bins.push_back(Bin{});
            freeRects.push_back({ Rect{0,0,binSize.w,binSize.h} });

            Rect pos = Insert(item.w, item.h, freeRects.back(), mt);
            if (pos.h == 0) {
                // 理论上不应发生（新 bin 也放不下），这里直接停止
                break;
            }
            bins.back().rects.push_back(pos);
        }
    }

    // score1：未装入面积
    double leftover = 0.0;
    for (int k = j; k < (int)parts.size(); ++k) {
        leftover += parts[k].w * parts[k].h;
    }

    // score2：最后一个 bin 的面积（用来打平）
    double lastBinArea = 0.0;
    if (!bins.empty()) {
        double minArea = std::numeric_limits<double>::infinity();
        for (const auto& b : bins) {
            if (b.rects.empty()) continue;              // 若仍希望跳过空 bin（推荐，否则最小可能变成 0）
            const double a = sumArea(b.rects);
            if (a < minArea) minArea = a;
        }
        if (minArea != std::numeric_limits<double>::infinity())
            lastBinArea = minArea;                      // 这里就是“最小的那个”
    }

    out.bins   = std::move(bins);
    out.score1 = leftover;
    out.score2 = lastBinArea;
    return out;
}
static PackResult pack_with_order_LLM(
        std::vector<Rect> parts,        // 已按规则排序/打乱
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt,
        int method
) {
    PackResult out;
    std::vector<Rect> parts_sorted=parts;
    std::unordered_map<int, PairInfo> pairs;
    std::uniform_real_distribution<double> urand(0.0, 1.0);
    if (urand(mt) <0) {
        // 不 sort：保持 parts_sorted 原顺序；pairs 给一个默认初始化（避免后续使用时报空/未定义）
        pairs.reserve(parts_sorted.size());
        for (size_t i = 0; i < parts_sorted.size(); ++i) {
            PairInfo info;
            // info.rectangles 默认空即可；如需默认包含自身，可 push_back(parts_sorted[i])
            info.position = 0;
            pairs.emplace(static_cast<int>(i), std::move(info));
        }
    } else {
        // 正常：LLM 排序/配对逻辑（假设会在函数内对 parts_sorted 做处理）
        pairs = Sort_by_LLM(parts_sorted, binSize, mt, method);
    }

    if (destroy_num <= 0) destroy_num = 1;

    std::vector<Bin> bins;
    bins.reserve(std::min<int>(destroy_num, (int)parts.size()));
    bins.push_back(Bin{});

    std::vector<std::vector<Rect>> freeRects;
    freeRects.reserve(bins.capacity());
    freeRects.push_back({ Rect{0,0,binSize.w,binSize.h} });

    int j = 0;
    std::uniform_int_distribution<int> algoDist(0, 3);
    int method_ = algoDist(mt);
    for (; j < (int)parts_sorted.size(); ++j) {
        const Rect& item = parts_sorted[j];
        bool placed = false;

        for (size_t k = 0; k < bins.size(); ++k) {
            vector<Rect> pos = Insert_by_LLM(item.w, item.h, freeRects[k], mt, pairs[j],method_);
            if (pos[0].h != 0) {
                for(int i=0;i<pos.size();i++){
                    bins[k].rects.push_back(pos[i]);
                }
                placed = true;
                break;
            }
        }

        if (!placed) {
            if ((int)bins.size() >= destroy_num) {
                break; // 达到 bin 上限且放不下：停止
            }
            // 新开一个 bin 再放
            bins.push_back(Bin{});
            freeRects.push_back({ Rect{0,0,binSize.w,binSize.h} });

            vector<Rect> pos = Insert_by_LLM(item.w, item.h, freeRects.back(), mt, pairs[j],method_);
            if (pos[0].h == 0) {
                // 理论上不应发生（新 bin 也放不下），这里直接停止
                break;
            }
            for(int i=0;i<pos.size();i++){
                bins.back().rects.push_back(pos[i]);
            }
        }
    }

    // score1：未装入面积
    double leftover = 0.0;
    for (int k = j; k < (int)parts.size(); ++k) {
        leftover += parts[k].w * parts[k].h;
    }

    // score2：最后一个 bin 的面积（用来打平）
    double lastBinArea = 0.0;
    if (!bins.empty()) {
        double minArea = std::numeric_limits<double>::infinity();
        for (const auto& b : bins) {
            if (b.rects.empty()) continue;              // 若仍希望跳过空 bin（推荐，否则最小可能变成 0）
            const double a = sumArea(b.rects);
            if (a < minArea) minArea = a;
        }
        if (minArea != std::numeric_limits<double>::infinity())
            lastBinArea = minArea;                      // 这里就是“最小的那个”
    }

    out.bins   = std::move(bins);
    out.score1 = leftover;
    out.score2 = lastBinArea;
    return out;
}
// ==== LLM_SORT_BEGIN ====

std::unordered_map<int, PairInfo> Sort_by_LLM(std::vector<Rect>& parts_sorted, const Rect& binSize, std::mt19937& mt, int method) {
    std::unordered_map<int, PairInfo> pairs;
    // --- Algorithm 0: Area sort, pair smallest horizontally where fits ---
    auto algo_area_pair_smallest = [&](std::vector<Rect>& rects, std::unordered_map<int, PairInfo>& out_pairs) {
        std::sort(rects.begin(), rects.end(), cmpArea);
        // Pair smallest (back of sorted) horizontally if they fit
        std::vector<bool> paired(rects.size(), false);
        std::vector<Rect> out_rects;
        std::vector<PairInfo> pairinfos;
        // Pair from back (smallest)
        int n = rects.size();
        int left = 0, right = n-1;
        while (left <= right) {
            if (left == right) {
                // Only one left unpaired
                PairInfo info;
                info.rectangles = {rects[left]};
                info.position = -1;
                out_rects.push_back(rects[left]);
                pairinfos.push_back(info);
                left++;
            } else {
                // Try to pair rects[right] and rects[right-1] horizontally
                double wsum = rects[right].w + rects[right-1].w;
                double hmax = std::max(rects[right].h, rects[right-1].h);
                if (wsum <= binSize.w && hmax <= binSize.h) {
                    Rect comb = {0,0,wsum,hmax};
                    PairInfo info;
                    info.rectangles = {rects[right-1], rects[right]};
                    info.position = 0; // horizontal
                    out_rects.push_back(comb);
                    pairinfos.push_back(info);
                    right -= 2;
                } else {
                    // Only pair rects[right] as single
                    PairInfo info;
                    info.rectangles = {rects[right]};
                    info.position = -1;
                    out_rects.push_back(rects[right]);
                    pairinfos.push_back(info);
                    right -= 1;
                }
            }
        }
        parts_sorted = out_rects;
        for (int i = 0; i < pairinfos.size(); ++i) {
            out_pairs[i] = pairinfos[i];
        }
    };
    // --- Algorithm 1: Long side sort, pair for best bounding box compactness (horizontal or vertical) ---
    auto algo_longside_compactpair = [&](std::vector<Rect>& rects, std::unordered_map<int, PairInfo>& out_pairs) {
        std::sort(rects.begin(), rects.end(), cmpLongSide);
        std::vector<bool> used(rects.size(), false);
        std::vector<Rect> out_rects;
        std::vector<PairInfo> pairinfos;
        for (int i = 0; i < (int)rects.size(); ++i) {
            if (used[i]) continue;
            double best_area = std::numeric_limits<double>::max();
            int best_j = -1, best_pos = 0;
            Rect best_bbox = {0,0,0,0};
            for (int j = i+1; j < (int)rects.size(); ++j) {
                if (used[j]) continue;
                // Try horizontal pair
                double w_h = rects[i].w + rects[j].w;
                double h_h = std::max(rects[i].h, rects[j].h);
                if (w_h <= binSize.w && h_h <= binSize.h) {
                    double area = w_h * h_h;
                    if (area < best_area) {
                        best_area = area;
                        best_j = j;
                        best_pos = 0;
                        best_bbox = {0,0,w_h,h_h};
                    }
                }
                // Try vertical pair
                double w_v = std::max(rects[i].w, rects[j].w);
                double h_v = rects[i].h + rects[j].h;
                if (w_v <= binSize.w && h_v <= binSize.h) {
                    double area = w_v * h_v;
                    if (area < best_area) {
                        best_area = area;
                        best_j = j;
                        best_pos = 1;
                        best_bbox = {0,0,w_v,h_v};
                    }
                }
            }
            if (best_j != -1) {
                used[i] = used[best_j] = true;
                PairInfo info;
                info.rectangles = {rects[i], rects[best_j]};
                info.position = best_pos;
                out_rects.push_back(best_bbox);
                pairinfos.push_back(info);
            } else {
                used[i] = true;
                PairInfo info;
                info.rectangles = {rects[i]};
                info.position = -1;
                out_rects.push_back(rects[i]);
                pairinfos.push_back(info);
            }
        }
        parts_sorted = out_rects;
        for (int i = 0; i < pairinfos.size(); ++i) {
            out_pairs[i] = pairinfos[i];
        }
    };
    // --- Algorithm 2: Shortest side sort, pair same height or same width where fits ---
    auto algo_shortside_samesizepair = [&](std::vector<Rect>& rects, std::unordered_map<int, PairInfo>& out_pairs) {
        std::sort(rects.begin(), rects.end(), cmpShortSide);
        std::vector<bool> paired(rects.size(), false);
        std::vector<Rect> out_rects;
        std::vector<PairInfo> pairinfos;
        int n = rects.size();
        for (int i=0; i<n; ++i) {
            if (paired[i]) continue;
            bool madepair = false;
            for (int j=i+1; j<n; ++j) {
                if (paired[j]) continue;
                // Horizontal (if exact same height)
                if (std::abs(rects[i].h-rects[j].h)<1e-7 && rects[i].w + rects[j].w <= binSize.w) {
                    Rect comb = {0,0,rects[i].w+rects[j].w, rects[i].h};
                    PairInfo info;
                    info.rectangles = {rects[i],rects[j]};
                    info.position = 0;
                    out_rects.push_back(comb);
                    pairinfos.push_back(info);
                    paired[i]=paired[j]=true;
                    madepair = true;
                    break;
                }
                // Vertical (if exact same width)
                if (std::abs(rects[i].w-rects[j].w)<1e-7 && rects[i].h + rects[j].h <= binSize.h) {
                    Rect comb = {0,0,rects[i].w, rects[i].h+rects[j].h};
                    PairInfo info;
                    info.rectangles = {rects[i],rects[j]};
                    info.position = 1;
                    out_rects.push_back(comb);
                    pairinfos.push_back(info);
                    paired[i]=paired[j]=true;
                    madepair = true;
                    break;
                }
            }
            if (!madepair) {
                PairInfo info;
                info.rectangles = {rects[i]};
                info.position = -1;
                out_rects.push_back(rects[i]);
                pairinfos.push_back(info);
                paired[i]=true;
            }
        }
        parts_sorted = out_rects;
        for (int i = 0; i < pairinfos.size(); ++i) {
            out_pairs[i] = pairinfos[i];
        }
    };

    int select_method;
    if (method >= 0 && method <=2)
        select_method = method;
    else {
        std::uniform_int_distribution<int> dist(0,2);
        select_method = dist(mt);
    }
    std::vector<Rect> input_copy = parts_sorted;
    std::unordered_map<int, PairInfo> out_pairs;
    if (select_method == 0) {
        algo_area_pair_smallest(input_copy, out_pairs);
    } else if (select_method == 1) {
        algo_longside_compactpair(input_copy, out_pairs);
    } else {
        algo_shortside_samesizepair(input_copy, out_pairs);
    }
    parts_sorted = input_copy;
    pairs = out_pairs;
    return pairs;
}

// ==== LLM_SORT_END ====
// ==== LLM_PLACE_BEGIN ====

Rect FindPosition_by_LLM(double width, double height, std::vector<Rect> freeRectangles, std::mt19937& mt, int method) {
    Rect node{0,0,0,0};

    // Method 0: Best Short Side Fit  
    auto best_short_side_fit = [&]() -> Rect {
        Rect bestNode{0,0,0,0};
        double bestShortSideFit = std::numeric_limits<double>::max();
        double bestLongSideFit = std::numeric_limits<double>::max();
        std::uniform_real_distribution<double> dt(0.0, 1.0);
        bool bestFound = false;
        for(const auto& fr : freeRectangles) {
            if(fr.w >= width && fr.h >= height) {
                double leftoverHoriz = std::abs(fr.w - width);
                double leftoverVert  = std::abs(fr.h - height);
                double shortSideFit  = std::min(leftoverHoriz, leftoverVert);
                double longSideFit   = std::max(leftoverHoriz, leftoverVert);
                if (shortSideFit < bestShortSideFit) {
                    bestShortSideFit = shortSideFit;
                    bestLongSideFit  = longSideFit;
                    bestNode = Rect{fr.x, fr.y, width, height};
                    bestFound = true;
                } else if (shortSideFit == bestShortSideFit) {
                    if (longSideFit < bestLongSideFit) {
                        bestLongSideFit = longSideFit;
                        bestNode = Rect{fr.x, fr.y, width, height};
                        bestFound = true;
                    } else if (longSideFit == bestLongSideFit) {
                        // Mild randomness in tie
                        if(dt(mt) < 0.5) {
                            bestNode = Rect{fr.x, fr.y, width, height};
                            bestFound = true;
                        }
                    }
                }
            }
        }
        if(!bestFound) return Rect{0,0,0,0};
        return bestNode;
    };

    // Method 1: First Area Fit (requires area >= 2 * needed)
    auto first_area_fit = [&]() -> Rect {
        for(const auto& fr : freeRectangles) {
            if(fr.w >= width && fr.h >= height) {
                double freeArea = fr.w*fr.h;
                double reqArea = width*height;
                if(freeArea >= reqArea*2) { // area at least twice
                    return Rect{fr.x, fr.y, width, height};
                }
            }
        }
        // If none are at least double-area, fallback: just pick first fit
        for(const auto& fr : freeRectangles) {
            if(fr.w >= width && fr.h >= height) {
                return Rect{fr.x, fr.y, width, height};
            }
        }
        return Rect{0,0,0,0};
    };

    // Method 2: Top-Left-Most (break ties by min unused perimeter)
    auto top_left_most = [&]() -> Rect {
        double bestY = std::numeric_limits<double>::max();
        double bestX = std::numeric_limits<double>::max();
        double bestUnusedPeri = std::numeric_limits<double>::max();
        bool found = false;
        Rect bestNode{0,0,0,0};
        for(const auto& fr : freeRectangles) {
            if(fr.w >= width && fr.h >= height) {
                double unusedPeri = (fr.w-fr.x-width)+(fr.h-fr.y-height);
                if(fr.y < bestY || (fr.y == bestY && fr.x < bestX)) {
                    bestY = fr.y;
                    bestX = fr.x;
                    bestUnusedPeri = unusedPeri;
                    bestNode = Rect{fr.x, fr.y, width, height};
                    found = true;
                } else if (fr.y == bestY && fr.x == bestX) {
                    if(unusedPeri < bestUnusedPeri) {
                        bestUnusedPeri = unusedPeri;
                        bestNode = Rect{fr.x, fr.y, width, height};
                        found = true;
                    }
                }
            }
        }
        if(!found) return Rect{0,0,0,0};
        return bestNode;
    };

    int select_method = method;
    if(select_method < 0 || select_method > 2) {
        std::uniform_int_distribution<int> md(0,2);
        select_method = md(mt);
    }
    if(select_method == 0) {
        node = best_short_side_fit();
    } else if(select_method == 1) {
        node = first_area_fit();
    } else if(select_method == 2) {
        node = top_left_most();
    }
    return node;
}

// ==== LLM_PLACE_END ====
std::vector<Rect> Insert_by_LLM(
        double width,
        double height,
        std::vector<Rect>& freeRectangles,
        std::mt19937& mt,
        PairInfo pair,
        int method
) {
    Rect newNode;
    if (method==0){
        newNode = FindPositionForNewNodeBestShortSideFit(width, height, freeRectangles, mt);
    }
    else {
        newNode = FindPosition_by_LLM(width, height, freeRectangles, mt,method);
    }
    if (newNode.h == 0) return {newNode};
    if (pair.rectangles.size()>1){
        Rect newNode_1;
        Rect newNode_2;
        for(int i=0;i<pair.rectangles.size();i++){
            if(pair.position==0){
                newNode_1=Rect{newNode.x,newNode.y,pair.rectangles[0].w,pair.rectangles[0].h};
                newNode_2=Rect{newNode.x+pair.rectangles[0].w,newNode.y,pair.rectangles[1].w,pair.rectangles[1].h};
            }
            else{
                newNode_1=Rect{newNode.x,newNode.y,pair.rectangles[0].w,pair.rectangles[0].h};
                newNode_2=Rect{newNode.x,newNode.y+pair.rectangles[0].h,pair.rectangles[1].w,pair.rectangles[1].h};
            }
        }
        for (size_t i = 0; i < freeRectangles.size(); ++i) {
            if (SplitFreeNode(freeRectangles[i], newNode_1, freeRectangles)) {
                freeRectangles.erase(freeRectangles.begin() + i);
                --i; // 继续检查当前位置（因为 erase 使后续左移）
            }
        }
        PruneFreeList(freeRectangles);
        for (size_t i = 0; i < freeRectangles.size(); ++i) {
            if (SplitFreeNode(freeRectangles[i], newNode_2, freeRectangles)) {
                freeRectangles.erase(freeRectangles.begin() + i);
                --i; // 继续检查当前位置（因为 erase 使后续左移）
            }
        }
        PruneFreeList(freeRectangles);
        return {newNode_1,newNode_2};
    }
    else{
        for (size_t i = 0; i < freeRectangles.size(); ++i) {
            if (SplitFreeNode(freeRectangles[i], newNode, freeRectangles)) {
                freeRectangles.erase(freeRectangles.begin() + i);
                --i; // 继续检查当前位置（因为 erase 使后续左移）
            }
        }
        PruneFreeList(freeRectangles);
        return {newNode};
    }

    // 这里不要反复调用 freeRectangles.size()（虽小，但免费优化）
}
// 结合 maxRectsSorts 的逻辑：遍历 6 种规则，选 (score1, score2) 最小
std::vector<Bin> regroup_4(
        std::vector<Rect> destroy_set,
        bool* success_flag,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
) {
    if (success_flag) *success_flag = false;
    if (destroy_set.empty()) return {};

    std::vector<Bin> bestBins;
    double bestScore1 = std::numeric_limits<double>::infinity();
    double bestScore2 = std::numeric_limits<double>::infinity();

    for (int rule = 0; rule < 5; ++rule) {
        // 每个 rule 从原始 destroy_set 重新开始（避免被上一次排序污染）
//        std::vector<Rect> parts = destroy_set;
//
//        switch (rule) {
//            case 0: std::sort(parts.begin(), parts.end(), cmpArea); break;
//            case 1: std::sort(parts.begin(), parts.end(), cmpPerimeter); break;
//            case 3: std::sort(parts.begin(), parts.end(), cmpLongSide); break;
//            case 4: std::sort(parts.begin(), parts.end(), cmpShortSide); break;
//            case 2: break;
////            default: std::shuffle(parts.begin(), parts.end(), mt); break; // rule 4/5
//        }
//
//        // 注意：pack_with_order 内部会继续调用 Insert（Insert 可能用随机打平），会消耗 mt
//        PackResult cand = pack_with_order(std::move(parts), binSize, destroy_num, mt);
        std::vector<Rect> parts = destroy_set;
        PackResult cand = pack_with_order_LLM(std::move(parts), binSize, destroy_num, mt,rule);

        if (cand.score1 < bestScore1 || (cand.score1 == bestScore1 && cand.score2 < bestScore2)) {
            bestScore1 = cand.score1;
            bestScore2 = cand.score2;
            bestBins = std::move(cand.bins);
        }

    }

    if (success_flag) *success_flag = (bestScore1 == 0.0);
//    int sumnum=0;
//    for (int i=0;i<bestBins.size();i++){
//        sumnum+=bestBins[i].rects.size();
//    }
//    if (sumnum>destroy_set.size()){
//        std::cout<<1;
//    }
    return bestBins;
}
//std::vector<Bin> regroup_4(
//        std::vector<Rect> destroy_set,   // 按值：我们要 sort/shuffle
//        bool* success_flag,
//        const Rect& binSize,             // 用 Rect 承载 w,h（x,y 忽略）
//        int destroy_num,
//        std::mt19937& mt
//) {
//    std::uniform_int_distribution<int> dist(0, 5);
//    const int rule = dist(mt);
//
//    switch (rule) {
//        case 0: std::sort(destroy_set.begin(), destroy_set.end(), cmpArea); break;
//        case 1: std::sort(destroy_set.begin(), destroy_set.end(), cmpPerimeter); break;
//        case 2: std::sort(destroy_set.begin(), destroy_set.end(), cmpLongSide); break;
//        case 3: std::sort(destroy_set.begin(), destroy_set.end(), cmpShortSide); break;
//        case 4: [[fallthrough]];
//        case 5: std::shuffle(destroy_set.begin(), destroy_set.end(), mt); break;
//    }
//
//    std::vector<Bin> bins;
//    bins.reserve(std::min<int>(destroy_num > 0 ? destroy_num : 1, (int)destroy_set.size()));
//    bins.push_back(Bin{});
//
//    // 每个 bin 对应一份 free-rectangles 列表
//    std::vector<std::vector<Rect>> freeRects;
//    freeRects.reserve(bins.capacity());
//    freeRects.push_back({Rect{0, 0, binSize.w, binSize.h}});
//
//    for (const auto& item : destroy_set) {
//        bool placed = false;
//
//        // 尝试放入已有 bins
//        for (size_t k = 0; k < bins.size(); ++k) {
//            Rect pos = Insert(item.w, item.h, freeRects[k], mt);
//            if (pos.h != 0) {
//                // 注意：原 destroy_set[j] 里可能还带“其他字段”，但你现在结构已经明确为 Rect
//                bins[k].rects.push_back(pos);
//                placed = true;
//                break;
//            }
//        }
//
//        if (!placed) {
//            // 如果你想启用“超过 destroy_num 就停止”的逻辑，可在这里恢复你注释掉的 break
//            bins.push_back(Bin{});
//            freeRects.push_back({Rect{0, 0, binSize.w, binSize.h}});
//
//            Rect pos = Insert(item.w, item.h, freeRects.back(), mt);
//            // 你原代码这里是 push_back(destroy_set[j])，会把“未放置坐标”的矩形塞进去，不严谨。
//            // 建议：能放就放 pos；不能放就仍然记录 item（或直接失败返回）。
//            if (pos.h != 0) bins.back().rects.push_back(pos);
//            else            bins.back().rects.push_back(item);
//        }
//    }
//
//    if (success_flag) *success_flag = (destroy_num <= 0) ? true : ((int)bins.size() <= destroy_num);
//    return bins;
//}

std::vector<Bin> maxRects(
        std::vector<Rect> parts,     // 按值：内部要 shuffle + erase
        double width,
        double height,
        int binNum,
        float* score1,
        float* score2,
        std::mt19937& mt             // 关键：按引用
) {
//    std::shuffle(parts.begin(), parts.end(), mt);

    // 保持你原来 max(binNum-1,1) 的逻辑
    const int initBins = std::max(binNum - 1, 1);

    std::vector<Bin> ret(initBins);
    std::vector<std::vector<Rect>> freeRectangles(initBins, std::vector<Rect>{ Rect{0, 0, width, height} });

    double bestScore1 = std::numeric_limits<double>::max();
    int big_length = 1; // 只在前 big_length 个 bin 中尝试插入（与你原逻辑一致）

    // 复用候选容器，避免每轮频繁分配
    std::vector<std::pair<int,int>> candidates;
    candidates.reserve(128);

    while (!parts.empty()) {
        bestScore1 = std::numeric_limits<double>::max();
        candidates.clear();

        const int activeBins = std::min<int>(big_length, (int)freeRectangles.size());

        // 1) 在当前 activeBins 中找 bestScore1 的候选 (binIndex, partIndex)
        for (int i1 = 0; i1 < (int)parts.size(); ++i1) {
            const double pw = parts[i1].w;
            const double ph = parts[i1].h;

            for (int j = 0; j < activeBins; ++j) {
                const double sc = Insertscore(pw, ph, freeRectangles[j], mt);

                if (sc < bestScore1) {
                    bestScore1 = sc;
                    candidates.clear();
                    candidates.emplace_back(j, i1);
                } else if (sc == bestScore1 && sc != std::numeric_limits<double>::max()) {
                    candidates.emplace_back(j, i1);
                }
            }
        }

        // 2) 若当前 activeBins 都放不下任何 part
        if (bestScore1 == std::numeric_limits<double>::max()) {
            if (big_length < binNum - 1) {
                // 扩大搜索到更多“已有 bin”
                ++big_length;
                continue;
            }

            // 否则：新增一个 bin（与你原逻辑一致：big_length 也加 1）
            ++big_length;
            freeRectangles.push_back({ Rect{0, 0, width, height} });
            ret.push_back(Bin{});

            // 在新 bin 中挑一个最容易插的 part
            int bestRectIndex = -1;
            double bestLocal = std::numeric_limits<double>::max();

            for (int i1 = 0; i1 < (int)parts.size(); ++i1) {
                const double sc = Insertscore(parts[i1].w, parts[i1].h, freeRectangles.back(), mt);
                if (sc < bestLocal) {
                    bestLocal = sc;
                    bestRectIndex = i1;
                }
            }

            if (bestRectIndex != -1 && bestLocal != std::numeric_limits<double>::max()) {
                // 真正插入：Insert 会更新 freeRectangles.back()
                Rect placed = Insert(parts[bestRectIndex].w, parts[bestRectIndex].h, freeRectangles.back(), mt);
                if (placed.h != 0) {
                    ret.back().rects.push_back(placed);
                    parts.erase(parts.begin() + bestRectIndex);
                }
            }

            continue;
        }

        // 3) 存在可行候选：在 candidates 中随机选一个（与你原逻辑一致）
        std::uniform_int_distribution<int> dis(0, (int)candidates.size() - 1);
        const int pick = dis(mt);
        const int bestBinIndex  = candidates[pick].first;
        const int bestRectIndex = candidates[pick].second;

        Rect placed = Insert(parts[bestRectIndex].w, parts[bestRectIndex].h, freeRectangles[bestBinIndex], mt);
        if (placed.h != 0) {
            ret[bestBinIndex].rects.push_back(placed);
            parts.erase(parts.begin() + bestRectIndex);
        } else {
            // 理论上 Insertscore 过滤后不应发生；保底：避免死循环
            parts.erase(parts.begin() + bestRectIndex);
        }
    }

    // 4) 计算 score1 / score2（保持你原始定义）
    if (score1) {
        *score1 = 0.0f;
        for (int i = 0; i < (int)ret.size(); ++i) {
            if (i >= binNum) {
                for (const auto& p : ret[i].rects) *score1 += (float)(p.w * p.h);
            }
        }
        // 原代码这里还加 parts 的面积；但 while 结束 parts 已空，这里等价为 0
    }

    if (score2) {
        *score2 = 0.0f;
        if (!ret.empty()) {
            for (const auto& p : ret.back().rects) *score2 += (float)(p.w * p.h);
        }
    }

    return ret;
}
double Insertscore(
        double width,
        double height,
        const std::vector<Rect>& freeRectangles,  // 关键：const&
        std::mt19937& mt                           // 关键：&
) {
    double score1 = std::numeric_limits<double>::max();
    double score2 = std::numeric_limits<double>::max();

    // newNode 若你后续完全不用，可不接收返回值
    Rect node=FindPositionForNewNodeBestShortSideFit_(width, height, score1, score2, freeRectangles, mt);

    return score1;
}
Rect FindPositionForNewNodeBestShortSideFit_(double width, double height,
                                       double &bestShortSideFit,
                                       double &bestLongSideFit,std::vector<Rect> freeRectangles,std::mt19937 mt) {
    Rect bestNode{0,0,0,0};

    // 仅在需要随机打破平局时再取随机数，避免每次构造 distribution
    std::uniform_real_distribution<double> dt(0.0, 1.0);

    for (const auto& fr : freeRectangles) {
        if (fr.w >= width && fr.h >= height) {
            const double leftoverHoriz = std::abs(fr.w - width);
            const double leftoverVert  = std::abs(fr.h - height);
            const double shortSideFit  = std::min(leftoverHoriz, leftoverVert);
            const double longSideFit   = std::max(leftoverHoriz, leftoverVert);

            if (shortSideFit < bestShortSideFit) {
                bestNode = Rect{fr.x, fr.y, width, height};
                bestShortSideFit = shortSideFit;
                bestLongSideFit  = longSideFit;
            } else if (shortSideFit == bestShortSideFit && longSideFit < bestLongSideFit) {
                // 你原先的 tie-break：50% 概率替换
                bestNode = Rect{fr.x, fr.y, width, height};
                bestShortSideFit = shortSideFit;
                bestLongSideFit = longSideFit;
            }
                else if (shortSideFit == bestShortSideFit && longSideFit == bestLongSideFit){
                    if (dt(mt) > 0.5) {
                        bestNode = Rect{fr.x, fr.y, width, height};
                        bestShortSideFit = shortSideFit;
                        bestLongSideFit  = longSideFit;
                    }
                }
            }
        }
    return bestNode;
}
Rect FindPositionForNewNodeBestShortSideFit(
        double width,
        double height,
        const std::vector<Rect>& freeRectangles,
        std::mt19937& mt
) {
    Rect bestNode{0,0,0,0};
    double bestShortSideFit = std::numeric_limits<double>::max();
    double bestLongSideFit  = std::numeric_limits<double>::max();

    // 仅在需要随机打破平局时再取随机数，避免每次构造 distribution
    std::uniform_real_distribution<double> dt(0.0, 1.0);

    for (const auto& fr : freeRectangles) {
        if (fr.w >= width && fr.h >= height) {
            const double leftoverHoriz = std::abs(fr.w - width);
            const double leftoverVert  = std::abs(fr.h - height);
            const double shortSideFit  = std::min(leftoverHoriz, leftoverVert);
            const double longSideFit   = std::max(leftoverHoriz, leftoverVert);

            if (shortSideFit < bestShortSideFit) {
                bestNode = Rect{fr.x, fr.y, width, height};
                bestShortSideFit = shortSideFit;
                bestLongSideFit  = longSideFit;
            } else if (shortSideFit == bestShortSideFit && longSideFit < bestLongSideFit) {
                // 你原先的 tie-break：50% 概率替换
                bestNode = Rect{fr.x, fr.y, width, height};
                bestShortSideFit = shortSideFit;
                bestLongSideFit = longSideFit;
            }
            else if (shortSideFit == bestShortSideFit && longSideFit == bestLongSideFit){
                if (dt(mt) > 0.5) {
                    bestNode = Rect{fr.x, fr.y, width, height};
                    bestShortSideFit = shortSideFit;
                    bestLongSideFit  = longSideFit;
                }
            }
        }
    }
    return bestNode;
}

Rect Insert(
        double width,
        double height,
        std::vector<Rect>& freeRectangles,
        std::mt19937& mt
) {
    Rect newNode = FindPositionForNewNodeBestShortSideFit(width, height, freeRectangles, mt);
    if (newNode.h == 0) return newNode;

    // 这里不要反复调用 freeRectangles.size()（虽小，但免费优化）
    for (size_t i = 0; i < freeRectangles.size(); ++i) {
        if (SplitFreeNode(freeRectangles[i], newNode, freeRectangles)) {
            freeRectangles.erase(freeRectangles.begin() + i);
            --i; // 继续检查当前位置（因为 erase 使后续左移）
        }
    }
    PruneFreeList(freeRectangles);
    return newNode;
}

static PackResult pack_with_order_and_remain(
        std::vector<Rect> parts,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt,
        std::vector<Rect>& remain   // <<< 新增：输出剩余矩形
) {
    PackResult out;
    remain.clear();

    if (destroy_num <= 0) destroy_num = 1;

    std::vector<Bin> bins;
    bins.reserve(std::min<int>(destroy_num, (int)parts.size()));
    bins.push_back(Bin{});

    std::vector<std::vector<Rect>> freeRects;
    freeRects.reserve(bins.capacity());
    freeRects.push_back({ Rect{0,0,binSize.w,binSize.h} });

    int j = 0;
    for (; j < (int)parts.size(); ++j) {
        const Rect& item = parts[j];
        bool placed = false;

        for (size_t k = 0; k < bins.size(); ++k) {
            Rect pos = Insert(item.w, item.h, freeRects[k], mt);
            if (pos.h != 0) {
                bins[k].rects.push_back(pos);
                placed = true;
                break;
            }
        }

        if (!placed) {
            if ((int)bins.size() >= destroy_num) {
                break;
            }
            bins.push_back(Bin{});
            freeRects.push_back({ Rect{0,0,binSize.w,binSize.h} });

            Rect pos = Insert(item.w, item.h, freeRects.back(), mt);
            if (pos.h == 0) break;
            bins.back().rects.push_back(pos);
        }
    }

    // 剩余矩形 r：parts[j..end)
    if (j < (int)parts.size()) {
        remain.insert(remain.end(), parts.begin() + j, parts.end());
    }

    // score1：未装入面积
    double leftover = 0.0;
    for (int k = j; k < (int)parts.size(); ++k) {
        leftover += parts[k].w * parts[k].h;
    }

    // score2：最小装载 bin 的面积（与你当前 pack_with_order 一致）
    double lastBinArea = 0.0;
    if (!bins.empty()) {
        double minArea = std::numeric_limits<double>::infinity();
        for (const auto& b : bins) {
            if (b.rects.empty()) continue;
            const double a = sumArea(b.rects);
            if (a < minArea) minArea = a;
        }
        if (minArea != std::numeric_limits<double>::infinity())
            lastBinArea = minArea;
    }

    out.bins   = std::move(bins);
    out.score1 = leftover;
    out.score2 = lastBinArea;
    return out;
}
std::vector<Rect> insertparts_1(
        std::vector<Rect>& bin,
        const Rect& part,
        const Rect& binSize,
        std::mt19937& mt,
        bool *flag
) {
    std::vector<Rect> parts;
    parts.reserve(bin.size() + 1);
    parts.insert(parts.end(), bin.begin(), bin.end());
    parts.push_back(part);

    // 50% 概率保持原顺序；50% 概率应用 rule（sort 或 shuffle）
    float best1 = 1e30f, best2 = 1e30f;
    std::vector<Rect> bestParts = parts;
    std::vector<Bin>  bestRet;
    std::mt19937      bestMt = mt;
    std::vector<Rect> bestRemain;
    auto try_one = [&](std::vector<Rect> cand) {
        std::mt19937 mtLocal = mt;  // 每个候选从同一 RNG 状态出发，公平比较

        std::vector<Rect> remain;
        PackResult pr = pack_with_order_and_remain(std::move(cand), binSize, /*destroy_num=*/1, mtLocal, remain);

        const double s1 = pr.score1;
        const double s2 = pr.score2;

        if (s1 < best1 || (s1 == best1 && s2 < best2)) {
            best1 = s1; best2 = s2;
            // 注意：cand 已经 move 进 pack_with_order_and_remain 了，如果你确实要保存“该候选顺序”，
            // 需要在调用前保留一份；否则可删掉 bestParts 这一行。
            bestRet    = std::move(pr.bins);
            bestRemain = std::move(remain);  // <<< 这就是你要的“剩余矩形 r”
            bestMt     = mtLocal;
        }
    };

    // 0) 原顺序
    try_one(parts);

    // 1) 四种排序
    { auto cand = parts; std::sort(cand.begin(), cand.end(), cmpArea);      try_one(std::move(cand)); }
    { auto cand = parts; std::sort(cand.begin(), cand.end(), cmpPerimeter); try_one(std::move(cand)); }
    { auto cand = parts; std::sort(cand.begin(), cand.end(), cmpLongSide);  try_one(std::move(cand)); }
    { auto cand = parts; std::sort(cand.begin(), cand.end(), cmpShortSide); try_one(std::move(cand)); }

    // 2) 两次 shuffle（对应你原来 case 4/5 都是 shuffle）
    {;        // 用一个滚动 RNG 生成两个不同 shuffle
        for (int t = 0; t < 2; ++t) {
            auto cand = parts;
            std::shuffle(cand.begin(), cand.end(), mt);
            try_one(std::move(cand));
        }
    }

    // 写回最佳结果
    parts = std::move(bestParts);

    float score1MaxRects = best1, score2MaxRects = best2;
    std::vector<Bin> ret = std::move(bestRet);

    if (score1MaxRects == 0.0f && !ret.empty()) {
        bin = std::move(ret[0].rects);
        *flag=true;
        return {};
    }
    else{
        if (sumArea(bin)<sumArea(ret[0].rects)){
            bin = std::move(ret[0].rects);
//            std::vector<Rect> out;
//            size_t cap = 0;
//            for (size_t i = 1; i < ret.size(); ++i) cap += ret[i].rects.size();
//            out.reserve(cap);
//
//            for (size_t i = 1; i < ret.size(); ++i) {
//                const auto& rs = ret[i].rects;
//                out.insert(out.end(), rs.begin(), rs.end());
//            }
            *flag=true;
            return bestRemain;
        }
    }
    *flag=false;
    return {part};
}

bool SplitFreeNode(Rect freeNode, const Rect &usedNode,std::vector<Rect> &freeRectangles) {
    // Test with SAT if the rectangles even intersect.
    if (usedNode.x >= freeNode.x + freeNode.w || usedNode.x + usedNode.w <= freeNode.x ||
        usedNode.y >= freeNode.y + freeNode.h || usedNode.y + usedNode.h <= freeNode.y)
        return false;

    if (usedNode.x < freeNode.x + freeNode.w && usedNode.x + usedNode.w > freeNode.x) {
        // New node at the top side of the used node.
        if (usedNode.y > freeNode.y && usedNode.y < freeNode.y + freeNode.h) {
            Rect newNode = freeNode;
            newNode.h = usedNode.y - newNode.y;
            freeRectangles.push_back(newNode);
        }

        // New node at the bottom side of the used node.
        if (usedNode.y + usedNode.h < freeNode.y + freeNode.h) {
            Rect newNode = freeNode;
            newNode.y = usedNode.y + usedNode.h;
            newNode.h = freeNode.y + freeNode.h - (usedNode.y + usedNode.h);
            freeRectangles.push_back(newNode);
        }
    }

    if (usedNode.y < freeNode.y + freeNode.h && usedNode.y + usedNode.h > freeNode.y) {
        // New node at the left side of the used node.
        if (usedNode.x > freeNode.x && usedNode.x < freeNode.x + freeNode.w) {
            Rect newNode = freeNode;
            newNode.w = usedNode.x - newNode.x;
            freeRectangles.push_back(newNode);
        }

        // New node at the right side of the used node.
        if (usedNode.x + usedNode.w < freeNode.x + freeNode.w) {
            Rect newNode = freeNode;
            newNode.x = usedNode.x + usedNode.w;
            newNode.w = freeNode.x + freeNode.w - (usedNode.x + usedNode.w);
            freeRectangles.push_back(newNode);
        }
    }

    return true;
}
void PruneFreeList(std::vector<Rect> &freeRectangles) {
    for (size_t i = 0; i < freeRectangles.size(); ++i)
        for (size_t j = i + 1; j < freeRectangles.size(); ++j) {
            if (IsContainedIn(freeRectangles[i], freeRectangles[j])) {
                freeRectangles.erase(freeRectangles.begin() + i);
                --i;
                break;
            }
            if (IsContainedIn(freeRectangles[j], freeRectangles[i])) {
                freeRectangles.erase(freeRectangles.begin() + j);
                --j;
            }
        }
}
bool IsContainedIn(const Rect &a, const Rect &b)
    {
        return a.x >= b.x && a.y >= b.y
               && a.x+a.w <= b.x+b.w
                  && a.y+a.h <= b.y+b.h;
    }
std::vector<Bin> acceptance(
        const std::vector<Bin>& temp_solution,
        double& solution_cost,
        const std::vector<Bin>& solution,
        int iter,
        int max_iter,
        const Rect& binSize,
        std::mt19937& mt,
        bool flag,
        int failiter

) {
    // 原逻辑：flag=false 直接拒绝
    if (!flag) return solution;

    // theta=0.07*(max_iter-iter)/max_iter
    // 避免整型除法问题，且防止 max_iter==0
    double theta = (max_iter > 0)
                         ? 0.2 * (double)(max_iter - iter) / (double)max_iter
                         : 0.0;
    if (failiter>max_iter/10&& theta<0.1 * (double)(max_iter - iter) / (double)max_iter){
        theta=theta+0.1 * (double)(max_iter - iter) / (double)max_iter;
    }
    const double new_cost = getcost(temp_solution, binSize);
//    int sumnum=0;
//    for (int i=0;i<temp_solution.size();i++){
//        sumnum+=temp_solution[i].rects.size();
//    }
//    if (sumnum>100){
//        std::cout<<1;
//    }
    // 原逻辑：new_solution_cost > solution_cost - theta 则接受
    std::uniform_real_distribution<> dis(0, 1);
    float prob = dis(mt);

    if (prob < exp((-solution_cost +new_cost) / theta)) {
        solution_cost = new_cost;
        std::vector<Bin> accepted = temp_solution;
        accepted.erase(
                std::remove_if(accepted.begin(), accepted.end(),
                               [](const Bin& b) { return b.rects.empty(); }),
                accepted.end()
        );
        return accepted;
    }
/*    if (new_cost > solution_cost - theta) {
        solution_cost = new_cost;

        // 等价于你原来的“接受后删空 bin”
        std::vector<Bin> accepted = temp_solution;
        accepted.erase(
                std::remove_if(accepted.begin(), accepted.end(),
                               [](const Bin& b) { return b.rects.empty(); }),
                accepted.end()
        );
        return accepted;
    }*/

    return solution;
}
double getcost(const std::vector<Bin>& solution, const Rect& binSize) {
    const double binfill = binSize.w * binSize.h;
    if (binfill <= 0.0) return 0.0;

    const double coef = 4.0;  // 等价于 4/binfill/binfill

    double cost = 0.0;

    int lastIndex = -1;
    double lastAvg = std::numeric_limits<double>::infinity();
    double lastSumAreaSq = 0.0; // 用于最后一次性扣掉“最小平均面积bin”的项

    for (size_t i = 0; i < solution.size(); ++i) {
        const auto& bin = solution[i];
        if (bin.rects.empty()) continue;   // 对应你原来的 temp_solution[i].size()==0

        cost -= 10.0;

        double sumArea = 0.0;
        double sumAreaSq = 0.0;
        int cnt = 0;

        for (const auto& r : bin.rects) {
            // 你原代码是通过 "size()>0" 判断是否有效；这里用 h==0 作为无效约定（与前面一致）
            if (r.h == 0.0 || r.w == 0.0) continue;

            const double a = r.w * r.h;
            ++cnt;
            sumArea += a;
        }
        sumAreaSq = sumArea * sumArea/ (binfill * binfill);

        if (cnt > 0) {
            cost += coef * sumAreaSq;

//            const double avg = sumArea / (double)cnt;
            const double avg = sumArea ;
            if (avg < lastAvg) {
                lastAvg = avg;
                lastIndex = (int)i;
                lastSumAreaSq = sumAreaSq;
            }
        }
    }

    // 如果所有 bin 都为空或都没有有效 rect，直接返回当前 cost
    if (lastIndex == -1) return cost;

    // 扣掉“平均面积最小的 bin”中所有矩形对应的项（与你原来的第二个 for 等价）
    cost -= coef * lastSumAreaSq;
    return cost;
}
void compute(Json::Value &result_list,string str) {

    vector<Rect> parts= toRectVector(Input(str));
    Rect binsize={0,0,Inputbinsize(str)[0],Inputbinsize(str)[1]};
    optimize(parts,binsize);

}
bool cmpArea(const Rect& a, const Rect& b) {
    return a.w*a.h > b.w*b.h;
}
bool cmpPerimeter(const Rect& a, const Rect& b)  {
    return (a.w+ a.h) > ( b.w+ b.h);
}

bool cmpLongSide(const Rect& a, const Rect& b) {
    return max(a.w, a.h) > max (b.w, b.h);
}

bool cmpShortSide(const Rect& a, const Rect& b)  {
    return min(a.w, a.h) > min( b.w, b.h);
}
int rand_with_linear_range(int iter, int max_iter, std::mt19937 &mt, int failiter) {
    if (max_iter <= 0) {
        std::uniform_int_distribution<int> dist_fallback(2, 4);
        return dist_fallback(mt);
    }

    iter = std::max(0, std::min(iter, max_iter));
//    int the_iter=std::min(20000,max_iter);
    int the_iter=max_iter/10;
    if(iter>the_iter){
        if((failiter*20/max_iter)%2==1){
            return 3;
        }
        return 2;
    }

    else{
        double t = static_cast<double>(iter) / static_cast<double>(the_iter); // t ∈ [0,1]

        // 线性插值端点

        double low_d  = 6.0 + (2.0 - 6.0) * t;   // 8 -> 2
        double high_d = 8.0 + (3.0 - 8.0) * t;  // 10 -> 4

        int low  = static_cast<int>(std::round(low_d));
        int high = static_cast<int>(std::round(high_d));

        // 保险起见，保证 low <= high，且落在 [2,10] 里
        if (low > high) std::swap(low, high);
        low  = std::max(2, low);
        high = std::min(10, high);

        std::uniform_int_distribution<int> dist(low, high);
        return dist(mt);
    }
}
double sumarea(vector<vector<vector<double>>> solution){
    double area=0;
    for (int i=0;i<solution.size();i++){
        for (int j=0;j<solution[i].size();j++){
            area+=   solution[i][j][2]*solution[i][j][3];
        }
    }
    return area;
}
static inline double binTotalArea(const Bin& b) {
    double sum = 0.0;
    for (const auto& r : b.rects) {
        if (r.w <= 0.0 || r.h <= 0.0) continue;  // 过滤无效块
        sum += r.w * r.h;
    }
    return sum;
}

// descending=true: 总面积大 -> 小（小的在最后）
// descending=false: 总面积小 -> 大（小的在最前）
void sortBinsByTotalArea(std::vector<Bin>& bins, bool descending = true) {
    std::sort(bins.begin(), bins.end(),
              [&](const Bin& a, const Bin& b) {
                  const double A = binTotalArea(a);
                  const double B = binTotalArea(b);
                  return descending ? (A > B) : (A < B);
              });
//    std::cout<<binTotalArea(bins[bins.size()-1])<<std::endl;
}
void optimize(vector<Rect> parts, Rect binsize) {

//    mt.seed(std::time(0u));
//    mt.seed(3);
    std::uniform_real_distribution<float> dt(0, 1);
    std::random_device rd;
    std::mt19937 mt(rd());
    bool success=true;
//    std::cout<<"optimize begin"<<std::endl;
    long long t=0;

    vector<Bin> bestsolution = regroup_2(parts, &success, binsize, 1,mt);
    vector<Bin> solution = bestsolution;
    int max_iter = 100000;

    // 累计时间（微秒）
    long long t_destroy1 = 0;
    long long t_destroy2 = 0;
    long long t_destroy3 = 0;
    long long t_repair3  = 0;
    long long t_repair4  = 0;
    long long t_accept   = 0;
    double current_cost=getcost(solution,binsize);
    double bestcost=current_cost;
    auto t_start = Clock::now();
    int failiter=0;
    for (int iter = 0; iter < max_iter; iter++) {

        // ---- 步骤 1：destroy_solution_1 + repair_solution4 + acceptance ----

        int destroynum = rand_with_linear_range(iter, max_iter, mt,failiter);
        vector<Bin> temp_solution = solution;
        auto t_d1_start = Clock::now();
        vector<Rect> destroy_set1 = destroy_solution_1(temp_solution, destroynum, mt);
        auto t_d1_end = Clock::now();
        t_destroy1 += std::chrono::duration_cast<microseconds>(t_d1_end - t_d1_start).count();

        // repair_solution4
        auto t_r4_start = Clock::now();
        bool flag = repair_solution4(temp_solution, destroy_set1, binsize, destroynum,mt);
        auto t_r4_end = Clock::now();
        t_repair4 += std::chrono::duration_cast<microseconds>(t_r4_end - t_r4_start).count();

        // acceptance
        auto t_acc_start = Clock::now();
        solution = acceptance(temp_solution, current_cost,solution, iter, max_iter, binsize, mt, flag,failiter);
        auto t_acc_end = Clock::now();
        t_accept += std::chrono::duration_cast<microseconds>(t_acc_end - t_acc_start).count();
        failiter+=1;
        if (bestcost < current_cost-0.0000001) {
            bestsolution = solution;
            sortBinsByTotalArea(bestsolution);
            failiter=0;
            bestcost=current_cost;
        }

        // ---- 步骤 2：每 10 轮做一次 destroy_solution_2 + repair_solution3 + acceptance ----
        if (iter % 5 == 1) {
            temp_solution = solution;

            auto t_d2_start = Clock::now();
            vector<Rect> destroy_set2 = destroy_solution_2(temp_solution);
            auto t_d2_end = Clock::now();
            t_destroy2 += std::chrono::duration_cast<microseconds>(t_d2_end - t_d2_start).count();

            auto t_r3_start = Clock::now();
            bool flag2 = repair_solution3(temp_solution, destroy_set2, binsize, 1,mt);
            auto t_r3_end = Clock::now();
            t_repair3 += std::chrono::duration_cast<microseconds>(t_r3_end - t_r3_start).count();

            auto t_acc2_start = Clock::now();
            solution = acceptance(temp_solution, current_cost,solution, iter, max_iter, binsize,  mt,flag2,failiter);
            auto t_acc2_end = Clock::now();
            t_accept += std::chrono::duration_cast<microseconds>(t_acc2_end - t_acc2_start).count();
            failiter+=1;
            if (bestcost < current_cost-0.000001) {
                bestsolution = solution;
                sortBinsByTotalArea(bestsolution);
                failiter=0;
                bestcost=current_cost;
            }
        }
        if (iter % 5 == 2) {
            temp_solution = solution;

            auto t_d2_start = Clock::now();
            vector<Rect> destroy_set2 = destroy_solution_5(temp_solution,mt);
            auto t_d2_end = Clock::now();
            t_destroy2 += std::chrono::duration_cast<microseconds>(t_d2_end - t_d2_start).count();

            auto t_r3_start = Clock::now();
            bool flag2 = repair_solution3(temp_solution, destroy_set2, binsize, 1,mt);
            auto t_r3_end = Clock::now();
            t_repair3 += std::chrono::duration_cast<microseconds>(t_r3_end - t_r3_start).count();

            auto t_acc2_start = Clock::now();
            solution = acceptance(temp_solution, current_cost,solution, iter, max_iter, binsize,  mt,flag2,failiter);
            auto t_acc2_end = Clock::now();
            t_accept += std::chrono::duration_cast<microseconds>(t_acc2_end - t_acc2_start).count();
            failiter+=1;
            if (bestcost < current_cost) {
                bestsolution = solution;
                sortBinsByTotalArea(bestsolution);
                failiter=0;
                bestcost=current_cost;
            }
        }
        if (iter % 5 == 3) {
            temp_solution = solution;

            auto t_d2_start = Clock::now();
            vector<Rect> destroy_set2 = destroy_solution_6(temp_solution);
            auto t_d2_end = Clock::now();
            t_destroy2 += std::chrono::duration_cast<microseconds>(t_d2_end - t_d2_start).count();

            auto t_r3_start = Clock::now();
            bool flag2 = repair_solution3(temp_solution, destroy_set2, binsize, 1,mt);
            auto t_r3_end = Clock::now();
            t_repair3 += std::chrono::duration_cast<microseconds>(t_r3_end - t_r3_start).count();

            auto t_acc2_start = Clock::now();
            solution = acceptance(temp_solution, current_cost,solution, iter, max_iter, binsize,  mt,flag2,failiter);
            auto t_acc2_end = Clock::now();
            t_accept += std::chrono::duration_cast<microseconds>(t_acc2_end - t_acc2_start).count();
            failiter+=1;
            if (bestcost < current_cost) {
                bestsolution = solution;
                sortBinsByTotalArea(bestsolution);
                bestcost=current_cost;
                failiter=0;
            }
        }
//        if (iter % 10 == 0) {
//            temp_solution = solution;
//
//            auto t_d2_start = Clock::now();
//            vector<Rect> destroy_set2 = destroy_solution_4(temp_solution);
//            auto t_d2_end = Clock::now();
//            t_destroy2 += std::chrono::duration_cast<microseconds>(t_d2_end - t_d2_start).count();
//
//            auto t_r3_start = Clock::now();
//            bool flag2 = repair_solution3(temp_solution, destroy_set2, binsize, 1,mt);
//            auto t_r3_end = Clock::now();
//            t_repair3 += std::chrono::duration_cast<microseconds>(t_r3_end - t_r3_start).count();
//
//            auto t_acc2_start = Clock::now();
//            solution = acceptance(temp_solution, current_cost,solution, iter, max_iter, binsize,  mt,flag2);
//            auto t_acc2_end = Clock::now();
//            t_accept += std::chrono::duration_cast<microseconds>(t_acc2_end - t_acc2_start).count();
//            if (bestcost < current_cost) {
//                bestsolution = solution;
//                sortBinsByTotalArea(bestsolution);
//                bestcost=current_cost;
//            }
//        }

        // ---- 步骤 3：每 10 轮、iter%10==1 时，再做一次 destroy3 + 两次 repair4 + 两次 acceptance ----
        if (iter % 300 == 299) {
//            temp_solution = solution;
//
//            auto t_d3_start = Clock::now();
//            vector<Rect> destroy_set3 = destroy_solution_3(temp_solution, 2);
//            auto t_d3_end = Clock::now();
//            t_destroy3 += std::chrono::duration_cast<microseconds>(t_d3_end - t_d3_start).count();
//
//            // 第一次 repair4 + acceptance
//            auto t_r4_1_start = Clock::now();
//            bool flag3 = repair_solution4(temp_solution, destroy_set3, binsize, 2,mt);
//            auto t_r4_1_end = Clock::now();
//            t_repair4 += std::chrono::duration_cast<microseconds>(t_r4_1_end - t_r4_1_start).count();
//
//            auto t_acc3_start = Clock::now();
//            solution = acceptance(temp_solution, current_cost,solution, iter, max_iter, binsize,mt, flag3);
//            auto t_acc3_end = Clock::now();
//            t_accept += std::chrono::duration_cast<microseconds>(t_acc3_end - t_acc3_start).count();

            // 第二次 repair4 + acceptance（你原来代码也有这两句，我保留）
            temp_solution = solution;

            auto t_d4_start = Clock::now();
            vector<Rect> destroy_set4 = destroy_solution_3(temp_solution, 2);
            auto t_d4_end = Clock::now();
            t_destroy3 += std::chrono::duration_cast<microseconds>(t_d4_end - t_d4_start).count();
            auto t_r4_2_start = Clock::now();
            bool flag4 = repair_solution4(temp_solution, destroy_set4, binsize, destroynum,mt);
            auto t_r4_2_end = Clock::now();
            t_repair4 += std::chrono::duration_cast<microseconds>(t_r4_2_end - t_r4_2_start).count();

            auto t_acc4_start = Clock::now();
            solution = acceptance(temp_solution, current_cost,solution, iter, max_iter, binsize,mt, flag4,failiter);
            auto t_acc4_end = Clock::now();
            t_accept += std::chrono::duration_cast<microseconds>(t_acc4_end - t_acc4_start).count();
            failiter+=1;
            if (bestcost < current_cost) {
                bestsolution = solution;
                sortBinsByTotalArea(bestsolution);
                bestcost=current_cost;
                failiter=0;
            }
        }


    }
    auto t_end = Clock::now();
    t += std::chrono::duration_cast<microseconds>(t_end - t_start).count();
    // 循环结束后，输出每个步骤的总耗时（秒）
//    std::cout << "Time destroy_solution_1: " << t_destroy1 / 1e6 << " s\n";
//    std::cout << "Time destroy_solution_2: " << t_destroy2 / 1e6 << " s\n";
//    std::cout << "Time destroy_solution_3: " << t_destroy3 / 1e6 << " s\n";
//    std::cout << "Time repair_solution3:   " << t_repair3  / 1e6 << " s\n";
//    std::cout << "Time repair_solution4:   " << t_repair4  / 1e6 << " s\n";
//    std::cout << "Time acceptance:         " << t_accept   / 1e6 << " s\n";
//    std::cout << "Time : " << t / 1e6 << " s\n";
    int binnum=0;
    bool res=validatePackingBool(bestsolution,parts,binsize,1e-5);
//    std::cout<<"feasibity="<<res<<std::endl;
    for(int i=0;i<bestsolution.size();i++)

{
        if (bestsolution[i].rects.size()>0) binnum+=1;
}
    fstream f;
    f.open("record.txt",ios::out|ios::app);
    f<<binnum<<std::endl;
    if (res){
        std::cout<<binnum<<std::endl;
    }
    else{std::cout<<99<<std::endl;}

    f.close();
//    }
}
std::vector<Bin> toBinVector(
        const std::vector<std::vector<std::vector<double>>>& vvv
) {
    std::vector<Bin> bins;
    bins.reserve(vvv.size());

    for (const auto& bin_vec : vvv) {
        Bin b;
        b.rects = toRectVector(bin_vec);  // 复用前面的转换函数
        bins.push_back(std::move(b));
    }
    return bins;
}
std::vector<Rect> toRectVector(const std::vector<std::vector<double>>& vvd) {
    std::vector<Rect> result;
    result.reserve(vvd.size());

    for (const auto& v : vvd) {
        if (v.size() < 4) {
            // 如果不满足 4 个元素，你可以选择：
            // 1) 跳过
            // 2) 或者抛异常
            // 这里先简单跳过
            continue;
        }
        Rect r;
        r.x = v[0];
        r.y = v[1];
        r.w = v[2];
        r.h = v[3];
        result.push_back(r);
    }

    return result;
}
std::vector<std::vector<double>> toVecVecDouble(const std::vector<Rect>& rects) {
    std::vector<std::vector<double>> result;
    result.reserve(rects.size());  // 预分配，稍微快一点

    for (const auto& r : rects) {
        // 每个 Rect 转成 [x, y, w, h]
        result.push_back({ r.x, r.y, r.w, r.h });
    }

    return result;
}
static inline double rectArea(const Rect& r) {
    return r.w * r.h;
}

static inline bool insideBin(const Rect& r, const Rect& binSize, double eps) {
    return (r.x >= -eps) &&
           (r.y >= -eps) &&
           (r.w >  eps) &&
           (r.h >  eps) &&
           (r.x + r.w <= binSize.w + eps) &&
           (r.y + r.h <= binSize.h + eps);
}

static inline bool overlapsPositiveArea(const Rect& a, const Rect& b, double eps) {
    const double ax2 = a.x + a.w, ay2 = a.y + a.h;
    const double bx2 = b.x + b.w, by2 = b.y + b.h;
    const double ox = std::min(ax2, bx2) - std::max(a.x, b.x);
    const double oy = std::min(ay2, by2) - std::max(a.y, b.y);
    return (ox > eps) && (oy > eps);
}

bool validatePackingBool(
        const std::vector<Bin>& bins,
        const std::vector<Rect>& inputParts,
        const Rect& binSize,
        double eps = 1e-9
) {
    // 输入统计
    const size_t inCnt = inputParts.size();
    double inArea = 0.0;
    for (const auto& p : inputParts) {
        inArea += p.w * p.h;
    }

    // 输出统计 + 可行性检查
    size_t outCnt = 0;
    double outArea = 0.0;

    for (const auto& b : bins) {
        const auto& rs = b.rects;

        // 边界 + 面积统计
        for (const auto& r : rs) {
            if (!insideBin(r, binSize, eps)) return false;
            ++outCnt;
            outArea += r.w * r.h;
        }

        // 重叠检查
        for (size_t i = 0; i < rs.size(); ++i) {
            for (size_t j = i + 1; j < rs.size(); ++j) {
                if (overlapsPositiveArea(rs[i], rs[j], eps)) return false;
            }
        }
    }

    // 个数一致
    if (outCnt != inCnt) return false;

    // 面积一致（相对容差）
    const double diff = std::fabs(outArea - inArea);
    const double scale = std::max(1.0, std::fabs(inArea));
    if (diff > eps * scale) return false;

    return true;
}