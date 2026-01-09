//
// Created by raiser on 2021/12/13.
//

#ifndef AMOPT_CARVINGMACHINE_H
#define AMOPT_CARVINGMACHINE_H

#include "amopt_pack_base.h"
#include <unordered_map>
//#include "Instance.h"
//#include "visitor.h"
#include <numeric>
#include "ctime"
#include <fstream>
//#include "marching/Graph.h"
//#include "marching/Matching.h"
//#include "marching/utils.h"
#include <random>
using namespace amopt;
using std::ofstream;
using std::make_shared;

struct Rect {
    double x, y, w, h;
};
struct Bin{
    vector<Rect> rects;
};
struct PairInfo {
    std::vector<Rect> rectangles; // 这里先用 int 表示 rectangle 的 id；如果你有 Rectangle 类型可替换
    int position;              // (x,y) 或者你需要的坐标含义
};
vector<vector<double>> Input(
        string str
);

void compute(Json::Value &result_list,string str);
Error computeBRKGA(Json::Value &result_list);




void combination();

void firstpack();
void optimize(vector<Rect> parts, Rect binsize);
void optimizeBRKGA();


void output(Json::Value &result_list);
std::vector<std::vector<float>> initial(int p, int num);
std::vector<float> Decode(std::vector<std::vector<float>> &population);
float Decode_(std::vector<std::vector<float>> population);
std::vector<std::vector<float>> Selectpe(int pe, std::vector<std::vector<float>> population);
std::vector<std::vector<float>> CrossOver(std::vector<std::vector<float>> populatione, std::vector<std::vector<float>> population, int p_, float rou);
std::vector<std::vector<float>> Mutiation(int pm, int num);
float calfit(std::vector<float> population);
float calfit_(std::vector<float> population);
std::vector<Bin> maxRects(
        std::vector<Rect> parts,     // 按值：内部要 shuffle + erase
        double width,
        double height,
        int binNum,
        float* score1,
        float* score2,
        std::mt19937& mt             // 关键：按引用
);
Rect FindPositionForNewNodeBestShortSideFit_(double width, double height,
                                             double &bestShortSideFit,
                                             double &bestLongSideFit,std::vector<Rect> freeRectangles,std::mt19937 mt);
Rect FindPositionForNewNodeBestShortSideFit(
        double width,
        double height,
        const std::vector<Rect>& freeRectangles,
        std::mt19937& mt
);
std::vector<Bin> regroup_1(
        std::vector<Rect> destroy_set,   // 按值：内部可直接传给 maxRects
        bool* success_flag,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
);
std::vector<Bin> regroup_2(
        std::vector<Rect> destroy_set,   // 按值：内部要 sort
        bool* success_flag,
        const Rect& binSize,             // 只用 binSize.w / binSize.h
        int destroy_num,
        std::mt19937& mt                 // 关键：按引用
);
std::vector<Bin> regroup_3(
        std::vector<Bin> solution,                 // 按值：返回新解（与原函数一致）
        const std::vector<Rect>& destroy_set,       // const&：避免拷贝
        bool* success_flag,
        const Rect& binSize,
        int /*destroy_num*/,                        // 原函数里没用到，保留占位
        std::mt19937& mt                            // 关键：按引用
);
std::vector<Bin> regroup_4(
        std::vector<Rect> destroy_set,   // 按值：我们要 sort/shuffle
        bool* success_flag,
        const Rect& binSize,             // 用 Rect 承载 w,h（x,y 忽略）
        int destroy_num,
        std::mt19937& mt
);
bool SplitFreeNode(Rect freeNode, const Rect &usedNode,std::vector<Rect> &freeRectangles);
void PruneFreeList(std::vector<Rect> &freeRectangles);
bool IsContainedIn(const Rect &a, const Rect &b);
double Insertscore(
        double width,
        double height,
        const std::vector<Rect>& freeRectangles,  // 关键：const&
        std::mt19937& mt                           // 关键：&
);
Rect Insert(
        double width,
        double height,
        std::vector<Rect>& freeRectangles,
        std::mt19937& mt
);
void repair_solution1(
        std::vector<Bin>& solution,
        const std::vector<Rect>& destroy_set,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
);
bool repair_solution3(
        std::vector<Bin>& solution,
        const std::vector<Rect>& destroy_set,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
);
bool repair_solution4(
        std::vector<Bin>& solution,
        const std::vector<Rect>& destroy_set,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
);
std::vector<Rect> insertparts_1(
        std::vector<Rect>& bin,
        const Rect& part,
        const Rect& binSize,
        std::mt19937& mt,
        bool *flag
);
void repair_solution2(
        std::vector<Bin>& solution,
        const std::vector<Rect>& destroy_set,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
);
bool cmpArea(const Rect& a, const Rect& b) ;
bool cmpPerimeter(const Rect& a, const Rect& b) ;
std::vector<Rect> destroy_solution_1(std::vector<Bin> &sol,
                                     int binnum,
                                     std::mt19937 &mt);
bool cmpLongSide(const Rect& a, const Rect& b) ;

bool cmpShortSide(const Rect& a, const Rect& b) ;
std::vector<Rect> destroy_solution_3(std::vector<Bin>& sol, int binnum);
double getcost(const std::vector<Bin>& solution, const Rect& binSize);
int rand_with_linear_range(int iter, int max_iter, std::mt19937 &mt,int failiter);
//vector<vector<vector<double>>> regroup_5(vector<vector<double>> destroy_set, bool *success_flag, vector<double> binsize, int destroy_num,std::mt19937 mt);
bool check_bins_no_overlap(const vector<vector<vector<double>>>& bins);
bool rectsOverlap(const vector<double>& a, const vector<double>& b);
double sumarea(vector<vector<vector<double>>> solution);
//bool repair_solution5(vector<vector<vector<double>>> &solution, vector<vector<double>> destroy_set,vector<double> binsize,int destroy_num);
std::vector<std::vector<double>> toVecVecDouble(const std::vector<Rect>& rects) ;
std::vector<Rect> toRectVector(const std::vector<std::vector<double>>& vvd);
std::vector<Bin> toBinVector(
        const std::vector<std::vector<std::vector<double>>>& vvv);
bool validatePackingBool(
        const std::vector<Bin>& bins,
        const std::vector<Rect>& inputParts,
        const Rect& binSize,
        double eps
);
std::vector<Rect> destroy_solution_4(std::vector<Bin>& sol);
std::vector<Rect> destroy_solution_5(std::vector<Bin>& sol, std::mt19937& mt);
std::vector<Rect> destroy_solution_6(std::vector<Bin>& sol);
double sumarea(vector<vector<vector<double>>> solution);
static inline double binTotalArea(const Bin& b);
void sortBinsByTotalArea(std::vector<Bin>& bins, bool descending);
static inline double sumArea(const std::vector<Rect>& rs);
std::vector<Bin> regroup_3_(
        std::vector<Bin> solution,                 // 按值：返回新解（与原函数一致）
        const std::vector<Rect>& destroy_set,       // const&：避免拷贝
        bool* success_flag,
        const Rect& binSize,
        int /*destroy_num*/,                        // 原函数里没用到，保留占位
        std::mt19937& mt                            // 关键：按引用
);
bool repair_solution5(
        std::vector<Bin>& solution,
        const std::vector<Rect>& destroy_set,
        const Rect& binSize,
        int destroy_num,
        std::mt19937& mt
);
std::vector<Rect> Insert_by_LLM(double width, double height, vector<Rect> &freeRectangles, std::mt19937 &mt, PairInfo info,int method);
Rect FindPosition_by_LLM(double width,double height,std::vector<Rect> freeRectangles,std::mt19937& mt,int method);
std::unordered_map<int, PairInfo> Sort_by_LLM(std::vector<Rect>& parts_sorted, const Rect& binSize, std::mt19937& mt,int method);
#endif //AMOPT_CARVINGMACHINE_H

