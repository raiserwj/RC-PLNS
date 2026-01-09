#ifndef PACKING_SOLUTION_CBFF_H
#define PACKING_SOLUTION_CBFF_H

#include<iostream>
#include <vector>
#include "MaxRectsBinPack.h"
#include <iostream>
#include <algorithm>
#include <string>
#include <random>
#include "amopt_pack_base.h"

using namespace rbp;
using namespace amopt;

Bin_Ptr_V regroup(Bin_Ptr_V bins_, bool *success_flag, std::mt19937 &mt, bool return_flag = false, bool junyun = false, bool sort = false);

Bin_Ptr_V regroup_minus(Bin_Ptr_V bins_, bool *success_flag, bool return_flag = false, double ratio = 0);
//vector<amopt_pack::Bin>
//regroup_new(std::vector<amopt_pack::Bin> bins_, bool success_flag);

Part_Ptr_V insert(Bin_Ptr bin, Part_Ptr part, std::mt19937 &mt, bool initial_flag = false, bool return_flag = false, bool sort = false);

bool output_pack(Bin_Ptr bin);

bool output_pack_last(Bin_Ptr bin, int index, bool return_flag);

bool output_pack_last_new(Bin_Ptr bin, int index, bool return_flag);

void last_insert(vector<Bin_Ptr> bins, int index, std::mt19937 &mt);

vector<Part_Ptr_V> maxRects(Part_Ptr_V parts, float width, float height, int binNum, std::mt19937 &mt, float *score1, float *score2, bool junyun);

vector<Part_Ptr_V> maxRectsSorts(Part_Ptr_V parts, float width, float height, int binNum, std::mt19937 &mt, float *score1, float *score2);

#endif //PACKING_SOLUTION_CBFF_H