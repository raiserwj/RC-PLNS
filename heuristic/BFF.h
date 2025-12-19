#include<iostream>
#include <vector>
#include "GuillotineBinPack.h"
#include <iostream>
#include <algorithm>

using namespace rbp;
class BFF
{

public:
    std::vector<std::vector<rbp::Rect > > pack(std::vector<std::vector<double> > rects_, double bin_width,double bin_height,int nBin,int rotate);
    int num_bins;
    int num_qiege;
};