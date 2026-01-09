//
// Created by dxy on 2021/12/16.
//

#ifndef PACKING_SOLUTION_AMOPT_PACK_BASE_H
#define PACKING_SOLUTION_AMOPT_PACK_BASE_H

#include <iostream>
#include<vector>
#include<set>
#include <string>
#include<json/json.h>
#include "iso646.h"
#include "algorithm_error.h"
//#include "heuristic/Rect.h"
#include "algorithm_error.h"
//#include "heuristic/Rect.h"
#include <memory>
#include "map"

using std::string;
using std::shared_ptr;
using std::vector;
using std::map;





        class Order {
        public:
            Order();

            ~Order();
        };

        class Configure {
            double pack_divide;
            double cut_divide;
            double minimum_width;
            vector<int> angles;
            int yongliao_score;
            int zhengkuai_score;
            int bianjiao_score;
            int waiqiege_score;
            int neiqiege_score;
            int chazuo_score;
            std::string pack_type;
            int zhipu_n;
            bool wall_flag;
            bool group_flag;
            double remnantsarea;
            double remnantsspace;
            double remnantsheight;
            double remnantswidth;
            bool surplusDir;
        public:
            int timeLimit;

            Configure(double remnantsarea_,
                      double remnantsspace_,
                      double remnantsheight_,
                      double remnantswidth_,
                      bool surplusDir_,
                      int timeLimit_) {
                remnantsarea = remnantsarea_;
                remnantsspace = remnantsspace_;
                remnantsheight = remnantsheight_;
                remnantswidth = remnantswidth_;
                surplusDir = surplusDir_;
                timeLimit = timeLimit_;
            }

        public:
            void show() {
                std::cout << timeLimit;
            }

            int SurplusDir() {
                if (surplusDir) {
                    return 2;
                } else { return 1; }
            }

            Configure();

            //~configure();

        };

        class Output {
            int startpoint;
            int angle;
            int yongliao;
            int zhengkuai;
            int bianjiao;
            int waiqiege;
            int neiqiege;
            int chazuo;
            int calculatetime;
            std::vector<std::vector<std::vector<double>>>
                    zheng;
            std::vector<std::vector<std::vector<double>>>
                    bian;
            std::vector<std::vector<std::vector<double>>>
                    cut;
        public:
            Output();

            ~Output();
        };

#endif //PACKING_SOLUTION_AMOPT_PACK_BASE_H
