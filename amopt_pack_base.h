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
#include "heuristic/Rect.h"
#include "heuristic/polygons.h"
#include <memory>
#include "map"
using std::string;
using std::shared_ptr;
using std::vector;
using std::map;

using namespace rbp;
using namespace amopt;
namespace amopt {
    namespace amopt_pack {
        class CombineTransform {
        public:
            int rotation;
            double xShift;
            double yShift;
        };

        class Part {
        public:
            Rect rect;
            Rect temp_rect;
            double height;
            double width;
            string id1;
            int id2;
            int type;

            vector <vector<float>> points;
            vector <vector<vector < float>>>
            pointsCombine;
            vector<float> boundingBox;
            CombineTransform transform;
            vector <shared_ptr<amopt_pack::Part>> combineParts;

            bool backfrontpriority;
            bool smallitem;
            bool rotate;
            int rotate_degree;
            vector<int> rotate_degrees;
            double size;
            float area;

            void
            Input(double height_, double width_, float area_, vector <vector<float>> points_, vector<float> boundingBox_, vector<int> rotate_degrees_, string id1_,
                  bool backFrontPriority_, bool smallItem_, bool rotate_, int rotate_degree_, int type_) {
                height = height_;
                width = width_;
                size = width_ * height_;
                area = area_;
                points={};
                for(int i=0;i<points_.size();i++){
                    vector<float> point=points_[i];
                    points.push_back(point);
                }
                pointsCombine = {points_};
                boundingBox = boundingBox_;
                rotate_degrees = rotate_degrees_;
                id1 = id1_;
                backfrontpriority = backFrontPriority_;

                smallitem = smallItem_;
                rotate = rotate_;
                rotate_degree = rotate_degree_;
                type = type_;
            }


            void show() {

                for (int i = 0; i < points.size(); i++) {
                    for (int j = 0; j < points[i].size(); j++) {
                        std::cout << points[i][j] << ",";
                    }
                    std::cout << "\n";
                }
                std::cout << "\n";
            }

            bool Showbackfrontpriority() {

                return backfrontpriority;
            }

            //part();

        };

        class Bin {
            std::vector<std::vector<std::vector<int>>>
                    floor;
            std::vector<std::vector<std::vector<double>>> result;
            double size;
        public:
            double averagearea;
            double fillarea;
            vector <shared_ptr<amopt_pack::Part>> parts;
            double height;
            double width;
            bool backornot;
            bool big_part;
            double unity;
            double fitness;
            double overlap;
            vector<vector<vector<vector<vector<vector<float>>>>>> nfp;//(0,0)to(0,0)

            void updatenfp(vector<int> partsindex,vector<vector<vector<vector<vector<vector<float>>>>>> nfp_){
//                nfp={};
//                int flag1=0,flag2=0;
//                int n=nfp_.size();
//                int i=0;
//                while(flag1<n){
//                    if(std::find(partsindex.begin(),partsindex.end(),flag1)!=partsindex.end())
//                    {
//                        flag2=0;
//                        nfp.push_back({});
//                        while(flag2<n){
//                            if(std::find(partsindex.begin(),partsindex.end(),flag2)!=partsindex.end()){
//                                nfp[i].push_back(nfp_[flag1][flag2]);
//                            }
//                            flag2++;
//                        }
//                        i++;
//                    }
//                    flag1++;
//                }
                nfp=nfp_;
            }
            double caloverlap(){
                overlap=0;
                for(int i=0;i<parts.size();i++){
                    for(int j=0;j<parts.size();j++){
                        if(i==j){
                            continue;
                        }
                        else{
                            vector<vector<float>> nfp_={};
                            for(int k=0;k<nfp[parts[i]->id2][parts[j]->id2][parts[i]->rotate_degree][parts[j]->rotate_degree].size();k++){
                                nfp_.push_back(nfp[parts[i]->id2][parts[j]->id2][parts[i]->rotate_degree][parts[j]->rotate_degree][k]);
                            }
                            for(int k=0;k<nfp_.size();k++){
                                nfp_[k][0]+=parts[j]->points[0][0];
                                nfp_[k][1]+=parts[j]->points[0][1];
                            }
                            if(pointInPolygon(nfp_,parts[i]->points[0])){
                                overlap+=Distance(nfp_,parts[i]->points[0])*Distance(nfp_,parts[i]->points[0]);
                            }
                        }
                    }
                    float x_max=0;
                    for(int j=0;j<parts[i]->points.size();j++){
                        if(parts[i]->points[j][0]>x_max){
                            x_max=parts[i]->points[j][0];
                        }
                    }
                    if(x_max>width){
                        overlap+=(x_max-width)*(x_max-width);
                    }
                    float y_max=0;
                    for(int j=0;j<parts[i]->points.size();j++){
                        if(parts[i]->points[j][1]>y_max){
                            y_max=parts[i]->points[j][1];
                        }
                    }
                    if(y_max>height){
                        overlap+=(y_max-height)*(y_max-height);
                    }
                    float x_min=0;
                    for(int j=0;j<parts[i]->points.size();j++){
                        if(parts[i]->points[j][0]<x_min){
                            x_min=parts[i]->points[j][0];
                        }
                    }
                    if(x_min<0){
                        overlap+=(x_min)*(x_min);
                    }
                    float y_min=0;
                    for(int j=0;j<parts[i]->points.size();j++){
                        if(parts[i]->points[j][1]<y_min){
                            y_min=parts[i]->points[j][1];
                        }
                    }
                    if(y_min<0){
                        overlap+=(y_min)*(y_min);
                    }
                }
                return overlap;
            }
            Bin(double width_, double height_, bool backornot_) {
                width = width_;
                height = height_;
                backornot = backornot_;
                size = width * height;
            };

            Bin(double width_, double height_, bool backornot_, vector <shared_ptr<amopt_pack::Part>> parts_) {
                width = width_;
                height = height_;
                backornot = backornot_;
                parts = parts_;
                size = width * height;
            };

            void update(vector <shared_ptr<amopt_pack::Part>> parts_) {
                parts = parts_;
            }

            void update_area_sum() {
                fillarea = 0;
                for (int i = 0; i < parts.size(); i++) {
                    fillarea += parts[i]->area;
                }
                if (fillarea==0){
                    return;
                }
                averagearea = fillarea / parts.size();
                unity = fillarea / size;
                for (int i = 0; i < parts.size(); i++) {
                    if (not parts[i]->smallitem) {
                        big_part = true;
                        break;
                    }
                }
                fitness = (backornot ? 2 : 1) * (4 * unity * unity - 10) + 3 * (big_part ? 1 : 0);
//                fitness = (backornot ? 2 : 1) * (4 * unity * unity - 50) + 3 * ((not parts[0]->smallitem) ? 1 : 0);
            }

            void updata_area_add(double area_, Part *part) {
                fillarea += area_;
                averagearea = fillarea / parts.size();
                unity = fillarea / size;
                big_part = big_part | (!part->smallitem);
            }

            void update_fitness() {
                fitness = (1 + backornot) * 4 * unity * unity + 3 * big_part;
            }
        };


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

    }
}

typedef vector<shared_ptr<amopt_pack::Bin>> Bin_Ptr_V;
typedef vector<shared_ptr<amopt_pack::Part>> Part_Ptr_V;
typedef shared_ptr<amopt_pack::Bin> Bin_Ptr;
typedef shared_ptr<amopt_pack::Part> Part_Ptr;

#endif //PACKING_SOLUTION_AMOPT_PACK_BASE_H
