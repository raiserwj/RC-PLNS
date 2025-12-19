//
// Created by raiser on 2021/12/13.
//

#ifndef AMOPT_CARVINGMACHINE_H
#define AMOPT_CARVINGMACHINE_H

#include "amopt_pack_base.h"
#include "PackSolution.h"
#include "CBFF.h"
#include "Destroy.h"
#include "Repair.h"
#include "Instance.h"
#include "visitor.h"
#include <numeric>
#include "ctime"
#include <fstream>
#include "marching/Graph.h"
#include "marching/Matching.h"
#include "marching/utils.h"

using namespace amopt;
using std::ofstream;
using std::make_shared;

namespace amopt {
    namespace amopt_pack {
        class CarvingMachine {
            std::vector<Part_Ptr> squareparts;
            std::vector<Part_Ptr> notsquareparts;
            std::vector<Part_Ptr> parts;
            std::vector<Bin_Ptr> bins;
            std::vector<Part_Ptr> backparts;
            std::vector<Part_Ptr> notbackparts;
            amopt_pack::Configure configure;
            PackSolution solution;
            PackSolution bestSolution;
            struct timeval start;
            struct timeval end;
            std::mt19937 mt;
            std::vector<vector<double>> fitlist;
            int numr;
            int numi;
        public:
            CarvingMachine();

            ~CarvingMachine();

            Error Input(
                    string str
            );

            Error compute(Json::Value &result_list);
            Error computeBRKGA(Json::Value &result_list);

            void numitems() {
                std::cout << "normal parts:" << squareparts.size() << "\n" << "abnormalparts:" << notsquareparts.size()
                          << "\n" << "bin size:" << bins.size() << "\n" << std::endl;
            }

        private:
            float ratio;

            void combination();

            void firstpack();

            void optimize();
            void optimizeBRKGA();

            void
            move(mlpalns::DestroyMethod<PackSolution> &destroy, mlpalns::RepairMethod<PackSolution> &repair,
                 std::mt19937 &mt,float tem,int &num);

            void output(Json::Value &result_list);
            std::vector<std::vector<float>> initial(int p, int num);
            std::vector<float> Decode(std::vector<std::vector<float>> &population);
            float Decode_(std::vector<std::vector<float>> population);
            std::vector<std::vector<float>> Selectpe(int pe, std::vector<std::vector<float>> population);
            std::vector<std::vector<float>> CrossOver(std::vector<std::vector<float>> populatione, std::vector<std::vector<float>> population, int p_, float rou);
            std::vector<std::vector<float>> Mutiation(int pm, int num);
            float calfit(std::vector<float> population);
            float calfit_(std::vector<float> population);
        };
    }
}


#endif //AMOPT_CARVINGMACHINE_H
