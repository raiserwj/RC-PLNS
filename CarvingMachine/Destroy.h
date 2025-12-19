//
// Created by dxy on 2021/12/13.
//

#ifndef PACKING_SOLUTION_DESTROY_H
#define PACKING_SOLUTION_DESTROY_H

#include "CarvingMachine.h"
#include "amopt_pack_base.h"
#include "PackSolution.h"
#include "PALNS/DestroyMethod.h"
#include "numeric"
#include "PALNS/PALNS.h"
#include "algorithm"
#include "random"

using std::max;
using std::min;
using std::make_unique;

using namespace mlpalns;
using namespace amopt_pack;
namespace amopt {
    namespace amopt_pack {
        struct DestroyReconstruct : public mlpalns::DestroyMethod<PackSolution> {
        public:
            DestroyReconstruct() = default;

            DestroyReconstruct(int type_, float frac_, int min_num_) {
                frac = frac_;
                type = type_;
                min_num = min_num_;
            }

            void set(float frac_, int min_num_) {
                frac = frac_;
                min_num = min_num_;
            }

            std::unique_ptr<mlpalns::DestroyMethod <
                            PackSolution>> clone() const override { return make_unique<DestroyReconstruct>(*this); }

            void destroy_solution(PackSolution &sol, std::mt19937 &mt) override {
                sol.type = type;
                std::uniform_real_distribution<float> dt(0, 1);
                if (type == 100) {
//                    std::shuffle(sol.bins_back.begin(), sol.bins_back.end(), mt);
                    int num = max(int(frac * (sol.bins_back.size() + sol.bins_normal.size())), min_num);
                    int size = sol.bins_back.size();
                    vector<int> indexs;
                    vector<int> out;
                    for (int i = size - 1; i >= 0; i--) {
                        if (sol.bins_back[i]->unity < sol.unity_rate3_back) {
                            indexs.push_back(i);
                        }
                    }
                    if (indexs.size() == 0) return;
                    if (indexs.size() < num) num = indexs.size();
                    std::sample(indexs.begin(), indexs.end(), std::back_inserter(out), num, mt);
                    for (int i = 0; i < out.size(); i++) {
                        sol.destroy_set.push_back(sol.bins_back[out[i]]);
                        sol.bins_back.erase(sol.bins_back.begin() + out[i]);
                    }
                } else if (type == 101) {
                    int num = max(int(frac * (sol.bins_back.size() + sol.bins_normal.size())), min_num);
                    int size = sol.bins_normal.size();
                    vector<int> indexs;
                    vector<int> out;
                    for (int i = size - 1; i >= 0; i--) {
                        if (sol.bins_normal[i]->unity < sol.unity_rate3_normal) {
                            indexs.push_back(i);
                        }
                    }
                    if (indexs.size() == 0) return;
                    if (indexs.size() < num) num = indexs.size();
                    std::sample(indexs.begin(), indexs.end(), std::back_inserter(out), num, mt);
                    for (int i = 0; i < out.size(); i++) {
                        sol.destroy_set.push_back(sol.bins_normal[out[i]]);
                        sol.bins_normal.erase(sol.bins_normal.begin() + out[i]);
                    }
                } else if (type == 102) {
                    std::shuffle(sol.bins_back.begin(), sol.bins_back.end(), mt);
                    sol.destroy_set.insert(sol.destroy_set.begin(), sol.bins_back.begin(),
                                           sol.bins_back.begin() + min_num);
                    sol.bins_back.erase(sol.bins_back.begin(),
                                        sol.bins_back.begin() + min_num);
                } else if (type == 103) {
                    std::shuffle(sol.bins_normal.begin(), sol.bins_normal.end(), mt);
                    sol.destroy_set.insert(sol.destroy_set.begin(), sol.bins_normal.begin(),
                                           sol.bins_normal.begin() + min_num);
                    sol.bins_normal.erase(sol.bins_normal.begin(),
                                          sol.bins_normal.begin() + min_num);
                }
            }

            ~DestroyReconstruct() override = default;

        private:
            float frac;
            int type;
            int min_num;
        };

        struct DestroyInsertAverage : public mlpalns::DestroyMethod<PackSolution> {
        public:
            DestroyInsertAverage() = default;

            DestroyInsertAverage(int type_) { type = type_; }

            std::unique_ptr<mlpalns::DestroyMethod <
                            PackSolution>> clone() const override { return make_unique<DestroyInsertAverage>(*this); }

            void destroy_solution(PackSolution &sol, std::mt19937 &mt) override {
                sol.type = type;
                if (type == 200) {
                    int last_average = sol.calculate_last_average(true, true);
                    if (last_average == -1) return;
                    sol.destroy_set.push_back(sol.bins_back[last_average]);
                    sol.bins_back.erase(sol.bins_back.begin() + last_average);
                } else if (type == 201) {
                    int last2_average = sol.calculate_last_average(true, false);
                    if (last2_average == -1) return;
                    sol.destroy_set.push_back(sol.bins_back[last2_average]);
                    sol.bins_back.erase(sol.bins_back.begin() + last2_average);
                } else if (type == 202) {
                    int last_average = sol.calculate_last_average(false, true);
                    if (last_average == -1) return;
                    sol.destroy_set.push_back(sol.bins_normal[last_average]);
                    sol.bins_normal.erase(sol.bins_normal.begin() + last_average);
                } else if (type == 203) {
                    int last2_average = sol.calculate_last_average(false, false);
                    if (last2_average == -1) return;
                    sol.destroy_set.push_back(sol.bins_normal[last2_average]);
                    sol.bins_normal.erase(sol.bins_normal.begin() + last2_average);
                }
            }

            ~DestroyInsertAverage() override = default;

        private:
            int type;
        };

        struct DestroyInsertFill : public mlpalns::DestroyMethod<PackSolution> {
        public:
            DestroyInsertFill() = default;

            DestroyInsertFill(int type_) { type = type_; }

            std::unique_ptr<mlpalns::DestroyMethod <
                            PackSolution>> clone() const override { return make_unique<DestroyInsertFill>(*this); }

            void destroy_solution(PackSolution &sol, std::mt19937 &mt) override {
                sol.type = type;
                if (type == 300) {
                    int last_fill = sol.calculate_last_fill(true, true);
                    sol.destroy_set.push_back(sol.bins_back[last_fill]);
                    sol.bins_back.erase(sol.bins_back.begin() + last_fill);
                } else if (type == 301) {
                    int last2_fill = sol.calculate_last_fill(true, false);
                    sol.destroy_set.push_back(sol.bins_back[last2_fill]);
                    sol.bins_back.erase(sol.bins_back.begin() + last2_fill);
                } else if (type == 302) {
                    int random_fill = sol.insert_fill_rate1(true, mt);
                    if (random_fill == -1) return;
                    sol.destroy_set.push_back(sol.bins_back[random_fill]);
                    sol.bins_back.erase(sol.bins_back.begin() + random_fill);
                } else if (type == 303) {
                    int last_fill = sol.calculate_last_fill(false, true);
                    sol.destroy_set.push_back(sol.bins_normal[last_fill]);
                    sol.bins_normal.erase(sol.bins_normal.begin() + last_fill);
                } else if (type == 304) {
                    int last2_fill = sol.calculate_last_fill(false, false);
                    sol.destroy_set.push_back(sol.bins_normal[last2_fill]);
                    sol.bins_normal.erase(sol.bins_normal.begin() + last2_fill);
                } else if (type == 305) {
                    int random_fill = sol.insert_fill_rate1(false, mt);
                    if (random_fill == -1) return;
                    sol.destroy_set.push_back(sol.bins_normal[random_fill]);
                    sol.bins_normal.erase(sol.bins_normal.begin() + random_fill);
                }
            }

            ~DestroyInsertFill() override = default;

        private:
            int type;
        };

        struct DestroyReconstruct2 : public mlpalns::DestroyMethod<PackSolution> {
        public:
            DestroyReconstruct2() = default;

            DestroyReconstruct2(int type_) {
                type = type_;
            }

            std::unique_ptr<mlpalns::DestroyMethod <
                            PackSolution>> clone() const override { return make_unique<DestroyReconstruct2>(*this); }

            void destroy_solution(PackSolution &sol, std::mt19937 &mt) override {
                sol.type = type;
                if (type == 400) {
                    vector<int> indexs = sol.group_fill_rate2(true, mt);
                    if (indexs.size() != 2)return;
                    sol.destroy_set.push_back(sol.bins_back[indexs[0]]);
                    sol.destroy_set.push_back(sol.bins_back[indexs[1]]);
                    sol.bins_back.erase(sol.bins_back.begin() + max(indexs[0], indexs[1]));
                    sol.bins_back.erase(sol.bins_back.begin() + min(indexs[0], indexs[1]));
                } else if (type == 401) {
                    vector<int> indexs = sol.group_fill_rate2(false, mt);
                    if (indexs.size() != 2)return;
                    sol.destroy_set.push_back(sol.bins_normal[indexs[0]]);
                    sol.destroy_set.push_back(sol.bins_normal[indexs[1]]);
                    sol.bins_normal.erase(sol.bins_normal.begin() + max(indexs[0], indexs[1]));
                    sol.bins_normal.erase(sol.bins_normal.begin() + min(indexs[0], indexs[1]));
                } else if (type == 402) {
                    vector<int> indexs = sol.group_fill_last();
                    if (indexs.size() != 2)return;
                    sol.destroy_set.push_back(sol.bins_normal[indexs[0]]);
                    sol.destroy_set.push_back(sol.bins_normal[indexs[1]]);
                    sol.bins_normal.erase(sol.bins_normal.begin() + max(indexs[0], indexs[1]));
                    sol.bins_normal.erase(sol.bins_normal.begin() + min(indexs[0], indexs[1]));
                }
            }

            ~DestroyReconstruct2() override = default;

        private:
            int type;
        };

        struct DestroyBack : public mlpalns::DestroyMethod<PackSolution> {
        public:
            DestroyBack() = default;

            DestroyBack(int type_) {
                type = type_;
            }

            std::unique_ptr<mlpalns::DestroyMethod <
                            PackSolution>> clone() const override { return make_unique<DestroyBack>(*this); }

            void destroy_solution(PackSolution &sol, std::mt19937 &mt) override {
                sol.type = type;
                int index;
                if (type == 500) {
                    index = sol.calculate_last_average(false, true);
                } else {
                    index = sol.calculate_last_average(false, false);
                }
                if (index == -1) return;
                sol.destroy_set.push_back(sol.bins_normal[index]);
                sol.bins_normal.erase(sol.bins_normal.begin() + index);
            }
            ~DestroyBack() override = default;

        private:
            int type;
        };
    };
};


#endif //PACKING_SOLUTION_DESTROY_H
