//
// Created by dxy on 2021/12/13.
//

#ifndef PACKING_SOLUTION_PACKSOLUTION_H
#define PACKING_SOLUTION_PACKSOLUTION_H

#include <random>
#include "amopt_pack_base.h"
#include "set"
#include <algorithm>

using std::min;
namespace amopt {
    namespace amopt_pack {
        class PackSolution {
        public:
            double cost = 0;
            Bin_Ptr_V bins_normal;
            Bin_Ptr_V bins_back;
            Bin_Ptr_V destroy_set;
            Part_Ptr_V parts_influenced;
            int type;

            double getcost() {
                bool back_flag = (bins_normal.size() == 0) ? true : false;
                int last_index = calculate_last_average(back_flag, true);
                double fitness = 0;
                for (int i = 0; i < bins_normal.size(); i++) {
                    fitness += bins_normal[i]->fitness;
                }
                for (int i = 0; i < bins_back.size(); i++) {
                    fitness += bins_back[i]->fitness;
                }
                if (last_index != -1) {
                    if (back_flag == false)fitness -= bins_normal[last_index]->fitness;
                    else fitness -= bins_back[last_index]->fitness;
                }
                return -(fitness);
            }

            double getcost_last() {
                double fitness = 0;
                for (int i = 0; i < bins_normal.size(); i++) {
                    fitness += bins_normal[i]->parts[0]->smallitem ? 0 : 1;
                }
                for (int i = 0; i < bins_back.size(); i++) {
                    fitness += bins_back[i]->parts[0]->smallitem ? 0 : 1;
                }
                return -(fitness);
            }

            double getunity() {
                double last_unity = 1;
                for (int i = 0; i < bins_normal.size(); i++) {
                    if (bins_normal[i]->unity < last_unity) {
                        last_unity = bins_normal[i]->unity;
                    }
                }
                for (int i = 0; i < bins_back.size(); i++) {
                    if (bins_back[i]->unity < last_unity) {
                        last_unity = bins_back[i]->unity;
                    }
                }
                return (bins_normal.size() + bins_back.size() - 1 + last_unity);
            }

            double getareasum() {
                double areasum = 0;
                for (int i = 0; i < bins_normal.size(); i++) {
                    areasum += bins_normal[i]->fillarea;
                }
                for (int i = 0; i < bins_back.size(); i++) {
                    areasum += bins_back[i]->fillarea;
                }
                return areasum;
            }

            int getPartnum() {
                int areasum = 0;
                for (int i = 0; i < bins_normal.size(); i++) {
                    areasum += bins_normal[i]->parts.size();
                }
                for (int i = 0; i < bins_back.size(); i++) {
                    areasum += bins_back[i]->parts.size();
                }
                return areasum;
            }

            void update_area(int type) {
                if (type == 0) {
                    for (int i = 0; i < bins_back.size(); i++) {
                        bins_back[i]->update_area_sum();
                    }
                } else {
                    for (int i = 0; i < bins_normal.size(); i++) {
                        bins_normal[i]->update_area_sum();
                    }
                }
            }

            int calculate_last_average(bool back, bool last) const {
                int last_index = 0, last2_index = 0;
                const Bin_Ptr_V *temp_bins;
                if (back)
                    temp_bins = &bins_back;
                else temp_bins = &bins_normal;
                if (temp_bins->size() == 0) return -1;
                else if (temp_bins->size() == 1) return 0;
                if ((*temp_bins)[0]->averagearea < (*temp_bins)[1]->averagearea) {
                    last_index = 0;
                    last2_index = 1;
                } else {
                    last_index = 1;
                    last2_index = 0;
                };
                for (int i = 2; i < (*temp_bins).size(); i++) {
                    if ((*temp_bins)[i]->averagearea < (*temp_bins)[last_index]->averagearea) {
                        last2_index = last_index;
                        last_index = i;
                    } else if ((*temp_bins)[i]->averagearea > (*temp_bins)[last_index]->averagearea &&
                               (*temp_bins)[i]->averagearea > (*temp_bins)[last_index]->averagearea) {
                        last2_index = i;
                    }
                }
                assert(last_index != -1 && last2_index != -1);
                return last ? last_index : last2_index;
            }

            int calculate_last_fill(bool back, bool last) {
                int last_index = 0, last2_index = 0;
                Bin_Ptr_V *temp_bins;
                if (back) temp_bins = &bins_back;
                else temp_bins = &bins_normal;
                if (temp_bins->size() == 0) return -1;
                else if (temp_bins->size() == 1) return 0;

                if ((*temp_bins)[0]->fillarea < (*temp_bins)[1]->fillarea) {
                    last_index = 0;
                    last2_index = 1;
                } else {
                    last_index = 1;
                    last2_index = 0;
                };
                for (int i = 2; i < (*temp_bins).size(); i++) {
                    if ((*temp_bins)[i]->fillarea < (*temp_bins)[last_index]->fillarea) {
                        last2_index = last_index;
                        last_index = i;
                    } else if ((*temp_bins)[i]->fillarea < (*temp_bins)[last2_index]->fillarea &&
                            (*temp_bins)[i]->fillarea >= (*temp_bins)[last_index]->fillarea) {
                        last2_index = i;
                    }
                }
                assert(last_index != -1 && last2_index != -1);
                return last ? last_index : last2_index;
            }

            int insert_fill_rate1(bool back, std::mt19937 &mt) {
                Bin_Ptr_V *temp_bins;
                if (back) temp_bins = &bins_back;
                else temp_bins = &bins_normal;
                vector<int> sample_indexs;
                double unity_rate1 = back ? unity_rate1_back : unity_rate1_normal;
                for (int i = 0; i < (*temp_bins).size(); i++) {
                    if ((*temp_bins)[i]->unity < unity_rate1) {
                        sample_indexs.push_back(i);
                    }
                }
                if (sample_indexs.size() == 0) return -1;
                vector<int> output;
                std::sample(sample_indexs.begin(), sample_indexs.end(), std::back_inserter(output), 1, mt);
                return output[0];
            }

            vector<int> group_fill_rate2(bool back, std::mt19937 &mt) {
                Bin_Ptr_V *temp_bins;
                if (back) temp_bins = &bins_back;
                else temp_bins = &bins_normal;
                vector<int> sample_indexs;
                for (int i = 0; i < (*temp_bins).size(); i++) {
                    double unity_rate2 = back ? unity_rate2_back : unity_rate2_normal;
                    if ((*temp_bins)[i]->unity < unity_rate2) {
                        sample_indexs.push_back(i);
                    }
                }
                vector<int> output;
                if (sample_indexs.size() < 2) return output;
                std::sample(sample_indexs.begin(), sample_indexs.end(), std::back_inserter(output), 2, mt);
                return output;
            }

            vector<int> group_fill_last() {
                int last_index = 0, last2_index = 0;
                Bin_Ptr_V *temp_bins;
                temp_bins = &bins_normal;
                if (temp_bins->size() < 2) return {};
                if ((*temp_bins)[0]->fillarea < (*temp_bins)[1]->fillarea) {
                    last_index = 0;
                    last2_index = 1;
                } else {
                    last_index = 1;
                    last2_index = 0;
                };
                for (int i = 2; i < (*temp_bins).size(); i++) {
                    if ((*temp_bins)[i]->fillarea < (*temp_bins)[last_index]->fillarea) {
                        last2_index = last_index;
                        last_index = i;
                    } else if ((*temp_bins)[i]->fillarea > (*temp_bins)[last_index]->fillarea &&
                               (*temp_bins)[i]->fillarea > (*temp_bins)[last_index]->fillarea) {
                        last2_index = i;
                    }
                }
                assert(last_index != -1 && last2_index != -1);
                return {last_index, last2_index};
            }

            void update_rate(int flag) {
                float rate1, rate2, rate3;
                float frac;
                if (ratio >= 0.8) frac = 1;
                else if (ratio <= 0.15) frac = 0;
                else {
                    frac = (ratio - 0.15) / 0.66;
                }
                if (flag == 0) {
                    rate1 = 0.25 * frac;
                    rate2 = 1 * frac;
                    rate3 = 1 * frac;
                } else if (flag == 1) {
                    rate1 = 0.1 * frac;
                    rate2 = 0.33 * frac;
                    rate3 = 0.33 * frac;
                }
                vector<double> unity_vec;
                for (int i = 0; i < bins_back.size(); i++) {
                    unity_vec.push_back(bins_back[i]->unity);
                }
                for (int i = 0; i < bins_normal.size(); i++) {
                    unity_vec.push_back(bins_normal[i]->unity);
                }
                std::sort(unity_vec.begin(), unity_vec.end());
                unity_rate1_back = unity_vec[min(int((1 - rate1) * unity_vec.size()), int(unity_vec.size() - 1))];
                unity_rate2_back = unity_vec[min(int((1 - rate2) * unity_vec.size()), int(unity_vec.size() - 1))];
                unity_rate3_back = unity_vec[min(int((1 - rate3) * unity_vec.size()), int(unity_vec.size() - 1))];

//                unity_vec.clear();
//                for (int i = 0; i < bins_normal.size(); i++) {
//                    unity_vec.push_back(bins_normal[i]->unity);
//                }
//                std::sort(unity_vec.begin(), unity_vec.end());
                unity_rate1_normal = unity_vec[min(int((1 - rate1) * unity_vec.size()), int(unity_vec.size() - 1))];
                unity_rate2_normal = unity_vec[min(int((1 - rate2) * unity_vec.size()), int(unity_vec.size() - 1))];
                unity_rate3_normal = unity_vec[min(int((1 - rate3) * unity_vec.size()), int(unity_vec.size() - 1))];
            }

            double unity_rate1_back = 1;
            double unity_rate2_back = 1;
            double unity_rate3_back = 1;
            double unity_rate1_normal = 1;
            double unity_rate2_normal = 1;
            double unity_rate3_normal = 1;
            double ratio;
            bool junyun;
        private:;
        };

    };
};


#endif //PACKING_SOLUTION_PACKSOLUTION_H
