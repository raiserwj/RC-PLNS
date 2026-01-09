//
// Created by dxy on 2021/12/13.
//

#ifndef PACKING_SOLUTION_REPAIR_H
#define PACKING_SOLUTION_REPAIR_H

#include "CarvingMachine.h"
#include "amopt_pack_base.h"
#include "PackSolution.h"
#include "PALNS/RepairMethod.h"
#include "numeric"
#include "PALNS/PALNS.h"
#include "algorithm"
#include "random"
#include "CBFF.h"

using std::set;
using namespace mlpalns;
using namespace amopt_pack;
namespace amopt {
    namespace amopt_pack {
        struct Repair : public mlpalns::RepairMethod<PackSolution> {
            Repair() = default;

            ~Repair() override = default;

            std::unique_ptr<mlpalns::RepairMethod<PackSolution>> clone() const override {
                return std::make_unique<Repair>(*this);
            }

            void repair_solution(PackSolution &solution, std::mt19937 &mt) override {
                if (solution.destroy_set.size() == 0) return;
                if (solution.type == 100 || solution.type == 101 || solution.type == 102 || solution.type == 103 ||
                    solution.type == 400 || solution.type == 401 || solution.type == 402) {
                    bool flag;
                    Bin_Ptr_V *temp_bins;
                    if (solution.type == 100 || solution.type == 102 || solution.type == 400)
                        temp_bins = &solution.bins_back;
                    else temp_bins = &solution.bins_normal;
                    Bin_Ptr_V bins = regroup(solution.destroy_set, &flag, mt, true,
                                             solution.junyun, solution.ratio < 0.7);
                    if (flag) {
                        for (int i = 0; i < bins.size(); i++) {
                            bins[i]->update_area_sum();
                            assert(bins[i]->parts.size() != 0);
                        }
                        (*temp_bins).insert((*temp_bins).end(), bins.begin(), bins.end());
                    } else {
                        (*temp_bins).insert((*temp_bins).end(), solution.destroy_set.begin(),
                                            solution.destroy_set.end());
                    }
                } else if ((solution.type >= 200 && solution.type <= 203) ||
                           (solution.type >= 300 && solution.type <= 305)) {
                    Bin_Ptr_V *temp_bins;
                    Bin_Ptr destroy_bin = solution.destroy_set[0];
                    if (solution.type <= 201 || (solution.type >= 300 && solution.type <= 302))
                        temp_bins = &solution.bins_back;
                    else temp_bins = &solution.bins_normal;
                    if (temp_bins->size() == 0) {
                        (*temp_bins).push_back(destroy_bin);
                        solution.destroy_set.clear();
                        return;
                    }
                    Part_Ptr_V parts(destroy_bin->parts);
                    std::shuffle(parts.begin(), parts.end(), mt);

                    Part_Ptr_V remain_parts;
                    map<int, Bin_Ptr> destroy_bins;

                    for (int i = 0; i < parts.size(); i++) {
                        std::uniform_int_distribution<int> dt(0, (*temp_bins).size() - 1);
                        int index = dt(mt);
                        if (destroy_bins.find(index) == destroy_bins.end()) {
                            Bin_Ptr temp_bin(
                                    new amopt_pack::Bin((*temp_bins)[index]->width, (*temp_bins)[index]->height,
                                                        (*temp_bins)[index]->backornot, (*temp_bins)[index]->parts));
                            destroy_bins[index] = temp_bin;
                        }
                        Part_Ptr_V parts_ret = insert(destroy_bins[index], parts[i], mt, false, true, solution.ratio < 0.7);
                        remain_parts.insert(remain_parts.end(), parts_ret.begin(), parts_ret.end());
                    }

                    if (remain_parts.size() != 0) {
                        Bin_Ptr temp_bin(
                                new amopt_pack::Bin(destroy_bin->width, destroy_bin->height, destroy_bin->backornot,
                                                    remain_parts));
                        bool flag;
                        Bin_Ptr_V return_bin = regroup({temp_bin}, &flag, mt, true, solution.junyun, solution.ratio < 0.7);
                        if (flag) {
                            std::map<int, Bin_Ptr>::iterator it;
                            for (it = destroy_bins.begin(); it != destroy_bins.end(); ++it) {
                                (*temp_bins)[it->first] = it->second;
                                assert(it->second->parts.size() != 0);
                                it->second->update_area_sum();
                            }
                            return_bin[0]->update_area_sum();
                            (*temp_bins).push_back(return_bin[0]);
                        } else {
                            (*temp_bins).push_back(destroy_bin);
                            assert(destroy_bin->parts.size() != 0);
                        }
                    } else {
                        std::map<int, Bin_Ptr>::iterator it;
                        for (it = destroy_bins.begin(); it != destroy_bins.end(); ++it) {
                            (*temp_bins)[it->first] = it->second;
                            assert(it->second->parts.size() != 0);
                            it->second->update_area_sum();
                        }
                    }
                } else if (solution.type == 500 || solution.type == 501) {
                    if (solution.bins_back.size() == 0) {
                        solution.bins_normal.push_back(solution.destroy_set[0]);
                    } else {
                        Bin_Ptr destroy_bin = solution.destroy_set[0];
                        Bin_Ptr_V *temp_bins;
                        std::shuffle(solution.bins_back.begin(), solution.bins_back.end(), mt);
                        temp_bins = &solution.bins_back;

                        Part_Ptr_V parts(destroy_bin->parts);
                        std::shuffle(parts.begin(), parts.end(), mt);

                        Part_Ptr_V remain_parts;
                        map<int, Bin_Ptr> destroy_bins;

                        for (int i = 0; i < parts.size(); i++) {
                            int index ;
                            const int n = static_cast<int>((*temp_bins).size());
                            std::vector<int> idx(n);
                            std::iota(idx.begin(), idx.end(), 0);
                            Part_Ptr_V parts_ret;
                            if (n > 10) {
                                std::shuffle(idx.begin(), idx.end(), mt);
                                idx.resize(10);
                            }
                            for(int k=0;k<idx.size();k++){
                                bool success=false;
                                index=idx[k];
                                if (destroy_bins.find(index) == destroy_bins.end()) {
                                    Part_Ptr_V base_parts;
                                    base_parts.reserve((*temp_bins)[index]->parts.size());
                                    for (auto &p : (*temp_bins)[index]->parts) {
                                        if (p != parts[i]) base_parts.push_back(p);  // 关键：过滤同地址 Part
                                    }

                                    Bin_Ptr temp_bin(new amopt_pack::Bin((*temp_bins)[index]->width,
                                                                         (*temp_bins)[index]->height,
                                                                         (*temp_bins)[index]->backornot,
                                                                         base_parts));
                                    destroy_bins[index] = temp_bin;
                                }
                                parts_ret = insert(destroy_bins[index], parts[i], mt, true, true, solution.ratio < 70);
                                double bestarea=0;
                                for (int j=0;j<parts_ret.size();j++){
                                    bestarea+=parts_ret[j]->area;
                                }
                                if(bestarea<1){
                                    success=true;
                                    break;
                                }
                                float score1MaxRects=9999;
                                float score2MaxRects=9999;
                                int successnum=0;
                                int setarea=0;
                                std::vector<int> successset={};
                                for (int j=0;j<destroy_bins[index]->parts.size();j++){
                                    if(destroy_bins[index]->parts[j]->area<bestarea && destroy_bins[index]->parts[j]->backfrontpriority==false){
                                        setarea+=destroy_bins[index]->parts[j]->area;
                                        successset.push_back(j);
                                        successnum+=1;
                                        Part_Ptr_V base_parts;
                                        base_parts.reserve((destroy_bins)[index]->parts.size());
                                        for (auto &p : (destroy_bins)[index]->parts) {
                                            if (p != parts[i]) base_parts.push_back(p);  // 关键：过滤同地址 Part
                                        }
                                        Bin_Ptr temp_bin(new amopt_pack::Bin(destroy_bins[index]->width,
                                                                             destroy_bins[index]->height,
                                                                             destroy_bins[index]->backornot,
                                                                             base_parts));
                                        temp_bin->parts.erase(temp_bin->parts.begin() + j);
                                        temp_bin->parts.push_back(parts[i]);
                                        vector<Part_Ptr_V> retMaxRects = maxRects(temp_bin->parts, temp_bin->width, temp_bin->height, 1, mt, &score1MaxRects, &score2MaxRects, false);
                                        if (score1MaxRects==0){
                                            success=true;
                                            parts_ret={destroy_bins[index]->parts[j]};
                                            destroy_bins[index]->parts=retMaxRects[0];
                                            break;
                                        }
                                    }
                                    if(successnum==2){
                                        successnum+=1;
                                        if(setarea<bestarea){
                                            Part_Ptr_V base_parts;
                                            base_parts.reserve((destroy_bins)[index]->parts.size());
                                            for (auto &p : (destroy_bins)[index]->parts) {
                                                if (p != parts[i]) base_parts.push_back(p);  // 关键：过滤同地址 Part
                                            }
                                            Bin_Ptr temp_bin(new amopt_pack::Bin(destroy_bins[index]->width,
                                                                                 destroy_bins[index]->height,
                                                                                 destroy_bins[index]->backornot,
                                                                                 base_parts));
                                            temp_bin->parts.erase(temp_bin->parts.begin() + successset[1]);
                                            temp_bin->parts.erase(temp_bin->parts.begin() + successset[0]);
                                            temp_bin->parts.push_back(parts[i]);
                                            vector<Part_Ptr_V> retMaxRects = maxRects(temp_bin->parts, temp_bin->width, temp_bin->height, 1, mt, &score1MaxRects, &score2MaxRects, false);
                                            if (score1MaxRects==0){
                                                success=true;
                                                parts_ret={destroy_bins[index]->parts[successset[1]],destroy_bins[index]->parts[successset[0]]};
                                                destroy_bins[index]->parts=retMaxRects[0];
                                                break;
                                            }
                                        }
                                    }
                                    if (successnum>=3){
                                        break;
                                    }
                                }
                                if(success){
                                    break;
                                }
                            }
                            remain_parts.insert(remain_parts.end(), parts_ret.begin(), parts_ret.end());

                        }

                        if (remain_parts.size() != 0) {
                            Bin_Ptr temp_bin(
                                    new amopt_pack::Bin(destroy_bin->width, destroy_bin->height, destroy_bin->backornot,
                                                        remain_parts));
                            bool flag;
                            Bin_Ptr_V return_bin = regroup({temp_bin}, &flag, mt, true, solution.junyun, solution.ratio < 70);
                            if (flag) {
                                std::map<int, Bin_Ptr>::iterator it;
                                for (it = destroy_bins.begin(); it != destroy_bins.end(); ++it) {
                                    (*temp_bins)[it->first] = it->second;
                                    assert(it->second->parts.size() != 0);
                                    it->second->update_area_sum();
                                }
                                return_bin[0]->update_area_sum();
                                solution.bins_normal.push_back(return_bin[0]);
                            } else {
                                solution.bins_normal.push_back(destroy_bin);
                                assert(destroy_bin->parts.size() != 0);
                            }
                        } else {
                            std::map<int, Bin_Ptr>::iterator it;
                            for (it = destroy_bins.begin(); it != destroy_bins.end(); ++it) {
                                (*temp_bins)[it->first] = it->second;
                                assert(it->second->parts.size() != 0);
                                it->second->update_area_sum();
                            }
                        }
                    }
                }
                solution.destroy_set.clear();
            }
        };
    };
};


#endif //PACKING_SOLUTION_REPAIR_H
