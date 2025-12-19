//
// Created by raiser on 2021/12/29.
//

#ifndef PACKING_SOLUTION_VISITOR_H
#define PACKING_SOLUTION_VISITOR_H

#include "PALNS/AlgorithmVisitor.h"
#include "PackSolution.h"

struct PackVisitor : public mlpalns::AlgorithmVisitor<PackSolution> {
    PackVisitor(){};

    void on_algorithm_start(std::vector<mlpalns::DestroyMethod<PackSolution>*>& destroy,
                            std::vector<mlpalns::RepairMethod<PackSolution>*>& repair,
                            const std::vector<std::string>& dnames,
                            const std::vector<std::string>& rnames) {};

    void on_prerun_end(std::vector<mlpalns::DestroyMethod<PackSolution>*>& destroy,
                       std::vector<mlpalns::RepairMethod<PackSolution>*>& repair) {};

    void on_iteration_end(mlpalns::AlgorithmStatus<PackSolution>& alg_status) {};

    void on_many_iters_without_improvement(
            std::vector<mlpalns::DestroyMethod<PackSolution>*>& destroy,
            std::vector<mlpalns::RepairMethod<PackSolution>*>& repair
    ) {};

    ~PackVisitor() {};
};

#endif //PACKING_SOLUTION_VISITOR_H
