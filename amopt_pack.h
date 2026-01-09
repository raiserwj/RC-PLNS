#ifndef AMOPT_PACK_H
#define AMOPT_PACK_H

#include <iostream>
#include<vector>
#include <string>
#include<json/json.h>
#include "iso646.h"
#include "algorithm_error.h"
#include "CarvingMachine/CarvingMachine.h"
#include "IrregularPacking/IrregularPacking.h"
namespace amopt {
    namespace amopt_pack {
        Error Compute(string str, Json::Value &result_list,string filaname);
        Error ComputeIr(string str, Json::Value &result_list,string filaname);
    }
}

#endif