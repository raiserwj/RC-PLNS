//
// Created by raiser on 2021/12/13.
//

#include "CarvingMachine.h"
#include <chrono>
#include <sys/time.h>
#include<math.h>
#include <chrono>
using Clock = std::chrono::steady_clock;
using microseconds = std::chrono::microseconds;
amopt_pack::CarvingMachine::CarvingMachine() {
    std::cout << "carving machine initial" << std::endl;
}

amopt_pack::CarvingMachine::~CarvingMachine() {
    std::cout << "carving machine delete";
}

bool check(vector<float> a, vector<float> b, vector<float> c) {
    return ((b[0] == a[0] && b[1] == c[1]) || (b[1] == a[1] && b[0] == c[0]));
}

Error amopt_pack::CarvingMachine::Input(string str) {
    Json::CharReaderBuilder rbuilder;
    Json::CharReader *reader(rbuilder.newCharReader());
//    rbuilder["collectComments"] = false;
    Json::Value root_group;
    JSONCPP_STRING errs;
    if (!reader->parse(str.data(), str.data() + str.size(), &root_group, &errs)) {
        return Error::INPUT_FORMAT_ERROR;
    }
    delete reader;

    int num_small = 0;
    for (int i = 0; i < root_group["plates"].size(); i++) {
        for (int j = 0; j < root_group["plates"][i]["number"].asInt(); j++) {
            Bin_Ptr bin_(new amopt_pack::Bin(root_group["plates"][i]["width"].asDouble(),
                                             root_group["plates"][i]["height"].asDouble(), true));
            bins.push_back(bin_);
        }
    }
    for (int i = 0; i < root_group["items"].size(); i++) {
//        if(root_group["items"][i]["BackFrontPriority"].asBool()==false)
//        if(root_group["items"][i]["BackFrontPriority"].asBool()==true)
//        {
//            continue;
//        }
        float x_min = 100000, x_max = -100000, y_min = 100000, y_max = -100000;
        std::vector<std::vector<float>> points;
        for (int j = 0; j < root_group["items"][i]["points"].size(); j++) {
            std::vector<float> point;
            point.push_back(root_group["items"][i]["points"][j][0].asFloat());
            point.push_back(root_group["items"][i]["points"][j][1].asFloat());
            points.push_back(point);
            if (point[0] < x_min) x_min = point[0];
            if (point[0] > x_max) x_max = point[0];
            if (point[1] < y_min) y_min = point[1];
            if (point[1] > y_max) y_max = point[1];
        }
        int type = 1 - (points.size() == 4 &&
                        (check(points[0], points[1], points[2]) && check(points[1], points[2], points[3]) &&
                         check(points[2], points[3], points[0])));

        bool rotate_bool = false;
        int rotate_degree = 0;
        vector<int> rotationDegrees;
        bool feasibleflag = false;
        for (int j = 0; j < root_group["items"][i]["rotate"].size(); j++) {
            rotationDegrees.push_back(root_group["items"][i]["rotate"][j].asInt());
            for (int k = 0; k < bins.size(); k++) {
                if ((-x_min + x_max) <= bins[k]->width and (-y_min + y_max) <= bins[k]->height) {
                    feasibleflag = true;
                    break;
                }
            }

            if (root_group["items"][i]["rotate"][j].asInt() == 90 ||
                root_group["items"][i]["rotate"][j].asInt() == 270) {
                for (int k = 0; k < bins.size(); k++) {
                    if ((-x_min + x_max) <= bins[k]->height and (-y_min + y_max) <= bins[k]->width) {
                        feasibleflag = true;
                        break;
                    }
                }
                rotate_bool = true;
                rotate_degree = root_group["items"][i]["rotate"][j].asInt();
                break;
            }
        }
        if (feasibleflag == false) {
            continue;
        }
        rotate_bool=false;
        Part_Ptr part_(new amopt_pack::Part());
        part_->Input(y_max - y_min,
                     x_max - x_min, polygonarea(points), points, {x_min, x_max, y_min, y_max}, rotationDegrees, root_group["items"][i]["id"].asString(),
                     root_group["items"][i]["BackFrontPriority"].asBool(),
                     false, rotate_bool, rotate_degree, type);
//            part_->Input(y_max - y_min,
//                     x_max - x_min, polygonarea(points), points, {x_min, x_max, y_min, y_max}, rotationDegrees, root_group["items"][i]["id"].asString(),
//                     false,
//                         root_group["items"][i]["BackFrontPriority"].asBool(), rotate_bool, rotate_degree, type);
//        part_->Input(y_max - y_min,
//                     x_max - x_min, polygonarea(points), points, {x_min, x_max, y_min, y_max}, {0, 90, 180, 270},
//                     root_group["items"][i]["id"].asString(),
//                     root_group["items"][i]["BackFrontPriority"].asBool(),
//                     root_group["items"][i]["smallItem"].asBool(), rotate_bool, rotate_degree, type);
        parts.push_back(part_);
        points.clear();
        if (part_->smallitem) {
            num_small += 1;
        }
    }
    ratio = (num_small + 0.0) / (parts.size());
    std::cout << ratio << std::endl;
    solution.ratio = ratio;
//    if (squareparts.size() != 0) {
//        ratio = (num_small + 0.0) / (squareparts.size());
//        solution.ratio = ratio;
//    } else {
//        solution.ratio = 0;
//    }
    configure = amopt_pack::Configure(root_group["prameters"]["remnantsSet"]["area"].asDouble(),
                                      root_group["prameters"]["remnantsSet"]["externSpace"].asDouble(),
                                      root_group["prameters"]["remnantsSet"]["height"].asDouble(),
                                      root_group["prameters"]["remnantsSet"]["width"].asDouble(),
                                      root_group["prameters"]["surplusDir"].asBool(),
                                      root_group["prameters"]["timeLimit"].asInt());
//    solution.junyun = root_group["prameters"]["uniform"].asBool();
    solution.junyun = false;
    return Error::INPUT_FORMAT_ERROR;
}

void amopt_pack::CarvingMachine::combination() {
    int successnum = 0;
    int successnum_ = 0;
    std::cout << "combine" << std::endl;
    double sumwaste = 0;
    double sumrect = 0;
    for (int i_ = 0; i_ <= 1; i_++) {
        Part_Ptr_V *tempParts = (i_ == 0) ? &backparts : &notbackparts;
        vector<int> usedlist;
        for (int j_ = 1; j_ <= 2; j_++) {


            vector<int> squareIndexs;
            for (int j = 0; j < (*tempParts).size(); j++) {
                if (tempParts->operator[](j)->type == 0) squareIndexs.push_back(j);
            }

            vector<int> notsquareIndexs;
            for (int j = 0; j < (*tempParts).size(); j++) {
                if (tempParts->operator[](j)->type == j_) notsquareIndexs.push_back(j);

            }

            sort(notsquareIndexs.begin(), notsquareIndexs.end());
            int notsquarenum = notsquareIndexs.size();
            std::cout << "notsquarenum=" << notsquarenum << std::endl;
            int totalSize = notsquarenum + squareIndexs.size();
            if (j_ == 1) {
                for (int i = 0; i < notsquarenum; i++) {
                    double polygonarea = 0;
                    int iIndex = notsquareIndexs[i];
                    Part_Ptr partI = (*tempParts)[iIndex];
                    float minxA = 100000, maxxA = -100000, minyA = 100000, maxyA = -100000;
                    float minxATemp = 100000, maxxATemp = -100000, minyATemp = 100000, maxyATemp = -100000;
                    for (int i1 = 0; i1 < partI->points.size(); i1++) {
                        if (partI->points[i1][0] < minxATemp) minxATemp = partI->points[i1][0];
                        if (partI->points[i1][0] > maxxATemp) maxxATemp = partI->points[i1][0];
                        if (partI->points[i1][1] < minyATemp) minyATemp = partI->points[i1][1];
                        if (partI->points[i1][1] > maxyATemp) maxyATemp = partI->points[i1][1];
                    }

                    sumwaste += (maxxATemp - minxATemp) * (maxyATemp - minyATemp);
                    sumrect += (maxxATemp - minxATemp) * (maxyATemp - minyATemp);
                    for (int j = 0; j < partI->points.size() - 1; j++) {
                        polygonarea = polygonarea + (partI->points[j][0] * partI->points[j + 1][1] -
                                                     partI->points[j][1] * partI->points[j + 1][0]);
                    }
                    polygonarea = polygonarea + (partI->points[partI->points.size() - 1][0] * partI->points[0][1] -
                                                 partI->points[partI->points.size() - 1][1] * partI->points[0][0]);
                    sumwaste = sumwaste - abs(polygonarea) / 2;
                    if ((maxxATemp - minxATemp) * (maxyATemp - minyATemp) - abs(polygonarea) / 2 <
                        0.05 * (maxxATemp - minxATemp) * (maxyATemp - minyATemp)) {
                        successnum_ += 1;
                    }
                    std::cout << "sumwaste=" << sumwaste << std::endl;
                }
            }
            std::cout << "sumrect=" << sumrect << std::endl;
            for (int j = 0; j < squareIndexs.size(); j++) {
                notsquareIndexs.push_back(squareIndexs[j]);
            }
            Graph G(2 * totalSize);
            vector<double> cost;

            map<int, vector<vector<float>>> marchingRet;
            for (int i = 0; i < notsquarenum; i++) {
                int iIndex = notsquareIndexs[i];
                Part_Ptr partI = (*tempParts)[iIndex];
                if (j_ == 1) {
                    for (int j = i + 1; j < notsquarenum; j++) {
                        int jIndex = notsquareIndexs[j];
                        Part_Ptr partJ = (*tempParts)[jIndex];
                        vector<vector<float>> ret = group_two_(partI->pointsCombine, partJ->points,
                                                               partI->rotate_degrees, partJ->rotate_degrees,
                                                               bins[0]->width, bins[0]->height);
                        if (ret.size() > 0 && ret[0][0] > 0) {
                            G.AddEdge(i, j);
                            cost.push_back(-int(ret[0][0]));
                            G.AddEdge(i + totalSize, j + totalSize);
                            cost.push_back(-int(ret[0][0]));
                            marchingRet[i * 10000 + j] = ret;
                        }
                    }
                }

                for (int j = 0; j < squareIndexs.size(); j++) {
                    int squareIndex = squareIndexs[j];
                    if (1.2 * (partI->size - partI->area) >= (*tempParts)[squareIndex]->size) {
                        Part_Ptr partJ = (*tempParts)[squareIndex];
                        vector<vector<float>> ret = group_two_(partI->pointsCombine, partJ->points,
                                                               partI->rotate_degrees, partJ->rotate_degrees,
                                                               bins[0]->width, bins[0]->height);
                        if (ret.size() > 0 && ret[0][0] > 0) {
                            G.AddEdge(i, j + notsquarenum);
                            cost.push_back(-int(ret[0][0]));
                            G.AddEdge(i + totalSize, j + notsquarenum + totalSize);
                            cost.push_back(-int(ret[0][0]));
                            marchingRet[i * 10000 + j + notsquarenum] = ret;
                        }
                    }
                }
                G.AddEdge(i, i + totalSize);
                cost.push_back(0);
            }
            for (int j = 0; j < squareIndexs.size(); j++) {
                G.AddEdge(j + notsquarenum, j + notsquarenum + totalSize);
                cost.push_back(0);
            }
            Matching M(G);
            pair<list<int>, double> solution = M.SolveMinimumCostPerfectMatching(cost);
            std::cout << "reduce" << solution.second << std::endl;
            list<int> matching(solution.first);
            map<int, int> edge;
            for (list<int>::iterator it = matching.begin(); it != matching.end(); it++) {
                pair<int, int> e = G.GetEdge(*it);
                if (e.first < totalSize and e.second < totalSize) {
                    int min_index = min(e.first, e.second);
                    int max_index = max(e.first, e.second);
                    vector<vector<float>> ret = marchingRet[min_index * 10000 + max_index];

                    Part_Ptr part_(new amopt_pack::Part());
                    part_->width = ret[0][1];
                    part_->height = ret[0][2];
                    part_->size = ret[0][1] * ret[0][2];
                    part_->type = j_ + 1;
                    part_->area = 0;
                    part_->rotate = false;
                    part_->rotate_degrees = {0};

                    if ((*tempParts)[notsquareIndexs[min_index]]->combineParts.size() == 0) {
                        (*tempParts)[notsquareIndexs[min_index]]->transform.rotation = ret[0][5];
                        (*tempParts)[notsquareIndexs[min_index]]->transform.xShift = ret[1][0];
                        (*tempParts)[notsquareIndexs[min_index]]->transform.yShift = ret[1][1];
                        (*tempParts)[notsquareIndexs[min_index]]->type = -1;
                        part_->combineParts.push_back((*tempParts)[notsquareIndexs[min_index]]);
                        part_->pointsCombine = {
                                rotate_polygon((*tempParts)[notsquareIndexs[min_index]]->points, ret[0][5], 0, 0)};
                        part_->area += (*tempParts)[notsquareIndexs[min_index]]->area;
                    } else {
                        (*tempParts)[notsquareIndexs[min_index]]->type = -1;
                        for (int j = 0; j < (*tempParts)[notsquareIndexs[min_index]]->combineParts.size(); j++) {
//                            (*tempParts)[notsquareIndexs[min_index]]->combineParts[j]->transform.rotation = ret[0][5];
                            (*tempParts)[notsquareIndexs[min_index]]->combineParts[j]->transform.xShift = ret[j + 1][0];
                            (*tempParts)[notsquareIndexs[min_index]]->combineParts[j]->transform.yShift = ret[j + 1][1];
                            part_->combineParts.push_back((*tempParts)[notsquareIndexs[min_index]]->combineParts[j]);
                            part_->area += (*tempParts)[notsquareIndexs[min_index]]->combineParts[j]->area;
                        }
                        rotate_polygons_(part_->pointsCombine, ret[0][5]);
                    }
                    part_->area += (*tempParts)[notsquareIndexs[max_index]]->area;
                    part_->combineParts.push_back((*tempParts)[notsquareIndexs[max_index]]);
                    part_->pointsCombine.push_back(
                            rotate_polygon((*tempParts)[notsquareIndexs[max_index]]->points, ret[0][6], ret[0][3],
                                           ret[0][4]));
                    (*tempParts).push_back(part_);

                    (*tempParts)[notsquareIndexs[max_index]]->transform.rotation = ret[0][6];
                    (*tempParts)[notsquareIndexs[max_index]]->transform.xShift = ret.back()[0];
                    (*tempParts)[notsquareIndexs[max_index]]->transform.yShift = ret.back()[1];
                    (*tempParts)[notsquareIndexs[max_index]]->type = -1;

                    std::cout << min_index << "  " << max_index << "success" << std::endl;
                    std::cout << (*tempParts)[notsquareIndexs[min_index]]->id1 << "   "
                              << (*tempParts)[notsquareIndexs[max_index]]->id1 << std::endl;
                    std::cout << "1111" << std::endl;
                    successnum += 1;
//                    if(min_index <notsquarenum){
//                        vector<int>::iterator it;
//                        it = find(usedlist.begin(), usedlist.end(), min_index);
//                        if (it == usedlist.end())
//                        {
//                            successnum_+=1;
//                            usedlist.push_back(min_index);
//                        }
//
//                    }
//                    if(max_index <notsquarenum and j_==1){
//                        vector<int>::iterator it;
//                        it = find(usedlist.begin(), usedlist.end(),max_index);
//                        if (it == usedlist.end())
//                        {
//                            successnum_+=1;
//                            usedlist.push_back(min_index);
//                        }
//                    }
                }
            }
        }
    }
    std::cout << "successnum=" << successnum << std::endl;
    std::cout << "successnum_=" << successnum_ << std::endl;
}

Error amopt_pack::CarvingMachine::compute(Json::Value &result_list) {
    std::cout << "compute start" << std::endl;
    auto t_d1_start = Clock::now();
    firstpack();
    optimize();
    auto t_d1_end = Clock::now();
    std::cout << "Time:         " <<  std::chrono::duration_cast<microseconds>(t_d1_end - t_d1_start).count()   / 1e6 << " s\n";
    output(result_list);
    Json::Value temp;
//    temp["back_num"] = int(solution.bins_back.size());
//    temp["unity"] = solution.getunity();
    std::cout << "compute finish" << std::endl;
    return COMPUTE_NO_ERROR;
}

void amopt_pack::CarvingMachine::firstpack() {
    backparts = {};
    notbackparts = {};
    std::vector<std::vector<std::vector<double> > > result_;
    std::vector<amopt_pack::Bin> bins_ = {};
    std::vector<amopt_pack::Part *> parts_ = {};
    for (int i = 0; i < parts.size(); i++) {
        if (parts[i]->Showbackfrontpriority()) {
            backparts.push_back(parts[i]);
        } else {
            notbackparts.push_back(parts[i]);
        }
    }
//    gettimeofday(&start, NULL);
    combination();
//    gettimeofday(&end, NULL);
//    std::cout<<"time="<<(end.tv_sec - start.tv_sec)<<std::endl;
//    exit(0);
    std::cout << backparts.size() << "  " << notbackparts.size() << std::endl;

    double area_sum = 0;
    for (int i = 0; i < backparts.size(); i++) {
        if (backparts[i]->type == -1) continue;
        area_sum += backparts[i]->size;
    }
    std::cout << "back area:" << area_sum << std::endl;

    double notback_area_sum = 0;
    for (int i = 0; i < notbackparts.size(); i++) {
        if (notbackparts[i]->type == -1) continue;
        notback_area_sum += notbackparts[i]->size;
    }
    std::cout << "notback area:" << notback_area_sum << std::endl;
    std::cout << "total area:" << notback_area_sum + area_sum << std::endl;
    for (int i = 0; i < backparts.size(); i++) {
        for (int j = i; j < backparts.size(); j++) {
            if (backparts[i]->size < backparts[j]->size) {
                swap(backparts[i], backparts[j]);
            }
        }
    }
    for (int i = 0; i < backparts.size(); i++) {
        if (backparts[i]->type == -1) continue;
        bool check = false;
        for (int j = 0; j < solution.bins_back.size(); j++) {
            Part_Ptr_V parts = insert(solution.bins_back[j], backparts[i], mt, true, true);
            if (parts.size() == 0) {
                check = true;
                break;
            }
        }
        if (not check) {
            Bin_Ptr p_bin(new amopt_pack::Bin(bins[solution.bins_back.size()]->width,
                                              bins[solution.bins_back.size()]->height, true));
            Part_Ptr_V parts = insert(p_bin, backparts[i], mt, false, true);
            assert(parts.size() == 0);
            solution.bins_back.push_back(p_bin);
        }
    }
    if (backparts.size() != 0) {
        solution.update_area(0);
        auto repair = Repair();
        std::cout << "initial cost  " << solution.bins_back.size() << "  " << solution.getcost() << "  "
                  << solution.cost << "  " << solution.getunity()
                  << std::endl;
        ofstream fout1;
        Json::Value result_list1;
        output(result_list1);
        fout1.open("output_back_begin.json");
        Json::StyledWriter writer1;
        fout1 << writer1.write(result_list1) << std::endl;
        fout1.close();
        PackSolution best_solution = solution;
        auto destroy0 = DestroyReconstruct(100, 0.1, 7);
        gettimeofday(&start, NULL);
        while(true) {
            destroy0.destroy_solution(solution, mt);
            repair.repair_solution(solution, mt);
            if (solution.getcost() < best_solution.getcost()) {
                best_solution = solution;
            } else {
                solution = best_solution;
            }
            gettimeofday(&end, NULL);
            if (end.tv_sec - start.tv_sec > (60 * 0.2)) break;
        }
//        std::cout << solution.getareasum() << "  " << solution.getunity() << std::endl;
        solution = best_solution;
        std::cout << solution.getareasum() << "  " << solution.getunity() << std::endl;
        int bin_num = solution.bins_back.size() ;
        std::cout << "bin_back number:" << bin_num << std::endl;
        ofstream fout2;
        Json::Value result_list2;
        output(result_list2);
        fout2.open("output_back_over.json");
        Json::StyledWriter writer2;
        fout2 << writer2.write(result_list2) << std::endl;
        fout2.close();
    }
//    exit(0);
    for (int i = 0; i < notbackparts.size(); i++) {
        for (int j = i; j < notbackparts.size(); j++) {
            if (notbackparts[i]->size < notbackparts[j]->size) {
                swap(notbackparts[i], notbackparts[j]);
            }
        }
    }
    for (int i = 0; i < notbackparts.size(); i++) {
        if (notbackparts[i]->type == -1) continue;
        bool check = false;
        for (int j = 0; j < solution.bins_back.size(); j++) {
            Part_Ptr_V parts = insert(solution.bins_back[j], notbackparts[i], mt, true, true);
            if (parts.size() == 0) {
                check = true;
                break;
            }
        }

        if (not check) {
            check = false;
            for (int j = 0; j < solution.bins_normal.size(); j++) {
                Part_Ptr_V parts = insert(solution.bins_normal[j], notbackparts[i], mt, true, true);
                if (parts.size() == 0) {
                    check = true;
                    break;
                }
            }
            if (check) continue;

            Bin_Ptr p_bin(new amopt_pack::Bin(bins[solution.bins_back.size() + solution.bins_normal.size()]->width,
                                              bins[solution.bins_back.size() + solution.bins_normal.size()]->height,
                                              false));
            Part_Ptr_V parts = insert(p_bin, notbackparts[i], mt, true, true);
            if (parts.size() != 0) {
                std::cout << notbackparts[i]->width << "  " << notbackparts[i]->height << std::endl;
                std::cout << bins[solution.bins_back.size() + solution.bins_normal.size()]->width << "  "
                          << bins[solution.bins_back.size() + solution.bins_normal.size()]->height << std::endl;
                std::cout << "initial insert error" << std::endl;
//                exit(0);
            }
            solution.bins_normal.push_back(p_bin);
        }
    }
    solution.update_area(0);
    solution.update_area(1);
    std::cout << solution.getareasum() << "  " << solution.getunity() << std::endl;
    bestSolution = solution;

    for (int i = 0; i < solution.bins_normal.size(); i++) output_pack(solution.bins_normal[i]);

    int bin_num = solution.bins_back.size() + solution.bins_normal.size();
    std::cout << "bin number:" << bin_num << std::endl;
//    ofstream fout3;
//    Json::Value result_list3;
//    output(result_list3);
//
//    fout3.open("output_normal_begin.json");
//    Json::StyledWriter writer3;
//    fout3 << writer3.write(result_list3) << std::endl;
//    fout3.close();
}
void amopt_pack::CarvingMachine::move(mlpalns::DestroyMethod<PackSolution> &destroy,
                                      mlpalns::RepairMethod<PackSolution> &repair,
                                      std::mt19937 &mt, float tem)  {
    PackSolution temp_solution = solution;
    destroy.destroy_solution(temp_solution, mt);
    repair.repair_solution(temp_solution, mt);
    std::uniform_real_distribution<> dis(0, 1);
    Json::Value result_list1;
    Json::Value result_list2;
    float prob = dis(mt);
    if (solution.getunity() < bestSolution.getunity()) {
        bestSolution = solution;
    }
    double x1=solution.getcost();
    double x2=temp_solution.getcost();
    if (prob < exp((solution.getcost() - temp_solution.getcost()) / tem)) {
//        if( temp_solution.getunity() < bestSolution.getunity()){
//            std::cout<<1;
//        }
//       if( solution.bins_normal.size()>temp_solution.bins_normal.size()){
//            std::cout<<1;
//        }

        solution = temp_solution;
    }
//    ofstream fout1;
//    ofstream fout2;
//    if ((solution.getcost() - temp_solution.getcost())>1) {
//        std::cout<<solution.getcost()<<","<<temp_solution.getcost()<<std::endl;
//        output(result_list1);
//        fout1.open("output_1.json");
//        Json::StyledWriter writer1;
//        fout1 << writer1.write(result_list1) << std::endl;
//        fout1.close();
//        solution = temp_solution;
//        output(result_list2);
//        fout2.open("output_2.json");
//        Json::StyledWriter writer2;
//        fout2 << writer2.write(result_list2) << std::endl;
//        fout2.close();
//        exit(0);
//    }
}

void amopt_pack::CarvingMachine::optimize() {
    std::cout << "optimize start" << std::endl;
    mt.seed(std::time(0u));
//    mt.seed(3);
    std::uniform_real_distribution<float> dt(0, 1);
//    double frac_1 = 0.1, frac_2 = 0.05, frac_3 = 0.033;
    double frac_1 = 0.15, frac_2 = 0.1, frac_3 = 0.05;
    int num1 = 7, num2 = 5, num3 = 3;

    auto repair = Repair();
    //3类不同规模的重组
    auto regroup00 = DestroyReconstruct(100, frac_1, num1);
    auto regroup01 = DestroyReconstruct(100, frac_2, num2);
    auto regroup02 = DestroyReconstruct(100, frac_3, num3);
    auto regroup10 = DestroyReconstruct(101, frac_1, num1);
    auto regroup11 = DestroyReconstruct(101, frac_2, num2);
    auto regroup12 = DestroyReconstruct(101, frac_3, num3);
    auto regroup13 = DestroyReconstruct(101, frac_3, 2);
    vector<DestroyReconstruct> regroups = {regroup00, regroup01, regroup02, regroup10, regroup11, regroup12};

//    auto regrouplast0 = DestroyReconstruct(102, 0, num3);
//    auto regrouplast1 = DestroyReconstruct(103, 0, num3);
//    vector<DestroyReconstruct> regrouplasts = {regrouplast0, regrouplast1};

    //2类平均面积插入
    auto insertAverage0 = DestroyInsertAverage(200);
    auto insertAverage1 = DestroyInsertAverage(201);
    auto insertAverage2 = DestroyInsertAverage(202);
    auto insertAverage3 = DestroyInsertAverage(203);
    vector<DestroyInsertAverage> inserts = {insertAverage0, insertAverage1, insertAverage2, insertAverage3};

    //1类根据unity rate1的插入
    auto insertFill0 = DestroyInsertFill(302);
    auto insertFill1 = DestroyInsertFill(305);
    vector<DestroyInsertFill> insert_fills = {insertFill0, insertFill1};

    //1类根据unity rate2的重组
    auto regroup20 = DestroyReconstruct2(400);
    auto regroup21 = DestroyReconstruct2(401);
    vector<DestroyReconstruct2> regroup2 = {regroup20, regroup21};

    //1类填充面积最小两个重组
    auto regroup_last = DestroyReconstruct2(402);

    //2类普通板插入正反面板
    auto backinsert0 = DestroyBack(500);
    auto backinsert1 = DestroyBack(501);
    vector<DestroyBack> backs = {backinsert0, backinsert1};

    int step = 0;
    int index;
    double last_cost = solution.getcost();
    solution.update_rate(0);
//    while ((double) (clock() - start) / CLOCKS_PER_SEC < 60*configure.timeLimit) {
    float tem = 0.2;
    int timeflag=0;
    gettimeofday(&start, NULL);
    while (true) {
//        if (step>=200000) break;
        bool normal;
        gettimeofday(&end, NULL);
        double frac=-(end.tv_sec - start.tv_sec)/30*0.10+0.15;
        float back_ratio =
                (solution.bins_back.size() / 3) / (0.0 + solution.bins_back.size() + solution.bins_normal.size());
        if (dt(mt) < back_ratio) normal = false;
        else normal = true;
//        if (step < 10000) index = 0;
//        else if (step < 20000) index = 1;
        if (frac>0){
            auto regroup13 = DestroyReconstruct(101, frac, 2);
        }
//        else index = 2;
//        move(regroups[normal * 3 + index], repair, mt, tem); //rate3
        else{
            index = 2;
            move(regroups[normal * 3 + index], repair, mt, tem);
            move(regroup13, repair, mt, tem);
        }
        if (step % 5 == 0) {
            (dt(mt) < 0.5) ? index = 0 : index = 1;

            move(inserts[normal * 2 + index], repair, mt, tem);
        }

        if (step % 10 == 0) {
            move(insert_fills[normal], repair, mt, tem); // rate1
        }
        if (step % 300 == 0) {
            move(regroup_last, repair, mt, tem);
        }
        if (false || step % 500 == 149) {
            (dt(mt) < 0.5) ? index = 0 : index = 1;
            move(backs[index], repair, mt, tem);
        }
        if (step % 300 == 299) {
            move(regroup_last, repair, mt, tem);
        }
        if (step % 500 == 0) {
//            std::cout << step << "  " << tem << "  " << solution.getunity()
//                      << "  " << solution.bins_back.size() << std::endl;
            if (bestSolution.getcost() != last_cost) {
                last_cost = bestSolution.getcost();
                regroups[normal * 3 + 2].set(0.0, 2);
                solution.update_rate(0);
            } else {
                solution.update_rate(1);
                regroups[normal * 3 + 2].set(0.0, 3);
            }
        }
        if (end.tv_sec - start.tv_sec > (60 * timeflag)){
            std::cout << step << "  " << tem << "  " << solution.getunity();
            timeflag+=1;
        }
        step += 1;
        gettimeofday(&end, NULL);
        if (end.tv_sec - start.tv_sec > (60 )) break;//终点到0.05
        tem = 0.2 * exp(-(end.tv_sec - start.tv_sec + 0.0) / (60 ) * 2);
    }
    int bin_num = bestSolution.bins_back.size() + bestSolution.bins_normal.size();
    std::cout << "bin number:" << bin_num << std::endl;
    std::ofstream fout("record.txt", std::ios::app);
    if (!fout) return;
    fout << bin_num<<'\n';
}

void amopt_pack::CarvingMachine::output(Json::Value &result_list) {
//    result_list["back_num"] = int(bestSolution.bins_back.size());
//    result_list["unity"] = bestSolution.getunity();
//    return;
//    solution = bestSolution;

    bool flag;
//    int smallitemnum = 0;
//    for (int i = 0; i < solution.bins_normal.size(); i++) {
//        for(int j = 0; j < solution.bins_normal[i]->parts.size(); j++)
//        {
//            if (solution.bins_back[i]->parts[j]->smallitem == true) {
//                smallitemnum++;
//            }
//        }
//    }
//    std::cout << "output back=" << smallitemnum << std::endl;
//    std::cout << "output back=" << solution.bins_back.size() << std::endl;


    int smallitemnum = 0;
    for (int i = 0; i < solution.bins_back.size(); i++) {
        if (solution.bins_back[i]->parts[0]->smallitem == true) {
            smallitemnum++;
        }
    }
    for (int i = 0; i < solution.bins_normal.size(); i++) {
        if (solution.bins_normal[i]->parts[0]->smallitem == true) {
            smallitemnum++;
        }
    }
    std::cout << "output smallitem=" << smallitemnum << std::endl;
    for (int i = 0; i < solution.bins_back.size(); i++) {
        result_list["solutions"][i]["id"] = "";
        int j = 0;
        output_pack(solution.bins_back[i]);
        int j0 = 0;
        int index = 0;
        vector<int> list;
        vector<int>::iterator it;
        list.clear();
        for (j0 = 0; j0 < solution.bins_back[i]->parts.size(); j0++) {
            index = -1;
            for (int j2 = 0; j2 < solution.bins_back[i]->parts.size(); j2++) {
                it = std::find(list.begin(), list.end(), j2);
                if (it == list.end()) {
                    if (index == -1) index = j2;
                    else if (solution.bins_back[i]->parts[j2]->rect.x < solution.bins_back[i]->parts[index]->rect.x) {
                        index = j2;
                    }
                }
            }
            double distan = solution.bins_back[i]->parts[index]->rect.x;
            list.push_back(index);
            for (int j2 = 0; j2 < solution.bins_back[i]->parts.size(); j2++) {
                if ((solution.bins_back[i]->parts[j2]->rect.x + solution.bins_back[i]->parts[j2]->rect.width) <=
                    solution.bins_back[i]->parts[index]->rect.x and
                    solution.bins_back[i]->parts[j2]->rect.y <
                    (solution.bins_back[i]->parts[index]->rect.y + solution.bins_back[i]->parts[index]->rect.height) and
                    (solution.bins_back[i]->parts[j2]->rect.y + solution.bins_back[i]->parts[j2]->rect.height) >
                    solution.bins_back[i]->parts[index]->rect.y) {
                    distan = min(distan, solution.bins_back[i]->parts[index]->rect.x -
                                         solution.bins_back[i]->parts[j2]->rect.x -
                                         solution.bins_back[i]->parts[j2]->rect.width);
                }
            }
            solution.bins_back[i]->parts[index]->rect.x = solution.bins_back[i]->parts[index]->rect.x - distan;
        }
        list.clear();
        for (j0 = 0; j0 < solution.bins_back[i]->parts.size(); j0++) {
            index = -1;
            for (int j2 = 0; j2 < solution.bins_back[i]->parts.size(); j2++) {
                it = std::find(list.begin(), list.end(), j2);
                if (it == list.end()) {
                    if (index == -1) index = j2;
                    else if (solution.bins_back[i]->parts[j2]->rect.y < solution.bins_back[i]->parts[index]->rect.y) {
                        index = j2;
                    }
                }
            }
            double distan = solution.bins_back[i]->parts[index]->rect.y;
            list.push_back(index);
            for (int j2 = 0; j2 < solution.bins_back[i]->parts.size(); j2++) {
                if ((solution.bins_back[i]->parts[j2]->rect.y + solution.bins_back[i]->parts[j2]->rect.height) <=
                    solution.bins_back[i]->parts[index]->rect.y and
                    solution.bins_back[i]->parts[j2]->rect.x <
                    (solution.bins_back[i]->parts[index]->rect.x + solution.bins_back[i]->parts[index]->rect.width) and
                    (solution.bins_back[i]->parts[j2]->rect.x + solution.bins_back[i]->parts[j2]->rect.width) >
                    solution.bins_back[i]->parts[index]->rect.x) {
                    distan = min(distan, solution.bins_back[i]->parts[index]->rect.y -
                                         solution.bins_back[i]->parts[j2]->rect.y -
                                         solution.bins_back[i]->parts[j2]->rect.height);
                }
            }
            solution.bins_back[i]->parts[index]->rect.y = solution.bins_back[i]->parts[index]->rect.y - distan;
        }
        for (int partIndex = 0; partIndex < solution.bins_back[i]->parts.size(); partIndex += 1) {
            Part_Ptr part = solution.bins_back[i]->parts[partIndex];
            if (part->type == 0 || part->type == 1) {
                result_list["solutions"][i]["items"][j]["id"] = part->id1;
                result_list["solutions"][i]["items"][j]["centPt"][0] = part->rect.x + part->rect.width / 2;
                result_list["solutions"][i]["items"][j]["centPt"][1] = part->rect.y + part->rect.height / 2;
                result_list["solutions"][i]["items"][j]["rotate"] = (part->width == part->rect.width &&
                                                                     part->height == part->rect.height) ? 0
                                                                                                        : part->rotate_degree;
                j++;
            } else {
                std::cout << "size " << part->combineParts.size() << std::endl;
                for (auto tempPart: part->combineParts) {
                    result_list["solutions"][i]["items"][j]["id"] = tempPart->id1;
                    result_list["solutions"][i]["items"][j]["centPt"][0] =
                            part->rect.x + part->rect.width / 2 + tempPart->transform.xShift;
                    result_list["solutions"][i]["items"][j]["centPt"][1] =
                            part->rect.y + part->rect.height / 2 + tempPart->transform.yShift;
                    result_list["solutions"][i]["items"][j]["rotate"] = tempPart->transform.rotation;
                    j++;
                }
            }
            for (int k = partIndex + 1; k < solution.bins_back[i]->parts.size(); k++) {
                Rect a = part->rect;
                Rect b = solution.bins_back[i]->parts[k]->rect;
                if (!(a.x + a.width <= b.x ||
                      b.x + b.width <= a.x ||
                      a.y + a.height <= b.y ||
                      b.y + b.height <= a.y))
                    std::cout << "error overlap" << solution.getunity() << "  " << solution.ratio << std::endl;
            }
        }
    }

    std::cout << "output normal" << std::endl;
    int lastbin = solution.calculate_last_fill(0, 1);
    Bin_Ptr lastb;
    if (lastbin != -1) {
        last_insert(solution.bins_normal, lastbin, mt);
        if (solution.bins_normal[lastbin]->parts.size() == 0) {
            solution.bins_normal.erase(solution.bins_normal.begin() + lastbin);
            lastbin = -1;
        }
    }
    std::cout << "lastBin  " << lastbin << std::endl;
    for (int i = 0; i < solution.bins_normal.size(); i++) {
        int j = 0;
        if (i == lastbin) {
            if (!output_pack_last_new(solution.bins_normal[i], 1, false)) {
                std::cout << "lastBin defeat" << std::endl;
                output_pack(solution.bins_normal[i]);
            }
            std::cout << "lastBin success" << std::endl;
        } else output_pack(solution.bins_normal[i]);

        int j0 = 0;
        int index = 0;
        vector<int> list;
        vector<int>::iterator it;
        list.clear();
        for (j0 = 0; j0 < solution.bins_normal[i]->parts.size(); j0++) {
            index = -1;
            for (int j2 = 0; j2 < solution.bins_normal[i]->parts.size(); j2++) {
                it = std::find(list.begin(), list.end(), j2);
                if (it == list.end()) {
                    if (index == -1) index = j2;
                    else if (solution.bins_normal[i]->parts[j2]->rect.x <
                             solution.bins_normal[i]->parts[index]->rect.x) {
                        index = j2;
                    }
                }

            }
            double distan = solution.bins_normal[i]->parts[index]->rect.x;
            list.push_back(index);
            for (int j2 = 0; j2 < solution.bins_normal[i]->parts.size(); j2++) {
                if ((solution.bins_normal[i]->parts[j2]->rect.x + solution.bins_normal[i]->parts[j2]->rect.width) <=
                    solution.bins_normal[i]->parts[index]->rect.x and
                    solution.bins_normal[i]->parts[j2]->rect.y < (solution.bins_normal[i]->parts[index]->rect.y +
                                                                  solution.bins_normal[i]->parts[index]->rect.height) and
                    (solution.bins_normal[i]->parts[j2]->rect.y + solution.bins_normal[i]->parts[j2]->rect.height) >
                    solution.bins_normal[i]->parts[index]->rect.y) {
                    distan = min(distan, solution.bins_normal[i]->parts[index]->rect.x -
                                         solution.bins_normal[i]->parts[j2]->rect.x -
                                         solution.bins_normal[i]->parts[j2]->rect.width);
                }
            }
            solution.bins_normal[i]->parts[index]->rect.x = solution.bins_normal[i]->parts[index]->rect.x - distan;
        }
        list.clear();
        for (j0 = 0; j0 < solution.bins_normal[i]->parts.size(); j0++) {
            index = -1;
            for (int j2 = 0; j2 < solution.bins_normal[i]->parts.size(); j2++) {
                it = std::find(list.begin(), list.end(), j2);
                if (it == list.end()) {
                    if (index == -1) index = j2;
                    else if (solution.bins_normal[i]->parts[j2]->rect.y <
                             solution.bins_normal[i]->parts[index]->rect.y) {
                        index = j2;
                    }
                }

            }
            double distan = solution.bins_normal[i]->parts[index]->rect.y;
            list.push_back(index);
            for (int j2 = 0; j2 < solution.bins_normal[i]->parts.size(); j2++) {
                if ((solution.bins_normal[i]->parts[j2]->rect.y + solution.bins_normal[i]->parts[j2]->rect.height) <=
                    solution.bins_normal[i]->parts[index]->rect.y and
                    solution.bins_normal[i]->parts[j2]->rect.x < (solution.bins_normal[i]->parts[index]->rect.x +
                                                                  solution.bins_normal[i]->parts[index]->rect.width) and
                    (solution.bins_normal[i]->parts[j2]->rect.x + solution.bins_normal[i]->parts[j2]->rect.width) >
                    solution.bins_normal[i]->parts[index]->rect.x) {
                    distan = min(distan, solution.bins_normal[i]->parts[index]->rect.y -
                                         solution.bins_normal[i]->parts[j2]->rect.y -
                                         solution.bins_normal[i]->parts[j2]->rect.height);
                }
            }
            solution.bins_normal[i]->parts[index]->rect.y = solution.bins_normal[i]->parts[index]->rect.y - distan;

        }
        result_list["solutions"][int(i + solution.bins_back.size())]["id"] = "";

        j = 0;
        for (int partIndex = 0; partIndex < solution.bins_normal[i]->parts.size(); partIndex += 1) {
            Part_Ptr part = solution.bins_normal[i]->parts[partIndex];
            if (part->type == 0 || part->type == 1) {
                result_list["solutions"][int(i + solution.bins_back.size())]["items"][j]["id"] = part->id1;
                result_list["solutions"][int(i + solution.bins_back.size())]["items"][j]["centPt"][0] =
                        part->rect.x + part->rect.width / 2;
                result_list["solutions"][int(i + solution.bins_back.size())]["items"][j]["centPt"][1] =
                        part->rect.y + part->rect.height / 2;
                result_list["solutions"][int(i + solution.bins_back.size())]["items"][j]["rotate"] = (part->width ==
                                                                                                      part->rect.width &&
                                                                                                      part->height ==
                                                                                                      part->rect.height)
                                                                                                     ? 0
                                                                                                     : part->rotate_degree;
                j++;
            } else {
                std::cout << "size " << part->combineParts.size() << std::endl;
                for (auto tempPart: part->combineParts) {
                    result_list["solutions"][int(i + solution.bins_back.size())]["items"][j]["id"] = tempPart->id1;
                    result_list["solutions"][int(i + solution.bins_back.size())]["items"][j]["centPt"][0] =
                            part->rect.x + part->rect.width / 2 + tempPart->transform.xShift;
                    result_list["solutions"][int(i + solution.bins_back.size())]["items"][j]["centPt"][1] =
                            part->rect.y + part->rect.height / 2 + tempPart->transform.yShift;
                    result_list["solutions"][int(
                            i + solution.bins_back.size())]["items"][j]["rotate"] = tempPart->transform.rotation;
                    j++;
                }
            }
            for (int k = partIndex + 1; k < solution.bins_normal[i]->parts.size(); k++) {
                Rect a = part->rect;
                Rect b = solution.bins_normal[i]->parts[k]->rect;
                if (!(a.x + a.width <= b.x + 0.1 ||
                      b.x + b.width <= a.x + 0.1 ||
                      a.y + a.height <= b.y + 0.1 ||
                      b.y + b.height <= a.y + 0.1) ||
                    a.x <= -0.1 ||
                    a.y <= -0.1 ||
                    a.height <= 0.1) {
                    std::cout << "error overlap" << solution.getunity() << "  " << solution.ratio << std::endl;
                    std::cout << i << std::endl;
                    std::cout << a.x << " " << a.width << " " << a.y << " " << a.height << "\n";
                    std::cout << b.x << " " << b.width << " " << b.y << " " << b.height;
                }
            }
        }
    }
    result_list["back_num"] = int(solution.bins_back.size());
    result_list["unity"] = solution.getunity();
    std::cout << "bin number:" << solution.getunity() << std::endl;
}