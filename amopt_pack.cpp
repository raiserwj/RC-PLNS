#include <sys/wait.h>
#include"amopt_pack.h"
#include <future>
namespace amopt {
    Error amopt_pack::Compute(string str, Json::Value &result_list, string filename) {
        std::cout << "pack compute" << std::endl;
        Json::CharReaderBuilder rbuilder;
        Json::CharReader *reader = rbuilder.newCharReader();
//        rbuilder["collectComments"] = false;
        Json::Value root_group;
        JSONCPP_STRING errs;
        if (!reader->parse(str.data(), str.data() + str.size(), &root_group, &errs)) {
            std::cout << "read error" << std::endl;
            return Error::INPUT_FORMAT_ERROR;
        }
        delete reader;

        std::cout << "CarvingMachine" << std::endl;
        CarvingMachine carvingmachine = CarvingMachine();
        carvingmachine.Input(str);
        carvingmachine.numitems();
        carvingmachine.compute(result_list);

        std::cout << "CarvingMachine Finish" << std::endl;
        return amopt::COMPUTE_NO_ERROR;
    }
    Error amopt_pack::ComputeIr(string str, Json::Value &result_list, string filename) {
        std::cout << "IrregularPacking compute" << std::endl;
        Json::CharReaderBuilder rbuilder;
        Json::CharReader *reader = rbuilder.newCharReader();
//        rbuilder["collectComments"] = false;
        Json::Value root_group;
        JSONCPP_STRING errs;
        if (!reader->parse(str.data(), str.data() + str.size(), &root_group, &errs)) {
            std::cout << "read error" << std::endl;
            return Error::INPUT_FORMAT_ERROR;
        }
        delete reader;

        std::cout << "IrregularPacking" << std::endl;
        IrregularPacking irregularpacking = IrregularPacking();
        irregularpacking.Input(str);
//        irregularpacking.numitems();
        irregularpacking.compute(result_list);

        std::cout << "IrregularPacking Finish" << std::endl;
        return amopt::COMPUTE_NO_ERROR;
    }
}


void Test(string filename, int i) {
    string s;
    std::stringstream str;
    fstream f;
    f.open(filename, ios::in);
    if (!f.is_open()) {
        cout << "Open json file error!" << endl;
        exit(0);
    }
    str << f.rdbuf();
    f.close();
    s = str.str();
    Json::Value result_list;
    amopt::Error err = amopt::amopt_pack::Compute(s, result_list, filename);
    if (err != amopt::Error::COMPUTE_NO_ERROR) {
        std::cout << "error" << err << std::endl;
    }
    std::cout << "finish calculate" << std::endl;
    ofstream fout;
    std::cout << "output_" + filename << std::endl;
    fout.open("output_" + filename);
    Json::StyledWriter writer1;
    fout << writer1.write(result_list) << std::endl;
    fout.close();
    ofstream of;
    of.open(filename+".txt",ios::out );
    of<<result_list["solutions"].size()<<std::endl;
    of.close();
}
void nfpc(string filename){
    Json::Value result_list;
    vector<vector<vector<int>>> rotation;
    string s;
    nfpClass nfpC;
    std::stringstream str;
    fstream f;
    f.open(filename, ios::in);
    if (!f.is_open()) {
        cout << "Open json file error!" << endl;
        exit(0);
    }
    str << f.rdbuf();
    f.close();
    s = str.str();
    Json::CharReaderBuilder rbuilder;
    Json::CharReader *reader = rbuilder.newCharReader();
//        rbuilder["collectComments"] = false;
    Json::Value root_group;
    JSONCPP_STRING errs;
    if (!reader->parse(s.data(), s.data() + s.size(), &root_group, &errs)) {
        std::cout << "read error" << std::endl;
    }
    delete reader;
    vector<vector<vector<vector<float>>>> datasets;
    for(int i1=0;i1<root_group["datasets"].size();i1++){
        datasets.push_back({});
        rotation.push_back({});
        for(int i2=0;i2<root_group["datasets"][i1].size();i2++){
            datasets[i1].push_back({});
            rotation[i1].push_back({0,180});
            for(int i3=0;i3<root_group["datasets"][i1][i2].size();i3++){
                datasets[i1][i2].push_back({{root_group["datasets"][i1][i2][i3][0].asFloat()}, {root_group["datasets"][i1][i2][i3][1].asFloat()}});
            }
        }
    }
    vector<vector<float>> nfpbetweenp;
    vector<vector<vector<vector<vector<vector<vector<float>>>>>>> nfp={};
    for(int dataset=0;dataset<datasets.size();dataset++){
        nfp.push_back({});
//        if(dataset!=10){
//            continue;
//        }
        for(int i=0;i<datasets[dataset].size();i++){
            nfp[dataset].push_back({});
            for(int j=0;j<datasets[dataset].size();j++){
                nfp[dataset][i].push_back({});
                for(int ri=0;ri<rotation[dataset][i].size();ri++){
                    nfp[dataset][i][j].push_back({});
                    for(int rj=0;rj<rotation[dataset][j].size();rj++){
                        nfpbetweenp=nfpC.nfpbetweentwo({rotate_polygon(datasets[dataset][i], rotation[dataset][i][ri], 0, 0)},rotate_polygon(datasets[dataset][j], rotation[dataset][j][rj], 0, 0),99999,99999);
                        nfp[dataset][i][j][ri].push_back(nfpbetweenp);
                    }
                }
            }
        }
        std::cout<<dataset<<std::endl;
    }
    result_list["nfp"]= {};
    for(int i1=0;i1<nfp.size();i1++){
        for(int i2=0;i2<nfp[i1].size();i2++){
            for(int i3=0;i3<nfp[i1][i2].size();i3++){
                for(int i4=0;i4<nfp[i1][i2][i3].size();i4++){
                    for(int i5=0;i5<nfp[i1][i2][i3][i4].size();i5++){
                        for(int i6=0;i6<nfp[i1][i2][i3][i4][i5].size();i6++){
                            result_list["nfp"][i1][i2][i3][i4][i5][i6][0]=nfp[i1][i2][i3][i4][i5][i6][0];
                            result_list["nfp"][i1][i2][i3][i4][i5][i6][1]=nfp[i1][i2][i3][i4][i5][i6][1];
                        }
                    }
                }
            }
        }
    }
    ofstream fout;
    std::cout << "2.txt" << std::endl;
    fout.open("2.txt" );
    Json::StyledWriter writer1;
    fout << writer1.write(result_list) << std::endl;
    fout.close();
}

void run_batch(const vector<string>& files) {
    // 启动一批子进程
    const int TIME_LIMIT_SEC = 150;  // 整批这 10 个实例允许的最大时间

    std::vector<pid_t> child_pids;
    child_pids.reserve(files.size());

    // 1) 启动所有子进程
    for (const auto& file : files) {
        pid_t pid = fork();
        if (pid == 0) {
            // 子进程：只做自己的事情
            Test(file, 0);
            _exit(0);   // 不要用 exit()
        } else if (pid < 0) {
            // fork 失败，清理已启动的子进程
            perror("fork");
            for (pid_t cpid : child_pids) {
                if (cpid > 0) {
                    kill(cpid, SIGKILL);
                }
            }
            return ;
        } else {
            // 父进程：记录子进程 pid
            child_pids.push_back(pid);
        }
    }

    // 2) 父进程轮询等待，带总时间上限
    auto start = std::chrono::steady_clock::now();
    size_t remaining = child_pids.size();

    while (remaining > 0) {
        // 尝试回收已经结束的子进程（非阻塞）
        for (pid_t &cpid : child_pids) {
            if (cpid <= 0) continue;  // 这个 pid 已经回收过了

            int status = 0;
            pid_t ret = waitpid(cpid, &status, WNOHANG);
            if (ret == cpid) {
                // 这个子进程结束了
                cpid = -1;
                --remaining;
            } else if (ret < 0) {
                if (errno == ECHILD) {
                    // 子进程已经不存在了（可能其它地方回收过）
                    cpid = -1;
                    --remaining;
                } else {
                    perror("waitpid");
                    cpid = -1;
                    --remaining;
                }
            }
        }

        if (remaining == 0) {
            break;  // 所有子进程都结束，正常退出
        }

        // 检查时间是否超过上限
        auto now = std::chrono::steady_clock::now();
        double elapsed =
                std::chrono::duration_cast<std::chrono::duration<double>>(now - start).count();

        if (elapsed > TIME_LIMIT_SEC) {
            std::cerr << "Timeout: killing unfinished children" << std::endl;

            // 杀掉仍然存活的子进程
            for (pid_t cpid : child_pids) {
                if (cpid > 0) {
                    kill(cpid, SIGKILL);
                }
            }

            // 尝试把它们都 wait 掉，避免僵尸进程
            int status = 0;
            while (waitpid(-1, &status, WNOHANG) > 0) {
                // no-op
            }

            break;  // 结束这批任务
        }

        // 避免忙等，稍微 sleep 一下
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
//    Test(files[0], 0);
}
static void Run10Parallel(const std::string& path, int baseSeed = 21) {
    std::vector<std::future<void>> futs;
    futs.reserve(10);

    for (int i = 0; i < 10; ++i) {
        const int seed = baseSeed + i;
        futs.emplace_back(std::async(std::launch::async, [path, seed]() {
            Test(path, seed);
        }));
    }

    // 等待这 10 个全部完成
    for (auto& f : futs) f.get();
}
int main() {
    Json::Value result_list;

    vector<string> batch1 = {
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_04.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_14.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_24.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_34.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_44.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_54.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_64.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_74.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_84.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_94.json"
    };

    vector<string> batch2 = {
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_03.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_13.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_23.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_33.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_43.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_53.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_63.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_73.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_83.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_93.json",
    };

    vector<string> batch3 = {
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_02.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_12.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_22.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_32.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_42.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_52.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_62.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_72.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_82.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_92.json"
    };
    vector<string> batch4 = {
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_01.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_11.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_21.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_31.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_41.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_51.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_61.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_71.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_81.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_91.json"
    };
    vector<string> batch5 = {
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_00.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_10.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_20.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_30.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_40.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_50.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_60.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_70.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_80.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_90.json"
    };

    // 先并行跑第一批 10 个，等全部完成

      run_batch(batch1);

//
//    // 再并行跑第二批 10 个，等全部完成
//    for (int i=0;i<30;i++){
//        run_batch(batch2);
////    }
//    run_batch(batch3);
//    run_batch(batch4);
//    run_batch(batch5);

//
//   for(int i=0; i<1;i++){
//       Test("D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//class_08_14.json",21+i);
//   }
//    for(int i=0; i<10;i++){
//    Test("test//E0级颗粒板E0级暖白色_18.json",21+i);
//}
//    for(int i=0; i<10;i++){
//        Test("test//E1级中纤板双面贴灰色三胺（森荣3131-J1123浅灰）_18.json",21+i);
//    }
//    for(int i=0; i<10;i++){
//        Test("test//E1级防潮板双面贴黑色XJ-951绒麻面三聚氰胺_16.json",21+i);
//    }
//    for(int i=0; i<10;i++){
//        Test("test//多层板G6063-311_18.json",21+i);
//    }
//    for(int i=0; i<10;i++){
//        Test("test//多层板双白胶板_19.json",21+i);
//    }
//    for(int i=0; i<10;i++){
//        Test("test//多层板安科纳胡桃_18.json",21+i);
//    }
//    for(int i=0; i<10;i++){
//        Test("test//实木多层板暖白麻面_18.json",21+i);
//    }
//    for(int i=0; i<10;i++){
//        Test("test//实木颗粒板珊瑚灰木纹_18.json",21+i);
//     }
//     for(int i=0; i<10;i++){
//         Test("test//实木颗粒板箭羽棕_18.json",21+i);
//     }
//    for(int i=0; i<10;i++){
//        Test("test//暖白多层_18.json",21+i);
//    }
//    for(int i=0; i<10;i++){
//        Test("test//进口杉木实芯板和信白_18.json",21+i);
//    }
//    for(int i=0; i<10;i++){
//        Test("tests// 18澳松板玄铁灰.json",21+i);
//    }
//    Run10Parallel("test//E0级颗粒板E0级暖白色_18.json");
//    Run10Parallel("test//E1级中纤板双面贴灰色三胺（森荣3131-J1123浅灰）_18.json");
//    Run10Parallel("test//E1级防潮板双面贴黑色XJ-951绒麻面三聚氰胺_16.json");
//    Run10Parallel("test//多层板G6063-311_18.json");
//    Run10Parallel("test//多层板双白胶板_19.json");
//    Run10Parallel("test//多层板安科纳胡桃_18.json");
//    Run10Parallel("test//实木多层板暖白麻面_18.json");
//    Run10Parallel("test//实木颗粒板珊瑚灰木纹_18.json");
//    Run10Parallel("test//实木颗粒板箭羽棕_18.json");
//    Run10Parallel("test//暖白多层_18.json");
//    Run10Parallel("test//进口杉木实芯板和信白_18.json");
//    for (int i = 1; i <= 30; ++i) {
//        char fname[16];
//        // 生成 TB004.json, TB004.json, TB004.json
//        std::sprintf(fname, "TN%03d.json", i);
//        calnfp(fname);
//    }
//    for(int i=0; i<1;i++){
//        Test("fu.json",21+i);
//    }
//    int value = 0;  // 比如计算出来的结果
//    std::cout << value << std::endl;
    return 0;
}