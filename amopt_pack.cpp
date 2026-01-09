#include <sys/wait.h>
#include"amopt_pack.h"
#include <thread>
#include "unistd.h"
using namespace std;
Error Compute(string str, Json::Value &result_list, string filename) {
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

    compute(result_list,str);
//        carvingmachine.compute(result_list);

    return amopt::COMPUTE_NO_ERROR;
}


void Test(string filename,int filenum) {
    string s;
    std::stringstream str;
    std::fstream f;
    f.open(filename, std::ios::in);
    if (!f.is_open()) {
        std::cout << "Open json file error!" << std::endl;
        exit(0);
    }
    str << f.rdbuf();
    f.close();
    s = str.str();
    Json::Value result_list;
    amopt::Error err = Compute(s, result_list, filename);
    if (err != amopt::Error::COMPUTE_NO_ERROR) {
        std::cout << "error" << err << std::endl;
    }
    ofstream fout;
    filename.erase(filename.length() - 5);
    filename.erase(0,6);
    fout.open("sol//"+filename+"_sol"+ std::to_string(filenum));
    Json::StyledWriter writer1;
//    fout << writer1.write(result_list) << std::endl;
    fout.close();
    ofstream of;
    of.open(filename+".txt",std::ios::out );
    of<<result_list["solutions"].size()<<std::endl;
    of.close();
}
void run_batch(const vector<string>& files) {
    // 启动一批子进程
    const int TIME_LIMIT_SEC = 120;  // 整批这 10 个实例允许的最大时间

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
int main() {
    Json::Value result_list;

    vector<string> batch1 = {
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_04.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_14.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_24.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_34.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_44.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_54.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_64.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_74.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_84.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_94.json",
    };

    vector<string> batch2 = {
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_03.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_13.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_23.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_33.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_43.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_53.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_63.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_73.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_83.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_93.json",
    };

    vector<string> batch3 = {
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_02.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_12.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_22.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_32.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_42.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_52.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_62.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_72.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_82.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_92.json",
    };
    vector<string> batch4 = {
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_01.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_11.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_21.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_31.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_41.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_51.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_61.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_71.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_81.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_91.json",
    };
    vector<string> batch5 = {
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_00.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_10.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_20.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_30.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_40.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_50.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_60.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_70.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_80.json",
            "D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_90.json",
    };

    // 先并行跑第一批 10 个，等全部完成
    run_batch(batch1);

    // 再并行跑第二批 10 个，等全部完成
    run_batch(batch2);
//    run_batch(batch3);
//    run_batch(batch4);
//    run_batch(batch5);

//
//   for(int i=0; i<1;i++){
//       Test("D://123456//LLMLNS//220914-CarvingMachine//cmake-build-release//2bp//Class_10_14.json",21+i);
//   }
//    for(int i=0; i<1;i++){
//    Test("test//E0级颗粒板E0级暖白色_18.json",21+i);
//}
//    for(int i=0; i<1;i++){
//        Test("test//E1级中纤板双面贴灰色三胺（森荣3131-J1123浅灰）_18.json",21+i);
//    }
//    for(int i=0; i<1;i++){
//        Test("test//E1级防潮板双面贴黑色XJ-951绒麻面三聚氰胺_16.json",21+i);
//    }
//    for(int i=0; i<1;i++){
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
//    for(int i=0; i<1;i++){
//        Test("tests// 18澳松板玄铁灰.json",21+i);
//    }
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
