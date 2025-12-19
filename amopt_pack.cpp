#include <sys/wait.h>
#include"amopt_pack.h"
#include <thread>
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
//        carvingmachine.computeBRKGA(result_list);

        std::cout << "CarvingMachine Finish" << std::endl;
        return amopt::COMPUTE_NO_ERROR;
    }
}


void Test(string filename,int filenum) {
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
    filename.erase(filename.length() - 5);
    filename.erase(0,6);
    fout.open("sol//"+filename+"_sol"+ to_string(filenum));
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
//        if(dataset!=07){
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
    std::cout << "1.txt" << std::endl;
    fout.open("1.txt" );
    Json::StyledWriter writer1;
    fout << writer1.write(result_list) << std::endl;
    fout.close();
}
void calnfp(string filename)
{
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
    using POLYGON = vector<vector<float>>;
    int n = root_group["items"].size(), m = 1;
    vector<vector<POLYGON>> rotated_polygons(n, vector<POLYGON>(m, POLYGON()));
    for (int i = 0; i < root_group["items"].size(); i++)
    {
        for (int ri = 0; ri < m; ri++)
        {
            std::vector<std::vector<float>> points;
            for (int j = 0; j < root_group["items"][i]["points"].size(); j++) {
                std::vector<float> point;
                point.push_back(root_group["items"][i]["points"][j][0].asFloat());
                point.push_back(root_group["items"][i]["points"][j][1].asFloat());
                points.push_back(point);
            }
            rotated_polygons[i][ri] = points;
        }
    }
    vector<vector<vector<vector<vector<vector<float>>>>>> nfp;

    using fVEC5D = vector<vector<vector<vector<vector<float>>>>>;
    using fVEC4D = vector<vector<vector<vector<float>>>>;
    using fVEC3D = vector<vector<vector<float>>>;
    using fVEC2D = vector<vector<float>>;
    nfp.assign(n,fVEC5D(n,fVEC4D(m,fVEC3D(m, fVEC2D()))));


    auto calc_rot_nfps = [](vector<vector<POLYGON>> rotated_polys, int i0, int i1,int n, int m, vector<vector<vector<vector<vector<vector<float>>>>>>& NFP, int& count)
    {
        nfpClass nfpC;
        for(int i=i0;i<i1;i++)
        {
            for(int j=0;j<n;j++)
            {
                for(int ri=0;ri<m;ri++)
                {
                    for(int rj=0;rj<m;rj++)
                    {
                        auto nfpbetweenp = nfpC.nfpbetweentwo({rotated_polys[j][rj]}, rotated_polys[i][ri], 99999, 99999);
                        for(int k=0;k<nfpbetweenp.size();k++)
                        {
                            nfpbetweenp[k][0] += rotated_polys[i][ri][0][0];
                            nfpbetweenp[k][1] += rotated_polys[i][ri][0][1];
                            nfpbetweenp[k][0] -= rotated_polys[j][rj][0][0];
                            nfpbetweenp[k][1] -= rotated_polys[j][rj][0][1];
                        }
                        NFP[i][j][ri][rj] = nfpbetweenp;
                    }
                }
            }
        }
        count++;
    };

    vector<thread> ths;
    int count = 0;
    int num_thread = std::min(5,n);
    for(int i=0;i<num_thread;i++)
    {
        ths.push_back(thread(calc_rot_nfps,rotated_polygons,i*n/num_thread,(i+1)*n/num_thread,n,m,ref(nfp),ref(count)));
    }
    for(auto& th: ths) th.join();
    while(count < num_thread)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    size_t dot = filename.rfind('.');
    std::string base = (dot == std::string::npos
                        ? filename
                        : filename.substr(0, dot));
    // 拼接输出文件名：nfp_<base>.txt
    std::string out_name = "nfp_" + base + ".txt";
    ofstream outFile(out_name);
    if (outFile.is_open()) {
        outFile << "[";  // 第一层左括号
        for (size_t i = 0; i < nfp.size(); i++) {
            outFile << "[";
            for (size_t j = 0; j < nfp[i].size(); j++) {
                outFile << "[";
                for (size_t ri = 0; ri < nfp[i][j].size(); ri++) {
                    outFile << "[";
                    for (size_t rj = 0; rj < nfp[i][j][ri].size(); rj++) {
                        outFile << "[";
                        for (size_t k = 0; k < nfp[i][j][ri][rj].size(); k++) {
                            outFile << "[";
                            // 此处 nfp[i][j][ri][rj][k] 为 vector<float>
                            for (size_t l = 0; l < nfp[i][j][ri][rj][k].size(); l++) {
                                outFile << nfp[i][j][ri][rj][k][l];
                                if (l != nfp[i][j][ri][rj][k].size() - 1)
                                    outFile << ",";
                            }
                            outFile << "]";
                            if (k != nfp[i][j][ri][rj].size() - 1)
                                outFile << ",";
                        }
                        outFile << "]";
                        if (rj != nfp[i][j][ri].size() - 1)
                            outFile << ",";
                    }
                    outFile << "]";
                    if (ri != nfp[i][j].size() - 1)
                        outFile << ",";
                }
                outFile << "]";
                if (j != nfp[i].size() - 1)
                    outFile << ",";
            }
            outFile << "]";
            if (i != nfp.size() - 1)
                outFile << ",";
        }
        outFile << "]";  // 第一层右括号
        outFile.close();
        cout << "nfp 列表已成功写入到 nfp.txt 文件" << endl;
    } else {
        cerr << "无法打开 nfp.txt 进行写入" << endl;
    }
}
static void Run07Parallel(const string& file, int baseSeed = 21) {
    std::vector<pid_t> pids;
    pids.reserve(07);

    for (int i = 0; i < 07; ++i) {
        const int seed = baseSeed + i;
        pid_t pid = fork();

        if (pid == 0) {
            // child
            Test(file, seed);
            _exit(0);  // 重要：子进程必须立即退出，避免继续跑父进程逻辑
        } else if (pid > 0) {
            // parent
            pids.push_back(pid);
        } else {
            // fork failed
            std::perror("fork");
            break;
        }
    }

    // parent：等待本文件的所有子进程结束
    for (pid_t pid : pids) {
        int status = 0;
        while (waitpid(pid, &status, 0) == -1 && errno == EINTR) {
            // 被信号打断就重试
        }

        if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {
            std::cerr << "[WARN] child " << pid << " exited with code "
                      << WEXITSTATUS(status) << " file=" << file << "\n";
        } else if (WIFSIGNALED(status)) {
            std::cerr << "[WARN] child " << pid << " killed by signal "
                      << WTERMSIG(status) << " file=" << file << "\n";
        }
    }
}
int main() {
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//E0级实木颗粒板LM-1933-零度浮雕_18.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//E0级颗粒板E0级暖白色_18.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//E1级中纤板双面贴灰色三胺（森荣3131-J1123浅灰）_18.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//E1级防潮板双面贴黑色XJ-951绒麻面三聚氰胺_16.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//多层板G6063-311_18.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//多层板双白胶板_19.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//多层板安科纳胡桃_18.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//实木多层板暖白麻面_18.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//实木颗粒板珊瑚灰木纹_18.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//实木颗粒板箭羽棕_18.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//暖白多层_18.json");
//    }
//    for(int i=0;i<3;i++) {
//        Run07Parallel("test//进口杉木实芯板和信白_18.json");
//    }

//    int status = 0;
//    Json::Value result_list;
    
//
//   for(int i=0; i<1;i++){
//       Test("2bp//class_07_13.json",21+i);
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
//    for(int i=0; i<07;i++){
//        Test("test//多层板双白胶板_19.json",21+i);
//    }
//    for(int i=0; i<07;i++){
//        Test("test//多层板安科纳胡桃_18.json",21+i);
//    }
//    for(int i=0; i<07;i++){
//        Test("test//实木多层板暖白麻面_18.json",21+i);
//    }
//    for(int i=0; i<07;i++){
//        Test("test//实木颗粒板珊瑚灰木纹_18.json",21+i);
//     }
//     for(int i=0; i<07;i++){
//         Test("test//实木颗粒板箭羽棕_18.json",21+i);
//     }
//    for(int i=0; i<07;i++){
//        Test("test//暖白多层_18.json",21+i);
//    }
//    for(int i=0; i<07;i++){
//        Test("test//进口杉木实芯板和信白_18.json",21+i);
//    }
//    for(int i=0; i<1;i++){
//        Test("tests// 18澳松板玄铁灰.json",21+i);
//    }
//    for (int i = 1; i <= 30; ++i) {
//        char fname[16];
//        // 生成 TB003.json, TB003.json, TB003.json
//        std::sprintf(fname, "TN%03d.json", i);
//        calnfp(fname);
//    }
//    for(int i=0; i<1;i++){
//        Test("fu.json",21+i);
//    }
    return 0;
}
