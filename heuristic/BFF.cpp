#include "BFF.h"
#include <limits>
std::vector<std::vector<rbp::Rect> > BFF::pack(std::vector<std::vector<double> > rects, double bin_width, double bin_height, int nBin_min,int rotate) {
    std::vector<rbp::GuillotineBinPack> detailpack = {};
    std::vector<std::vector<rbp::Rect>> detail = {};
    for (int nBin = nBin_min; nBin < 300; nBin++) {
        for (int i = 0; i < nBin; i++) {
            detailpack.push_back(GuillotineBinPack(bin_width, bin_height));
        }

        std::vector<std::vector<double> > rects_copy;
        rects_copy.assign(rects.begin(), rects.end());

        int nr = std::end(rects_copy) - std::begin(rects_copy);
        bool flag = true;
        for (int i = 0; i < nr; i++) {
            double bestScore1 = std::numeric_limits<double>::max();
            int bestRectIndex = 0;
            int bestBinIndex = 0;
            for (size_t i1 = 0; i1 < rects_copy.size(); ++i1) {
                for (int j = 0; j < nBin; j++) {
                    double score1 = detailpack[j].score(rects_copy[i1][0], rects_copy[i1][1],rotate, GuillotineBinPack::RectBestShortSideFit);
                    if (score1 < bestScore1) {
                        bestScore1 = score1;
                        bestRectIndex = i1;
                        bestBinIndex = j;
                    }
                }
            }
            if (bestScore1 == std::numeric_limits<double>::max()) {
                flag = false;
                break;
            } else {
                detailpack[bestBinIndex].Insert(rects_copy[bestRectIndex][0], rects_copy[bestRectIndex][1],rotate, int(rects_copy[bestRectIndex][2]), false,
                                                GuillotineBinPack::RectBestShortSideFit,
                                                GuillotineBinPack::SplitShorterAxis);
            }
            rects_copy.erase(rects_copy.begin() + bestRectIndex);
        }

        if (flag) {
            num_qiege = 0;
            for (int i = 0; i < nBin; i++) {
                detail.push_back(detailpack[i].GetUsedRectangles());
                num_qiege += detailpack[i].num_qiege;
            }
            num_bins = nBin;
            for (int i = 0; i < nBin; i++) {
                detailpack[i].free();
            }
            return detail;
        }
        for (int i = 0; i < nBin; i++) {
            detailpack[i].free();
        }
        detailpack.clear();
    }
	std::cout << "shibai" << std::endl;
    return detail;
}
