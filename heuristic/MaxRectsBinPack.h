/** @file MaxRectsBinPack.h
	@author Jukka Jylnki

	@brief Implements different bin packer algorithms that use the MAXRECTS data structure.

	This work is released to Public Domain, do whatever you want with it.
*/
#pragma once

#include <vector>

#include "Rect.h"

namespace rbp {

/** MaxRectsBinPack implements the MAXRECTS data structure and different bin packing algorithms that 
	use this structure. */
    class MaxRectsBinPack {
    public:
        /// Instantiates a bin of size (0,0). Call Init to create a new bin.
        MaxRectsBinPack();

        /// Instantiates a bin of the given size.
        /// @param allowFlip Specifies whether the packing algorithm is allowed to rotate the input rectangles by 90 degrees to consider a better placement.
        MaxRectsBinPack(double width, double height, bool allowFlip = true);

        /// (Re)initializes the packer to an empty bin of width x height units. Call whenever
        /// you need to restart with a new bin.
        void Init(double width, double height, bool allowFlip = true);

        /// Specifies the different heuristic rules that can be used when deciding where to place a new rectangle.
        enum FreeRectChoiceHeuristic {
            RectBestShortSideFit, ///< -BSSF: Positions the rectangle against the short side of a free rectangle into which it fits the best.
            RectBestLongSideFit, ///< -BLSF: Positions the rectangle against the long side of a free rectangle into which it fits the best.
            RectBestAreaFit, ///< -BAF: Positions the rectangle into the smallest free rect into which it fits.
            RectwidthRule,
            RectheightRule,///< -BL: Does the Tetris placement.
            RectContactPointRule,
            RectDFTRC///< -CP: Choosest the placement where the rectangle touches other rects as much as possible.
        };

        /// Inserts the given list of rectangles in an offline/batch mode, possibly rotated.
        /// @param rects The list of rectangles to insert. This vector will be destroyed in the process.
        /// @param dst [out] This list will contain the packed rectangles. The indices will not correspond to that of rects.
        /// @param method The rectangle placement rule to use when packing.
        double Insert(std::vector<std::vector<double>> rects, FreeRectChoiceHeuristic method);

        /// Inserts a single rectangle into the bin, possibly rotated.
        Rect Insert(double width, double height, int Flip, FreeRectChoiceHeuristic method);

        double Insertscore(double width, double height, int Flip, FreeRectChoiceHeuristic method);

        /// Computes the ratio of used surface area to the total bin area.
        double Occupancy() const;

        std::vector<std::vector<double> > return_location();

        void free();

    private:
        double binWidth;
        double binHeight;

        bool binAllowFlip;

        std::vector<Rect> usedRectangles;
        std::vector<Rect> freeRectangles;

        /// Computes the placement score for placing the given rectangle with the given method.
        /// @param score1 [out] The primary placement score will be outputted here.
        /// @param score2 [out] The secondary placement score will be outputted here. This isu sed to break ties.
        /// @return This struct identifies where the rectangle would be placed if it were placed.
        Rect ScoreRect(double width, double height, int Flip, FreeRectChoiceHeuristic method, double &score1,
                       double &score2) const;

        /// Places the given rectangle into the bin.
        void PlaceRect(const Rect &node);

        /// Computes the placement score for the -CP variant.
        double ContactPointScoreNode(double x, double y, double width, double height) const;

        Rect FindPositionForNewNodewidth(double width, double height, int Flip,
                                         double &bestY, double &bestX) const;

        Rect FindPositionForNewNodeheight(double width, double height, int Flip,
                                          double &bestY, double &bestX) const;

        Rect
        FindPositionForNewNodeBestShortSideFit(double width, double height, int Flip, double &bestShortSideFit,
                                               double &bestLongSideFit) const;

        Rect FindPositionForNewNodeBestLongSideFit(double width, double height, double &bestShortSideFit,
                                                   double &bestLongSideFit) const;

        Rect FindPositionForNewNodeBestAreaFit(double width, double height, double &bestAreaFit,
                                               double &bestShortSideFit) const;

        Rect FindPositionForNewNodeContactPoint(double width, double height, double &contactScore) const;

        Rect FindPositionForNewNodeDFTRC(double width, double height, double &bestScore) const;

        /// @return True if the free node was split.
        bool SplitFreeNode(Rect freeNode, const Rect &usedNode);

        /// Goes through the free rectangle list and removes any redundant entries.
        void PruneFreeList();
    };

}
