# Reproducibility Package (Statistical Testing Data)

This repository provides the data used for statistical analysis in our revision.

## 1) Build Environment (Code Compilation)

The C++ implementation is compiled under **Cygwin**.  
Please ensure a working Cygwin environment is installed, including a C++ toolchain (e.g., `g++`) and common build utilities (e.g., `make`/`cmake`, if applicable).

> Note: This package focuses on the reproducibility data for statistical testing. The compilation requirement is listed here to document the environment used in our experiments.

## 2) Released Results for Kruskal–Wallis (KW) Testing

We release the **final per-run objective results** used in the **Kruskal–Wallis (KW)** test.

- Each **line** in the provided results file corresponds to **one run result** (one value per line).
- The file contains results for three methods in the following fixed order:
  1. **Lines 1–360:** `RC-PLNS`
  2. **Lines 361–720:** `BRKGA`
  3. **Lines 721–1080:** `BRKGA+BSSF`

These values can be directly loaded to reproduce the KW test and subsequent post hoc comparisons (e.g., Dunn test with Holm correction), as reported in the manuscript.

## Contact

For questions regarding the released data format or statistical testing, please contact the corresponding author.
