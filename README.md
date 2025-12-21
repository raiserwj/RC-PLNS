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
# Code Usage Guide

## 1. Current Version
The current codebase is configured for **testing the KW test**.

## 2. How to Run
- **Main entry**: `amopt_pack.cpp`
- Place the `data` files into the **specified directory**, then run the program.
- **Output**: results will be written to `record_.txt` under the working directory.

## 3. Configuration for Other Experiments

### 3.1 Modify Input Attributes (Priority / Rotation Angle)
In `carvingmachine.cpp`:
- **Lines 98–102**: in the function inputs,
  - the **8th argument** corresponds to the **priority attribute**,
  - the **10th argument** corresponds to the **rotation-angle attribute**.
- You can adjust the problem input by modifying these two arguments.

### 3.2 Enable the Combination Module
In `carvingmachine.cpp`:
- **Line 379**: the `combination` function is the **combination module**.
- Uncomment it to enable this module.

### 3.3 Control Time Limits in the Packing Module
In `carvingmachine.cpp`:
- **Line 552**: controls the time limit for the **priority stage** in the `packing` module.
- **Lines 813, 829, 830**: control the time limits for the **ordinary stage** in the `packing` module.

### 3.4 Enable the Reproduced BRKGA Algorithm
In `amopt_pack.cpp`:
- **Line 24**: uncomment to enable the **reproduced BRKGA** algorithm.
