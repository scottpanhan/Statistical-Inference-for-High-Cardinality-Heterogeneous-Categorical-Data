# Codes for "Power-enhanced Testing and Conformalized-adaptive Selection for High-Cardinality Heterogeneous Categorical Data"

This repository contains the implementation code for the paper:  
**"Power-enhanced Testing and Conformalized-adaptive Selection for High-Cardinality Heterogeneous Categorical Data"**  
by Han Pan, Wei Xiong, and Youcheng Yu.

## 📁 `./codes/`

This directory includes the core scripts to reproduce the proposed methods and their competitors. Each file contains an independent simulation example for demonstration.

*   **`G-sum test for fixed R.R`**: Implements the G-sum test and the VM test. It also includes code for competitor methods: DC, HHG, SD, MINT, HSIC, and MV tests.
*   **`category-specific test for fixed R.R`**: Implements our proposed category-specific test.
*   **`power-enhanced test when R diverges`**: Implements our proposed power-enhanced test alongside the competitor SD test (both permutation and asymptotic versions).
*   **`FDR control with CASH procedure`**: Executes the proposed CASH procedure and compares it with established methods including Knockoff, BH, IHW, and AdaDetect.
