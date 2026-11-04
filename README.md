# In-Silico Prediction of Eukaryotic Secretory Signal Peptides 🧬

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![DOME Compliant](https://img.shields.io/badge/DOME-Compliant-brightgreen.svg)](https://www.nature.com/articles/s41592-021-01205-4)

This repository contains an end-to-end, **DOME-compliant** machine learning pipeline for the detection of eukaryotic secretory signal peptides (SPs) and the precise localization of their cleavage sites. 

This project was developed at the **Department of Pharmacy and Biotechnology (FaBiT), Alma Mater Studiorum – Università di Bologna**.

## 🚀 Project Overview

Accurate computational detection of SPs is essential for genome-scale functional annotation. However, it is heavily hindered by **Transmembrane (TM) mimicry** (the physicochemical similarity between SP hydrophobic cores and N-terminal TM $\alpha$-helices). 

In this study, we systematically implemented and compared three paradigms of increasing representational complexity:
1. **Rule-Based Motif Scanning:** The classic von Heijne Position-Specific Weight Matrix (PSWM).
2. **Feature-Engineered Machine Learning:** An RBF-SVM trained on 35 extracted physicochemical descriptors optimized via **Nested Cross-Validation**.
3. **Deep Representation Learning:** Deep MLP and Hybrid PCA-SVM architectures leveraging contextual embeddings from the pretrained **ESM-2** Protein Language Model (PLM).

## 📂 Pipeline Structure 

The pipeline is sequentially structured into the following steps:
* `step1_Data_Collection/`: Automated REST API retrieval from UniProtKB.
* `step2_data_preparation/`: MMseqs2 clustering and stratified 80/20 partitioning.
* `step3_data_analysis/`: Comprehensive Exploratory Data Analysis (EDA) & Sequence Logos.
* `step4_von_heijne/`: PSWM baseline implementation and threshold optimization.
* `step5_svm/`: Nested CV feature selection (RF Gini) and SVM training.
* `step6_deep_learning/`: ESM-2 feature extraction, PCA, and neural network training.
* `step7_bologna_advanced_pipeline/`: Dual-thresholding, RI calibration, Viterbi decoding, and exact McNemar's statistical testing.

## 📊 Key Results
On the completely held-out blind benchmarking set (N=2,037):
* **von Heijne PSWM:** MCC = 0.655
* **Physicochemical SVM:** MCC = 0.859
* **Hybrid ESM-2 PLM:** MCC = 0.943 | ROC-AUC = 99.37%

*Statistical error diagnostics (Fisher's Exact Test) proved that ESM-2 representations reduced N-terminal TM-helix false-positive rates from 60.6% (PSWM) down to 25.0%.*

## 📝 Final Manuscript
The final peer-reviewed manuscript (Version 13) has been generated outlining the complete methodology, results, and DOME compliance of this study.

## 👨‍💻 Author
**Mohammadamin Roshan**
* MSc Student in Bioinformatics, University of Bologna
* [LinkedIn](https://www.linkedin.com/in/mohammadaminroshan/) | [GitHub](https://github.com/mohammadaminroshan)