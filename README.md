# 🧬 In-Silico Prediction of Eukaryotic Secretory Signal Peptides

> **University of Bologna (Alma Mater Studiorum)**  
> *Department of Pharmacy and Biotechnology (FaBiT)*  
> **Author:** Mohammadamin Roshan (MSc Student in Bioinformatics)  

---

## 📌 Project Overview

In this project, we worked on a crucial problem in bioinformatics: **accurate identification of eukaryotic secretory signal peptides (SPs)**. 
Tackling this problem comes with a major challenge known as **Transmembrane (TM) mimicry**—the physicochemical similarity between the hydrophobic cores of signal peptides and N-terminal TM $\alpha$-helices, which complicates computational detection.

To systematically address this challenge, we implemented and evaluated **three paradigms of increasing representational complexity**:
1. **Rule-Based Motif Scanning:** The classic von Heijne Position-Specific Weight Matrix (PSWM).
2. **Feature-Engineered Machine Learning:** An RBF-SVM trained on 35 extracted physicochemical descriptors, optimized via Nested Cross-Validation.
3. **Deep Representation Learning:** Deep MLP and hybrid PCA-SVM architectures leveraging contextual embeddings from the pretrained **ESM-2 Protein Language Model**.

---

## 🗂️ Pipeline Structure

The project is organized in a clean, reproducible workflow from step 1 to step 7:

* `step1_Data_Collection/`: Automated data retrieval via UniProtKB REST API.
* `step2_data_preparation/`: MMseqs2 clustering and stratified 80/20 partitioning.
* `step3_data_analysis/`: Comprehensive Exploratory Data Analysis (EDA) and Sequence Logos.
* `step4_von_heijne/`: Baseline PSWM implementation and threshold optimization.
* `step5_svm/`: Feature selection (RF Gini) and SVM training.
* `step6_deep_learning/`: ESM-2 feature extraction, PCA dimensionality reduction, and neural network training.
* `step7_bologna_advanced_pipeline/`: Dual-thresholding, RI calibration, Viterbi decoding, and exact McNemar's statistical testing.

---

## 📊 Key Results & Visual Diagnostics

On the completely held-out blind benchmarking set ($N = 2,037$), our models achieved the following performance:
* **von Heijne PSWM Baseline:** $\text{MCC} = 0.655$
* **Physicochemical SVM:** $\text{MCC} = 0.859$
* **Hybrid ESM-2 PLM:** $\text{MCC} = 0.943$ | $\text{ROC-AUC} = 99.37\%$

### 1. Master Model Benchmark Comparison
![Master Model Benchmark Comparison](step6_deep_learning/figures/fig6_1_master_model_comparison.png)

### 2. SVM Feature Selection Diagnostics
![Nested CV Feature Selection](step5_svm/figures/fig5_1_nested_cv_feature_selection.png)

### 3. Advanced Pipeline & Error Analysis
![Bologna Advanced Pipeline Analysis](step7_bologna_advanced_pipeline/figures/fig7_1_bologna_advanced_analysis.png)

---

## 📄 Final Manuscript

* **Version 13 (`13.docx`):** The final peer-reviewed manuscript draft prepared for submission to *Oxford Bioinformatics* is stored within the repository.

---

## 📬 Author

**Mohammadamin Roshan**  
*MSc Student in Bioinformatics, University of Bologna*  
[LinkedIn Profile](https://www.linkedin.com/in/mohammadaminroshan/)