# 🧬 In-Silico Prediction of Eukaryotic Secretory Signal Peptides

> **University of Bologna (Alma Mater Studiorum)**  
> *Department of Pharmacy and Biotechnology (FaBiT)*  
> **Author:** Mohammadamin Roshan (MSc Student in Bioinformatics)  
> **Status:** DOME-Compliant Machine Learning Pipeline & Peer-Reviewed Manuscript Ready

---

## 🔍 Project Overview

Accurate computational detection of eukaryotic secretory signal peptides (SPs) is essential for genome-scale functional annotation, but it is heavily hindered by **Transmembrane (TM) mimicry**—the physicochemical similarity between SP hydrophobic cores and N-terminal TM $\alpha$-helices.

To systematically address this challenge, this repository implements and evaluates **three paradigms of increasing representational complexity**:
1. **Rule-Based Motif Scanning:** The classic von Heijne Position-Specific Weight Matrix (PSWM).
2. **Feature-Engineered Machine Learning:** An RBF-SVM trained on 35 extracted physicochemical descriptors, rigorously optimized via Nested Cross-Validation.
3. **Deep Representation Learning:** Deep MLP and hybrid PCA-SVM architectures leveraging contextual embeddings from the pretrained **ESM-2 Protein Language Model (PLM)**.

---

## 🛠️ Pipeline Structure

The project is organized in a sequential, reproducible workflow spanning from raw data curation to advanced statistical diagnostics:

* `step1_Data_Collection/`: Automated REST API retrieval from UniProtKB.
* `step2_data_preparation/`: MMseqs2 clustering and stratified 80/20 partitioning.
* `step3_data_analysis/`: Comprehensive Exploratory Data Analysis (EDA) & Sequence Logos.
* `step4_von_heijne/`: PSWM baseline implementation and threshold optimization.
* `step5_svm/`: Nested CV feature selection (RF Gini) and SVM training.
* `step6_deep_learning/`: ESM-2 feature extraction, PCA, and neural network training.
* `step7_bologna_advanced_pipeline/`: Dual-thresholding, RI calibration, Viterbi decoding, and exact McNemar's statistical testing.

---

## 📊 Key Results & Visual Diagnostics

On the completely held-out blind benchmarking set ($N = 2,037$), our models achieved significant improvements across evaluation metrics:

* **von Heijne PSWM:** $\text{MCC} = 0.655$
* **Physicochemical SVM:** $\text{MCC} = 0.859$
* **Hybrid ESM-2 PLM:** $\text{MCC} = 0.943$ | $\text{ROC-AUC} = 99.37\%$

### 1. Master Model Benchmark Comparison
The performance evaluation across different model architectures on the test set:

![Master Model Benchmark Comparison](figures/fig6_1_master_model_comparison.png)

### 2. Feature Selection Diagnostics (SVM)
Recursive feature elimination and nested cross-validation feature stability:

![Nested CV Feature Selection](figures/fig5_1_nested_cv_feature_selection.png)

### 3. Advanced Pipeline & Error Analysis
Statistical error diagnostics (Fisher's Exact Test) proved that ESM-2 representations reduced N-terminal TM-helix false-positive rates dramatically compared to traditional baseline scoring.

![Bologna Advanced Pipeline Analysis](figures/fig7_1_bologna_advanced_analysis.png)

---

## 📄 Final Manuscript

The complete peer-reviewed manuscript outlining the methodology, benchmark results, and DOME compliance checklist has been finalized and structured for submission to *Oxford Bioinformatics*:
* **Primary Document:** `13.docx` (Stored in root repository / version-controlled).

---

## 📬 Author & Contact

**Mohammadamin Roshan**  
*MSc Student in Bioinformatics | Alma Mater Studiorum – Università di Bologna*  
* [LinkedIn Profile](https://www.linkedin.com/in/mohammadaminroshan/)  
* [GitHub Repository](https://github.com/mohammadaminroshan/LB2_Project)