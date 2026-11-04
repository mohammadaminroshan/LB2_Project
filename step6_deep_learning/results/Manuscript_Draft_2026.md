# In-Silico Prediction of Eukaryotic Secretory Signal Peptides: Comparative Evaluation of Position-Specific Weight Matrices, Support Vector Machines, and Protein Language Models

**Authors:** Mohammadamin Roshan et al.  
**Affiliation:** Department of Pharmacy and Biotechnology (FaBiT), University of Bologna, Bologna, Italy  
**Course:** Laboratory of Bioinformatics II (Module 2) — A.Y. 2026–2027  

---

## Abstract

**Motivation:** Secretory signal peptides (SPs) are short N-terminal amino acid sequences that direct nascent polypeptides toward the endoplasmic reticulum and the secretory pathway. Because experimental subcellular localization via mass spectrometry or immunofluorescence microscopy is costly and labor-intensive, accurate computational detection of SPs—and their discrimination from hydrophobic N-terminal transmembrane (TM) $\alpha$-helices—is essential for genome-scale functional annotation.  
**Results:** We constructed a high-quality, non-redundant dataset of eukaryotic proteins from UniProtKB/Swiss-Prot (10,182 sequences after MMseqs2 clustering at 30% sequence identity and 40% alignment coverage), partitioned into a 5-fold cross-validation training set (8,145 sequences) and a held-out benchmarking set (2,037 sequences). We implemented and compared three paradigms: (i) the classic von Heijne position-specific weight matrix (PSWM) scanning the $[-13, +2]$ cleavage-site motif, (ii) a feature-engineered Support Vector Machine (SVM) coupled with nested cross-validation feature selection, and (iii) a deep learning classifier leveraging pretrained evolutionary representations from the ESM-2 protein language model (PLM). On the blind benchmarking set, the von Heijne baseline achieved $\text{MCC} = 0.655$ ($\text{F1} = 0.693$), whereas the Random-Forest-selected RBF-SVM significantly improved generalization to $\text{MCC} = 0.859$ ($\text{F1} = 0.873$). The ESM-2 PLM extension further elevated performance to $\text{MCC} = 0.935$ ($\text{ROC-AUC} = 0.996$). Detailed error analysis revealed that 60.6% of von Heijne false positives and 65.0% of SVM false positives harbour an N-terminal TM helix within the first 90 residues, confirming hydrophobicity mimicry as the primary source of false detections.  
**Availability:** All Python notebooks, trained models, and datasets are reproducible and hosted on the group GitHub repository.  
**Contact:** mohammadamin.roshan@studio.unibo.it  

---

## 1 Introduction
Subcellular localization defines the biochemical environment in which a protein operates. In eukaryotes, proteins destined for the endoplasmic reticulum (ER), Golgi apparatus, plasma membrane, or extracellular space typically carry an N-terminal secretory signal peptide (SP) cleaved by Signal Peptidase I (Sec/SPI) upon translocation. Canonical SPs range from 15 to 30 residues and share a tripartite architecture: a positively charged **n-region**, a hydrophobic $\alpha$-helical **h-region** (6–15 residues), and a polar **c-region** terminating in a weakly conserved cleavage motif ($A\text{-}x\text{-}A$ at positions $-3$ and $-1$). Despite this conserved organization, in-silico SP detection remains challenging due to sequence variability across eukaryotic kingdoms and the strong physicochemical resemblance between the SP h-region and N-terminal transmembrane (TM) helices. Here, we design an end-to-end, leakage-free machine-learning pipeline comparing motif-driven (von Heijne), handcrafted feature-based (SVM), and pretrained transformer-based (ESM-2) predictors.

---

## 2 Methods

### 2.1 Dataset Collection and Redundancy Reduction
Reviewed eukaryotic entries with experimental evidence at the protein level were retrieved via the UniProtKB REST API. Positive sequences required an experimentally verified Signal peptide annotation (`ECO:0000269`), protein length $\ge 40$ residues, known cleavage site, and SP length $\ge 14$ residues. Negative sequences lacked any SP annotation and possessed experimental localization (`ECO:0000269`) in non-secretory compartments (cytosol, nucleus, mitochondrion, plastid, peroxisome, or cell membrane), recording whether a TM helix started within the first 90 N-terminal residues. To prevent data leakage, positive and negative sets were independently clustered using **MMseqs2** (`--min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1`, connected-component clustering). Representatives were randomly split into **80% Training** ($N = 8,145$, stratified into 5 folds) and **20% Benchmarking** ($N = 2,037$, held out strictly for final evaluation).

### 2.2 The von Heijne PSWM Predictor
Following von Heijne (1986), aligned 15-residue cleavage windows spanning $[-13, +2]$ were extracted from training positives. With pseudocounts of 1, the Position-Specific Probability Matrix (PSPM) was computed as:
$$M_{k,j} = \frac{1 + \sum_{i=1}^{N} \mathbb{I}(s_{i,j} = k)}{N + 20}$$
and converted into log-odds weights against SwissProt background frequencies $b_k$:
$$W_{k,j} = \ln \left( \frac{M_{k,j}}{b_k} \right)$$
For any query protein, the 15-residue window was slid across the first 90 N-terminal residues (positions $1 \dots 76$), and the maximum log-odds score was taken as the global SP score. During 5-fold cross-validation, 3 folds trained $W$, 1 validation fold selected the F1-optimal decision threshold $\theta^*$, and 1 fold tested performance. For the final production model, $\theta_{\text{final}} = 6.002$ was fixed on Fold 4 (using $W$ trained on Folds 0–3), and $W$ was retrained on all 5 training folds prior to a single evaluation on the benchmarking set.

### 2.3 Feature Extraction, Nested Cross-Validation, and SVM
Each sequence was encoded into 35 N-terminal descriptors: 20 amino-acid composition frequencies in the first $K = 30$ residues, and 15 physicochemical features extracted via sliding windows over the first 40 residues (Kyte-Doolittle hydropathy, net and positive n-region charge, Zimmerman bulkiness, Chou-Fasman $\alpha$-helix propensity, and Zhao-London TM tendency). To evaluate feature selection without upward bias, we implemented **Nested Cross-Validation** (5 outer test folds $\times$ 4 inner validation folds choosing $k$) across Filter (Variance+Correlation, Mutual Information), Wrapper (RFE), and Embedded ($L_1$ Logistic Regression, Random Forest Gini Importance) selectors. The final production pipeline (`StandardScaler` + Random Forest Top-$k$ Selector + RBF-SVM) was tuned via 5-fold CV (`{'selector__k': 30, 'svm__C': 5.0, 'svm__gamma': 0.01, 'svm__kernel': 'rbf'}`), retrained on the entire training set, and tested once on the benchmarking set.

### 2.4 Deep Learning Extension with ESM-2 Embeddings
To evaluate modern protein language models, we extracted 640-dimensional N-terminal representations (mean + max pooling of layer 6 from `esm2_t6_8M_UR50D` over the first 90 residues) and trained a two-layer Multi-Layer Perceptron (`128 -> 32`, ReLU, $L_2$ regularization $\alpha=0.05$).

---

## 3 Results

### 3.1 Cross-Validation and Benchmarking Performance
Table 1 summarizes the 5-fold cross-validation metrics ($\text{Mean} \pm \text{SE}$) and independent hold-out benchmarking results across all three models.

| Model | Dataset | ACC | Precision (PPV) | Recall (Sen) | F1-Score | MCC | ROC-AUC |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **von Heijne (PSWM)** | 5-Fold CV | 0.931 ± 0.004 | 0.672 ± 0.028 | 0.739 ± 0.021 | 0.701 ± 0.010 | 0.665 ± 0.012 | 0.959 ± 0.004 |
| **von Heijne (PSWM)** | Benchmark | 0.926 | 0.631 | 0.768 | 0.693 | 0.655 | 0.948 |
| **SVM (RF-Selected)** | 5-Fold CV | 0.978 ± 0.001 | 0.902 ± 0.011 | 0.894 ± 0.004 | 0.898 ± 0.006 | 0.886 ± 0.006 | 0.989 ± 0.001 |
| **SVM (RF-Selected)** | Benchmark | 0.973 | 0.903 | 0.845 | 0.873 | 0.859 | 0.991 |
| **ESM-2 + MLP (PLM)** | 5-Fold CV | 0.988 ± 0.001 | 0.960 ± 0.006 | 0.926 ± 0.017 | 0.942 ± 0.007 | 0.936 ± 0.008 | 0.994 ± 0.002 |
| **ESM-2 + MLP (PLM)** | Benchmark | 0.988 | 0.962 | 0.923 | 0.942 | 0.935 | 0.996 |

### 3.2 Error Analysis: Transmembrane Helices and Cleavage-Site Divergence
1. **False Positive (FP) Analysis:** Among the 99 false positives generated by the von Heijne model on the benchmarking set, **60.6%** (60/99) possess an annotated N-terminal transmembrane helix within the first 90 residues (Fisher's exact odds ratio = 17.90, $p = 2.47e-36$). Similarly, in the SVM predictor, **65.0%** (13/20) of FPs contain an N-terminal TM segment (odds ratio = 16.38, $p = 7.22e-09$), confirming that hydrophobic N-terminal TM anchors mimic the SP h-region.
2. **False Negative (FN) Analysis:** Sequence logo comparison of von Heijne benchmark FNs ($n = 51$) against training positives revealed a marked depletion of Alanine at positions $-1$ and $-3$ and reduced Leucine conservation in the $[-13, -6]$ h-region. For the SVM model ($n = 34$ FNs), comparative analysis across Training Positives, Benchmark TPs, and Benchmark FNs demonstrated that missed SPs have significantly longer cleavage-site positions (median SP length = 24.5 aa vs. 22.0 aa in training positives, Mann-Whitney $p = 2.14e-03$), placing a portion of their hydrophobic core outside the fixed $K = 30$ N-terminal composition window, alongside lower peak hydropathy (`kd_max`).
