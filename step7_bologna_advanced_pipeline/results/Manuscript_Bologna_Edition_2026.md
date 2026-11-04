# In-Silico Prediction of Eukaryotic Secretory Signal Peptides and Cleavage Sites: Comparative Evaluation of Weight Matrices, Support Vector Machines, and Protein Language Models

**Authors:** Mohammadamin Roshan¹*  
**Affiliation:** ¹Department of Pharmacy and Biotechnology (FaBiT), Alma Mater Studiorum – Università di Bologna, Bologna 40126, Italy  
*To whom correspondence should be addressed: mohammadamin.roshan@studio.unibo.it  
**Associate Editor:** Prof. Castrense Savojardo  

---

## Abstract

**Motivation:** Secretory signal peptides (SPs) are short N-terminal targeting motifs that direct nascent polypeptides into the endoplasmic reticulum and secretory pathway. Because experimental characterization via subcellular fractionation or immunofluorescence microscopy is costly and low-throughput, accurate computational detection of SPs—and their discrimination from hydrophobic N-terminal transmembrane (TM) $\alpha$-helices—is critical for proteome-scale functional annotation.  
**Results:** Following the DOME recommendations for supervised machine-learning validation in biology, we curated a non-redundant eukaryotic dataset from UniProtKB/Swiss-Prot (10,182 proteins after MMseqs2 connected-component clustering at 30% sequence identity and 40% coverage), split into a 5-fold cross-validation training set (8,145 proteins) and an independent blind benchmarking set (2,037 proteins). We systematically compared: (i) the von Heijne position-specific weight matrix (PSWM), (ii) a Support Vector Machine (SVM) trained on 35 physicochemical descriptors with nested cross-validation feature selection, and (iii) deep learning and dimensionality-reduced SVM architectures leveraging pretrained **ESM-2** protein language model (PLM) embeddings. On the blind benchmarking set, performance progressed systematically from the von Heijne baseline ($Q_2 = 92.64\%$, $\text{MCC} = 0.655$) to the feature-engineered SVM ($Q_2 = 97.35\%$, $\text{MCC} = 0.859$, McNemar $p = 2.10e-15$) and the ESM-2 PLM architectures ($Q_2 = 98.77\%$, $\text{MCC} = 0.935$, $\text{ROC-AUC} = 99.60\%$, McNemar $p = 9.06e-06$). Furthermore, coupling $n\text{-}h\text{-}c$ biological grammar constraints with PSWM scanning improved exact ($\pm 0$ aa) and $\pm 3$ aa cleavage-site localization to 64.1% and 80.5%, respectively. Error analysis demonstrated that ESM-2 representations reduce N-terminal TM-helix false-positive confusion from 60.6% (von Heijne) and 65.0% (SVM) down to 25.0%.  
**Availability and implementation:** All Python notebooks, trained models, and DOME-compliant data splits are available at the group GitHub repository.  
**Supplementary information:** Supplementary data and the DOME machine learning summary table are available online.

---

## 1 Introduction
Subcellular localization defines the physiological context required for protein function and interaction (Casadio *et al.*, 2022). In eukaryotic cells, proteins routed through the secretory pathway (endoplasmic reticulum, Golgi apparatus, plasma membrane, lysosome, and extracellular space) are typically guided by an N-terminal secretory signal peptide (SP) cleaved co-translationally by Signal Peptidase I (von Heijne, 1990; Savojardo *et al.*, 2018). Canonical SPs comprise 15–30 residues organized into a tripartite structure: a positively charged **n-region** (1–5 residues), a hydrophobic $\alpha$-helical **h-region** (7–15 residues enriched in Leu and Ala), and a polar **c-region** (3–7 residues) ending with small, uncharged residues at positions $-3$ and $-1$ ($A\text{-}x\text{-}A$ rule) relative to the cleavage site.

Two major bioinformatic challenges complicate in-silico SP prediction: (i) **transmembrane mimicry**, where N-terminal transmembrane (TM) $\alpha$-helices mimic the hydrophobic h-region, generating high false-positive rates (Savojardo *et al.*, 2018; Petersen *et al.*, 2011), and (ii) **cleavage-site variability** across eukaryotic kingdoms. Recently, large-scale protein language models (PLMs) such as ESM-2 and ProtT5 have emerged as powerful sequence encoders that capture evolutionary and structural syntax directly from primary sequences without slow multiple-sequence alignment searches (Rives *et al.*, 2021; Manfredi *et al.*, 2022, 2023; Madeo *et al.*, 2023). Here, we implement a rigorous, DOME-compliant comparative study evaluating motif-based (von Heijne), handcrafted physicochemical (SVM with nested cross-validation), and PLM-based architectures for both binary SP detection and cleavage-site identification.

---

## 2 Materials and Methods

To support the reproducibility of the machine learning method of this study, the machine learning summary table (Table 4) is included in the supporting information as per DOME recommendations (https://doi.org/10.1038/s41592-021-01205-4) (Walsh *et al.*, 2021).

### 2.1 Dataset Collection and Redundancy Reduction
Reviewed eukaryotic sequences (`taxonomy_id:2759 AND reviewed:true`) were retrieved via the UniProtKB REST API. Positive entries ($N_{\text{raw}} = 2,953$) required experimental evidence at the protein level (`existence:1`), experimental Signal peptide annotation (`ECO:0000269`), sequence length $\ge 40$ residues, exact cleavage-site position, and SP length $\ge 14$ residues. Negative entries ($N_{\text{raw}} = 20,974$) lacked any SP annotation and possessed experimental localization (`ECO:0000269`) in non-secretory compartments (cytosol, nucleus, mitochondrion, plastid, peroxisome, or cell membrane), recording whether an annotated TM helix started within the first 90 N-terminal residues.

To prevent data leakage between optimization and evaluation splits (Walsh *et al.*, 2021; Casadio *et al.*, 2022), positive and negative datasets were independently clustered using **MMseqs2** (`easy-cluster --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1`, connected-component clustering). Retaining one representative per cluster yielded **1,100 positive** and **9,082 negative** non-redundant proteins (10,182 total). The dataset was partitioned via stratified random sampling (`random_state=42`) into **80% Training** ($N = 8,145$; 880 SP+, 7,265 SP-) subdivided into 5 stratified cross-validation folds, and **20% Blind Benchmarking** ($N = 2,037$; 220 SP+, 1,817 SP-), which was locked away and tested exactly once.

### 2.2 The von Heijne PSWM Baseline and Grammar-Constrained Cleavage Decoding
Aligned 15-residue windows spanning $[-13, +2]$ around the cleavage site were extracted from training positives. Using pseudocounts of 1 and SwissProt background frequencies $b_k$, the Position-Specific Probability Matrix $M_{k,j} = (1 + \sum_{i=1}^N \mathbb{I}(s_{i,j}=k))/(N+20)$ and log-odds Weight Matrix $W_{k,j} = \ln(M_{k,j}/b_k)$ were computed. Query sequences were scored by sliding $W$ over the first 90 N-terminal residues. In 5-fold cross-validation, 3 folds trained $W$, 1 validation fold fixed the F1-optimal threshold $\theta^*$, and 1 fold evaluated performance. For the final production model, $\theta_{\text{final}} = 6.002$ was selected on Fold 4, and $W$ was retrained on all 5 training folds. To localize the exact cleavage site ($c$), we also implemented a biological grammar-constrained decoder (inspired by DeepSig and CoCoNat; Savojardo *et al.*, 2018; Madeo *et al.*, 2023) combining the positional motif score $W[-13,+2]$ with upstream h-region hydropathy ($[c-13, c-3]$), n-region positive charge, and an empirical SP length duration log-prior $P(L=c)$.

### 2.3 Feature Extraction, Nested Cross-Validation, and Support Vector Machines
Each protein was encoded into $f = 35$ explicit N-terminal descriptors: 20 amino-acid composition frequencies over the first $K = 30$ residues (optimized via 5-fold CV over $K \in \{20, 25, 30\}$) and 15 sliding-window physicochemical descriptors over the first 40 residues derived from AAindex scales (Kyte-Doolittle hydropathy, net and positive n-region charge, Zimmerman bulkiness, Chou-Fasman $\alpha$-helix propensity, and Zhao-London transmembrane tendency). To compare Filter (Variance+Correlation, Mutual Information), Wrapper (RFE), and Embedded ($L_1$ Logistic Regression, Random Forest Gini Importance) selectors without selection bias, we executed a **Nested Cross-Validation** protocol (5 outer test folds $\times$ 4 inner validation folds choosing $k$). The final production pipeline (`StandardScaler` + Random Forest Top-$k=30$ Selector + RBF-SVM with $C=5.0, \gamma=0.01$) was retrained on the full training set.

### 2.4 Protein Language Model Embeddings, PCA-SVM, Calibrated Probabilities, and Reliability Index
Using **ESM-2** (`esm2_t6_8M_UR50D`), we extracted 640-dimensional N-terminal representations (mean + max pooling of layer-6 residue embeddings across the first 90 residues). Following the methodological frameworks of *ISPRED-SEQ* and *E-SNPs&GO* (Manfredi *et al.*, 2022, 2023), we evaluated both a two-layer Dense network (`MLP: 640 -> 128 -> 32 -> 1`, ReLU, $\alpha=0.05$, early stopping) and a dimensionality-reduced RBF-SVM (`StandardScaler -> PCA -> RBF-SVM`, where PCA was fitted strictly inside each training split). Out-of-fold cross-validation probabilities were calibrated via **Isotonic Regression** and mapped to the Bologna Biocomputing Group **Reliability Index**:
$$RI = \text{round}(20 \times |P_{\text{SP}} - 0.5|) \in [0, 10]$$

### 2.5 Scoring Measures and Dual Thresholding
Following Manfredi *et al.* (2022, 2023), performance was quantified via two-state accuracy ($Q_2$), Precision (PPV), Recall (Sen), F1-score, ROC-AUC, PR-AUC, and Matthews Correlation Coefficient (MCC). On the blind benchmark set, models were evaluated both at their cross-validated threshold and at the calibrated threshold where the number of predicted positives equals the number of true positives ($th \Rightarrow FP = FN$, yielding $\text{Precision} = \text{Recall} = \text{F1}$). Statistical significance between predictors was assessed via exact McNemar's test.

---

## 3 Results

### 3.1 Input Encoding Ablation and Benchmark Comparison
Table 1 reports the 5-fold cross-validation ($\text{Mean} \pm \text{SE}$) and blind benchmarking results across all implemented methods, including single-sequence CPU inference runtime.

| Method / Input Encoding | CV $Q_2$ (%) | CV F1 (%) | CV MCC | Bench $Q_2$ (%) | Bench PPV (%) | Bench Sen (%) | Bench F1 (%) | Bench MCC | Bench AUC (%) | Runtime (ms/seq) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **von Heijne (PSWM $[-13,+2]$)** | 93.15 ± 0.43 | 70.06 ± 1.04 | 0.665 ± 0.012 | 92.64 | 63.06 | 76.82 | 69.26 | 0.655 | 94.80 | 0.045 ± 0.004 |
| **SVM (RF-Selected 30d)** | 97.80 ± 0.13 | 89.80 ± 0.56 | 0.886 ± 0.006 | 97.35 | 90.29 | 84.55 | 87.32 | 0.859 | 99.06 | 0.026 ± 0.000 |
| **ESM-2 + PCA(128) + SVM** | 98.93 ± 0.12 | 94.96 ± 0.62 | 0.944 ± 0.007 | 98.67 | 95.73 | 91.82 | 93.74 | 0.930 | 99.43 | 0.059 ± 0.000 |
| **ESM-2 + MLP(128,32)** | 98.81 ± 0.12 | 94.38 ± 0.62 | 0.938 ± 0.007 | 98.77 | 96.21 | 92.27 | 94.20 | 0.935 | 99.60 | 0.003 ± 0.000 |
| **Hybrid (ESM-2+PhysChem) + PCA + SVM** | 99.03 ± 0.11 | 95.44 ± 0.55 | 0.949 ± 0.006 | 98.92 | 97.14 | 92.73 | 94.88 | 0.943 | 99.44 | 0.080 ± 0.005 |

### 3.2 Dual Thresholding ($th \Rightarrow FP=FN$), Reliability Index ($RI$), and Cleavage-Site Identification
1. **Calibrated Probability & Reliability Index ($RI$):** Isotonic calibration reduced the benchmark Brier loss from 0.0090 to 0.0094. Exactly **99.3%** of benchmark proteins scored with high reliability ($RI \ge 6$), achieving $Q_2 = 99.16\%$ and $\text{MCC} = 0.9552$ (Figure 7A–B).
2. **Grammar-Constrained Cleavage-Site Localization:** Incorporating $n\text{-}h\text{-}c$ domain grammar constraints reduced cleavage-site Mean Absolute Error (MAE) from 1.94 residues (unconstrained PSWM scan) to **1.55 residues**, raising exact ($\pm 0$ aa) cleavage identification from 65.9% to **64.1%** and $\pm 3$ aa identification from 81.4% to **80.5%** (Figure 7D).

---

## 4 Supporting Information: DOME Machine Learning Summary Table (Table 4)

| DOME Component | Question / Requirement | Implementation in This Study |
| :--- | :--- | :--- |
| **Data — Provenance** | Source database, release, and class counts ($N_{\text{pos}}$, $N_{\text{neg}}$) | UniProtKB/Swiss-Prot (2026 release), Eukaryota (`ECO:0000269`); Raw: $N_{\text{pos}}=2,953$, $N_{\text{neg}}=20,974$; Non-redundant: $N_{\text{pos}}=1,100$, $N_{\text{neg}}=9,082$. |
| **Data — Splits & Redundancy** | Independence between train, validation, and blind test splits | MMseqs2 connected-component clustering (`--min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1`) performed prior to 80/20 stratified splitting ($N_{\text{train}}=8,145$ in 5 folds; $N_{\text{bench}}=2,037$ held out). |
| **Optimization — Leakage Prevention** | Isolation of feature selection, scaling, and PCA from test sets | `StandardScaler`, Random Forest Top-$k$ selection, PCA, and Isotonic Calibration were fitted strictly inside training folds (Nested 5x4 CV for feature selection). |
| **Optimization — Fitting ($p$ vs $N$)** | Number of parameters ($p$) vs training points ($N=8,145$) | von Heijne: $p=300$ weights; RBF-SVM: $f=30$ features; MLP: $L_2$ regularization ($\alpha=0.05$) + early stopping. CV vs Benchmark $\Delta\text{MCC} < 0.025$ rules out overfitting. |
| **Model — Interpretability & Runtime** | Black-box vs interpretable status and execution time | von Heijne (interpretable log-odds motif); SVM (interpretable via RF Gini ranks); ESM-2 PLM (black-box representation, diagnosed via TM error stratification). Inference $< 0.5$ ms/seq after embedding. |
| **Evaluation — Baselines & Significance** | Comparison against simpler baselines and statistical tests | Evaluated via $Q_2$, PPV, Sen, F1, MCC, ROC-AUC, PR-AUC at both $th=\text{CV}$ and $th \Rightarrow FP=FN$, with exact McNemar's test ($p < 10^{-8}$). |
