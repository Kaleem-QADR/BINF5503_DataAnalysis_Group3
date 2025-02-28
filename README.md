# Cancer Transcriptomics Dataset Analysis

## Authors

Daniel Esguerra, Reethigha Uruthiran, Siddhi Naik, Kaleem Uddin Mohammed

## 1. Introduction

### Field Overview
The dataset originates from the field of **cancer transcriptomics**, a branch of molecular biology that studies gene expression through RNA sequencing (**RNA-seq**) in cancer cells. Transcriptomics provides insight into how genes are regulated in various conditions, such as **disease progression or drug response**.

In cancer research, **patient-derived organoids (PDOs)** have emerged as a promising model system. PDOs are **three-dimensional cultures** grown from tumor cells, maintaining the genetic and phenotypic characteristics of the patient’s tumor. In the case of **high-grade serous ovarian carcinoma (HGSOC)**, transcriptomics can reveal changes in gene expression that occur during culture and treatment, providing potential targets for **precision medicine** [Chen H., et al., 2020; Nero C., et al., 2021; Maenhoudt N., et al., 2020].

### Dataset Context
The dataset focuses on **RNA sequencing of PDOs** derived from **malignant effusions of HGSOC patients**. These effusions are **fluid accumulations** that often contain **chemotherapy-resistant multicellular spheroids**. By studying **gene expression changes at two time points (Day 0 and Day 6 of short-term culture)**, we aim to uncover **adaptive mechanisms** of the tumor cells transitioning from **in vivo environments to ex vivo conditions**, as well as potential biomarkers for **drug resistance**.

Previous studies, such as **Chen H., et al. (2020)**, have shown that **short-term organoid cultures from HGSOC malignant effusions** exhibit **upregulation of proliferation-related genes, epithelial-mesenchymal transition (EMT) markers, and KRAS signaling pathways**, providing insights into **drug resistance mechanisms** [Chen H., et al., 2020].

## 2. Research Question and Hypothesis

**Research Question:**
How does the gene expression profile change between **Day 0 and Day 6** in **HGSOC organoids**?

**Hypothesis:**
Transcriptomic profiles at **Day 6** will show **significant changes in stress response and metabolic pathways** due to culture conditions.

**Prediction:**
Genes associated with **cell survival, proliferation, and metabolism** will be **upregulated**.

## 3. Data Description

### 3.1 Dataset Overview
The dataset consists of **RNA-seq read counts** from **four HGSOC patient-derived organoids**, with samples collected at **two time points (Day 0 and Day 6)**. Each file includes **three columns**:
- **Ensembl Gene IDs**
- **Gene Names**
- **RNA-seq Read Counts**

#### **Files Provided:**
- `A778_D0.tsv`, `A778_D6.tsv`
- `A820_D0.tsv`, `A820_D6.tsv`
- `A870_D0.tsv`, `A870_D6.tsv`
- `A899_D0.tsv`, `A899_D6.tsv`

### 3.2 Data Wrangling and Exploratory Code

#### **Data Preprocessing Steps**
- **Loaded the data** using `tidyverse`’s `readr` package.
- **Checked for missing values** to ensure data completeness.
- **Normalized read counts** using the `DESeq2` package.
- **Filtered out low-expression genes, outliers, and ambiguous reads** to improve statistical power.
- **Merged samples from the same patient** for comparative analysis.

#### **Data Characteristics**
- Each file contains **58,725 rows representing genes** and three columns.
- Gene expression varies between **Day 0 and Day 6**, reflecting **biological changes during culture**.

## 4. Analysis Workflow

### **Step 1: Data Loading & Wrangling**
- Merged datasets into a structured format.
- Transformed data into a **count matrix** suitable for differential expression analysis.

### **Step 2: Differential Expression Analysis (DEA)**
- Used **DESeq2** for normalization.
- Set **Day 0** as the reference condition.
- Identified **differentially expressed genes (DEGs)** with `padj < 0.05`.

### **Step 3: Visualization**
- **Volcano Plot:** Displays significantly upregulated/downregulated genes.
- **Heatmap:** Visualizes expression patterns of the **top 30 DEGs**.
- **Summary Table:** Highlights the **top 10 DEGs** (5 upregulated, 5 downregulated).

## 5. Output Files

### **Generated Reports and Figures**
- `1_CancerTranscriptomics/reports/DESeq2_Results.csv` → Full list of DEGs.
- `1_CancerTranscriptomics/reports/Top_10_DEGs.csv` → Summary of top differentially expressed genes.
- `1_CancerTranscriptomics/reports/Volcano_Plot.png` → Visualization of significant genes.
- `1_CancerTranscriptomics/reports/Heatmap_Top_30_DEGs.png` → Heatmap of top DEGs.

## 6. Usage Instructions

1. **Ensure the dataset is placed in:** `1_CancerTranscriptomics/read_counts`
2. **Run the R script provided (`Das_Assignment_2_Code.R`).**
3. **Generated outputs will be saved in:** `1_CancerTranscriptomics/reports`

## 7. References

- Chen H, Gotimer K, De Souza C, Tepper CG, Karnezis AN, Leiserowitz GS, Chien J, Smith LH. *Short-term organoid culture for drug sensitivity testing of high-grade serous carcinoma.* Gynecol Oncol. 2020 Jun;157(3):783-792. doi: 10.1016/j.ygyno.2020.03.026.
- Nero C, Vizzielli G, Lorusso D, Cesari E, Daniele G, Loverro M, Scambia G, Sette C. *Patient-derived organoids and high-grade serous ovarian cancer: from disease modeling to personalized medicine.* J Exp Clin Cancer Res. 2021 Mar 31;40(1):116. doi: 10.1186/s13046-021-01917-7.
- Maenhoudt N, Defraye C, Boretto M, Jan Z, Heremans R, Boeckx B, Hermans F, Arijs I, Cox B, Van Nieuwenhuysen E, Vergote I, Van Rompuy AS, Lambrechts D, Timmerman D, Vankelecom H. *Developing Organoids from Ovarian Cancer as Experimental and Preclinical Models.* Stem Cell Reports. 2020 Apr 14;14(4):717-729. doi: 10.1016/j.stemcr.2020.03.004.
