# Systems Biology Project: Differential Pathway Analysis

This repository contains the scripts and documentation for a Systems Biology project focused on identifying and analyzing **differentially regulated biological pathways** between two distinct conditions (e.g., Disease vs. Control, Tumor vs. Normal).

The project utilizes gene expression data to perform statistical and topological pathway analyses, fulfilling the requirements for the final assessment (PowerPoint report submission).

---

##  Project Goal

The primary objective is to take a given gene expression **dataset** (featuring two comparative conditions) and analyze it using three main methodological approaches to extract meaningful biological insights:

1.  **Over-representation Analysis (ORA):** To identify pathways where a statistically significant number of differentially expressed genes are enriched.
2.  **Functional Class Scoring (FCS):** To evaluate pathway activity by integrating the expression changes of all pathway members (e.g., GSEA or similar methods).

---

##  Repository Structure

The project is organized into the following main directories and files:

###  Analysis Scripts

| File Name | Description | Status |
| :--- | :--- | :--- |
| `system_biology_preanalysis_lung_data.py` | Script for initial data cleaning, normalization, and preparation (e.g., Differential Gene Expression calculation). | Done |
| `ora_analysis.py` | Script implementing the **Over-representation Analysis** based on the list of differentially expressed genes. | Done |
| `deg_r_script.r` | Cleaning of the dataset to be performed before starting the analysis | Done |

---

## Getting Started

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/riccardocaccia/system_biology_project
    ```
2.  **Place Data:** Ensure your dataset (expression matrix, etc.) is placed in the `data/` directory.
3.  **Clean the dataset:** Assure that yoyr dataset and data are ready to be analized -> on Rstudio:
    ```bash
    deg_r_script.r
    ```
5.  **Run Pre-Analysis:** Execute the data preparation script.
    ```bash
    python system_biology_preanalysis_lung_data.py  
    ```
6.  **Run Pathway Analyses:** Execute the subsequent analysis scripts.
    ```bash
    exploration_plots.py
    python ora_analysis.py
    ```

---

## Authors

* [riccardocaccia]
* [GautieriGiuseppe]
