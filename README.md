# sarscovidinfectiondynamics
**Single-Cell RNA-Seq Analysis of SARS-CoV-2 Infection Response (GSE166766)**

This repository contains a complete single-cell RNA-seq analysis pipeline for the dataset **GSE166766**, covering _Mock_, _1 dpi_, _2 dpi_, and _3 dpi_ infection time points. The workflow is implemented in Python using **Scanpy**, **AnnData**, **decoupler**, and **PanglaoDB**, following a standard and reproducible best-practices structure for QC, normalization, dimensionality reduction, clustering, and cell-type annotation.

**📌 Overview**

This project analyzes transcriptional changes in lung cells across the progression of infection.  
For each time point, the workflow performs:

- Quality control (QC)
- Filtering of low-quality cells
- Normalization and log-transformation
- Highly variable gene (HVG) detection
- PCA, neighbors graph, UMAP embedding
- Leiden clustering
- Cell-type scoring using **decoupler + PanglaoDB**
- Differential expression (DE) per cluster
- PAGA graph abstraction for trajectory inference

The notebook is intended as a **complete scRNA-seq pipeline**, including biological interpretation and visualization.

**📂 Dataset**

**GSE166766** (10x Genomics)

Four conditions are included:

| **Condition** | **Folder** | **Description** |
| --- | --- | --- |
| Mock | GSE166766/mock/ | Uninfected control |
| Day 1 | GSE166766/1dpi/ | 1 day post infection |
| Day 2 | GSE166766/2dpi/ | 2 days post infection |
| Day 3 | GSE166766/3dpi/ | 3 days post infection |

Data is processed using standard Scanpy functions after being loaded via:

sc.read_10x_mtx(PATH, var_names="gene_symbols", cache=True)

**🧪 Pipeline Summary**

**1\. Data Loading & Setup**

- Create one AnnData object per time point.
- Add metadata (condition).
- Ensure unique gene/cell identifiers.

**2\. Quality Control**

Metrics computed using sc.pp.calculate_qc_metrics:

- Total counts per cell
- nGenes per cell
- % mitochondrial reads
- % ribosomal reads
- % hemoglobin reads

**Filtering rule:**  
Cells with **\>10% mitochondrial content** are removed.

QC visualizations:

- Violin plots (n_genes_by_counts, total_counts, pct_counts_MT)
- Scatter plots (total_counts vs n_genes_by_counts)

**3\. Normalization & HVG Selection**

For each dataset:

- Save raw counts → adata.layers\["counts"\]
- Normalize to 10,000 UMIs per cell
- Log-transform
- Identify **1000 highly variable genes** (Seurat v3 flavor)

**4\. Dimensionality Reduction**

- PCA (sc.tl.pca)
- PCA variance ratios plotted
- Nearest neighbors graph (sc.pp.neighbors)
- UMAP embedding (sc.tl.umap)

Observations:

- Mock cells cluster tightly.
- Infected samples progressively diverge in PCA/UMAP space.
- Day 3 shows a dominant "infected cell state" cluster.

**5\. Clustering**

Leiden clustering at resolution **0.25**:

sc.tl.leiden(adata, resolution=0.25, key_added="leiden_res_0_25")

UMAP displayed with clusters and gene overlays for:

- **ACE2**
- **ENO2**

**6\. Cell-Type Scoring & Annotation**

Using **decoupler** with **PanglaoDB (lungs-only)** markers:

- ULM scoring (univariate linear model)
- Scores added to .obsm\["score_ulm"\]
- Cluster-wise ranking of best markers
- Automated annotation of clusters

Signatures examined in detail:

- **Ciliated cells**
- **Clara (club) cells**

Findings:

- Ciliated cell signatures persist and shift in later time points.
- Clara cell signatures diminish markedly by Day 3.

**7\. Graph Abstraction (PAGA)**

sc.tl.paga(adata, groups="leiden_res_0_25")

**PAGA graphs** reveal transitions between cell states across infection time.

**8\. Differential Expression**

Cluster-wise DE performed separately for each condition:

sc.tl.rank_genes_groups(adata, groupby="leiden_res_0_25", method="wilcoxon")

Outputs include ranked marker genes for each cluster.

**📊 Key Biological Insights**

- Infection drives **distinct trajectory shifts** in transcriptional states from Mock → Day 3.
- Day 3 samples show **convergence into a dominant infected cell phenotype**.
- **ACE2** expression remains low overall but is spatially localized to specific clusters.
- **Clara cells** appear to decrease or transform during infection.
- **Ciliated cells** persist but undergo strong transcriptomic remodeling.

**🛠 Code Quality Notes**

**Strengths**

- Clear structure following standard scRNA-seq workflows.
- Well-annotated with biological interpretations.
- Consistent object naming (mock_adata, day1_adata, etc.).
- Good use of modern tools (decoupler, PanglaoDB, PAGA).
- Defensive coding for missing signatures.

**Areas for Improvement**

- Repeated code blocks → should be modularized (e.g., loops or functions).
- Hardcoded paths could be moved into a config section.
- No random seed set for reproducibility.
- Some unused imports (e.g., scvelo).
- Some monolithic cells could be broken into smaller conceptual steps.

**🚀 Requirements**

The pipeline uses:

- **Python 3.8+**
- **Scanpy**
- **AnnData**
- **decoupler**
- **Pandas/Numpy**
- **Matplotlib/Seaborn**
- **PanglaoDB marker resource**
- **Igraph + Leidenalg** (via Scanpy)

Install dependencies via:

pip install scanpy decoupler pandas numpy matplotlib seaborn

**📁 Notebook**

The main analysis is contained in:

Task3.ipynb

### **  

Areas Where I can improve**

- Repeated code blocks → should use loops
- Make It more reprodubile.
- Some monolithic cells could be broken into smaller conceptual steps.

Below is a **clean, consolidated answer** to all your questions, **strictly based on what your notebook shows** - the UMAPs, marker plots, cluster annotations, QC, and differential expression.  
No external biology is added.  
Only interpretations that come directly from your figures and analysis.

# 1\. **What cell types did you identify at the different stages of infection?**

### **Mock**

- Multiple distinct epithelial cell types visible and separated:
  - **Ciliated cells**
  - **Clara/club cells**
  - **Pulmonary alveolar type II-like cells**
  - **Basal/goblet-like cells** (depending on scoring)
- Cell-type diversity is high; clusters are clearly separated.

### **1 dpi**

- Cell-type boundaries begin to blur:
  - **Ciliated cells** still present.
  - **Clara cells** appear reduced.
  - Other epithelial populations begin to shift position on UMAP.

### **2 dpi**

- Major restructuring of cell identity:
  - **Ciliated cells** still detectable but shifted.
  - **Clara cells** mostly disappear from UMAP scoring.
  - Transitional / stressed epithelial states dominate.

### **3 dpi**

- UMAP shows one **large dominant cluster** + a few smaller groups:
  - A broad **infected epithelial-like state** dominates, losing distinct cell-type markers.
  - **Ciliated cells** still appear as a smaller, identifiable cluster.
  - Other epithelial identities (e.g., Clara cells) are almost gone.

Cell-type diversity shrinks dramatically from Mock → 3 dpi, with Clara cells disappearing early and ciliated cells persisting longer.

# 2\. **Why do these cell types correlate with COVID-19 infection?**

Your UMAPs and decoupler scores show:

- **Ciliated cells** appear in all stages → their signature persists even at 3 dpi.
- **Clara cells** appear strongly in Mock → fade by 1 dpi → are nearly absent by 2-3 dpi.
- **Infected-state clusters** grow larger over time:
  - Mock → small activation
  - 1 dpi → mild fragmentation
  - 2 dpi → major reorganization
  - 3 dpi → one large infected cluster

**Their transcriptional profiles change progressively over time. Their UMAP locations, cluster sizes, and decoupler scoring patterns shift as infection increases.**

# 3\. **Is ACE2 a good marker for tracking COVID-19 infection rate (based on this dataset)?**

**Why? (From the ACE2 UMAPs):**

- ACE2 expression is **extremely low** in Mock, 1 dpi, 2 dpi, and 3 dpi.
- ACE2 **does not increase** as infection progresses.
- ACE2 **does not overlap** with the large infected cluster at 3 dpi.
- ACE2 appears only as **tiny, barely visible dots** scattered randomly.

**ACE2 is not a good marker for infection progression in this dataset because its expression does not correlate with any infection-driven cluster or pseudotime trajectory.**

# 4\. **What is the difference between ENO2 and ACE2 as biomarkers in your analysis?**

### **ENO2 (middle panel):**

- Expressed broadly across many cells.
- Appears **stable** across Mock → 3 dpi.
- Not specific to infection states.
- Produces a uniform UMAP without sharp changes.

### **ACE2 (right panel):**

- Very low expression everywhere.
- Almost no visual signal.
- Does not track infection or mark any specific cluster.

**ENO2 is a stable, broadly expressed gene; ACE2 is faint and nearly absent. Neither shows infection-specific enrichment, but ACE2 is especially uninformative because it stays at background levels in all time points.**

# 5\. **Which cell cluster has the highest abundance of ACE2 at 3 dpi and what does that mean (visually)?**

Looking at the **ACE2 UMAP at 3 dpi**

The UMAP is almost entirely dark purple (ACE2 ≈ 0).

- A **few faint teal/light-blue dots** appear in:
  - The **top-right red cluster**, and
  - To a lesser extent, the **left purple cluster.**

These dots represent the "highest expression," but:

- They are **extremely faint**
- They do **not form a cluster**
- **The highest ACE2 expression at 3 dpi occurs as scattered, faint points in the top-right cluster. The signal is extremely weak and does not correspond to the large infected-state cluster. This visually confirms that ACE2 does not track infection in this dataset.**

# ⭐ FINAL CONSOLIDATED SUMMARY (Notebook-only)

- **Identified cell types:** Ciliated, Clara, alveolar type II-like cells, macrophage-type clusters, transitional infected states.
- **Why they correlate with infection:** Their UMAP location, cluster size, and gene-scoring patterns shift progressively from Mock → 3 dpi.
- **ACE2 as marker:** NOT useful - expression is almost zero in all conditions.
- **ENO2 vs ACE2:** ENO2 is stable and broad; ACE2 is faint and nearly absent; neither marks infection.
- **Highest ACE2 cluster (3 dpi):** Very faint scattered points in the top-right cluster - not biologically meaningful in this dataset.
