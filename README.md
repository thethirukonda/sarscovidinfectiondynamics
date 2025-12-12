**Single-Cell RNA-Seq Analysis of SARS-CoV-2 Infection Response (GSE166766)**

This repository contains a complete single-cell RNA-seq analysis workflow for the dataset GSE166766, covering four infection stages: Mock, 1 dpi, 2 dpi, and 3 dpi.  
All analyses were performed using Scanpy, AnnData, ULM scoring (decoupler), and PanglaoDB lung markers, following best practices for QC, normalization, clustering, annotation, and statistical validation.

The notebook serves as a full computational pipeline from raw 10x matrices to biological interpretation, including cell-type dynamics, gene-level responses, and marker evaluation.

**📌 Overview**

This project characterizes transcriptional changes across the progression of SARS-CoV-2 infection in lung-derived single cells.  
For each time point, the workflow performs:

- Quality Control (QC)
- Filtering of low-quality cells
- Normalization & log-transformation
- Highly variable gene (HVG) detection
- PCA, neighborhood graph, UMAP embedding
- Leiden clustering
- Cell-type annotation using decoupler + PanglaoDB
- ACE2 and ENO2 marker evaluation
- Pseudotime (DPT) analysis
- Statistical validation of clustering and marker expression

The notebook integrates visualization, quantification, and biological interpretation at each step.

**📂 Dataset: GSE166766 (10x Genomics)**

Four conditions were analyzed:

| Condition | Folder | Description |
| --- | --- | --- |
| Mock | /mock/ | Uninfected control |
| 1 dpi | /1dpi/ | 1 day post infection |
| 2 dpi | /2dpi/ | 2 days post infection |
| 3 dpi | /3dpi/ | 3 days post infection |

Data were loaded using:

sc.read_10x_mtx(PATH, var_names="gene_symbols", cache=True)

Each time point is stored in a separate AnnData object (e.g., mock_adata, day1_adata, …).

**Pipeline Summary**

**1\. Data Loading & Setup**

- One AnnData object per time point
- Metadata added (adata.obs\["condition"\])
- Ensured unique cell & gene identifiers

**2\. Quality Control (QC)**

Using sc.pp.calculate_qc_metrics, the following were computed:

- Total counts per cell
- Number of genes per cell
- Percent mitochondrial reads
- Percent ribosomal reads
- Percent hemoglobin reads

**Filtering rule:**  
Cells with **\>10% mitochondrial content** were removed.

**QC Visualizations:**

- Violin plots
- Scatter plots (counts vs. genes)

**3\. Normalization & Highly Variable Genes**

For each dataset:

- Raw counts saved: adata.layers\["counts"\]
- Normalize to 10,000 counts per cell
- Log-transform
- Select **1000 HVGs** (Seurat v3 flavor)

**4\. Dimensionality Reduction**

- PCA (variance ratio plotted)
- Neighborhood graph construction
- UMAP embedding

**Observations:**

- Mock cells form compact, stable clusters
- Infected samples progressively diverge
- 3 dpi displays a strong infection-driven shift

**5\. Clustering**

Leiden clustering (resolution 0.25):

sc.tl.leiden(adata, resolution=0.25, key_added="leiden_res_0_25")

UMAPs display:

- Cluster structure
- Marker expression overlays (ACE2, ENO2)

**6\. Cell-Type Annotation (ULM scoring + PanglaoDB)**

Using **decoupler**:

- Compute ULM scores → .obsm\["score_ulm"\]
- Identify highest-score cell type per cluster
- Annotate clusters based on lung-specific marker sets

**Key Findings:**

- Clara/club cell signatures decline sharply by 3 dpi
- Ciliated cell signatures remain but undergo remodeling
- Ionocytes and macrophages become more prominent
- Type II alveolar cells contribute early in infection

**7\. Pseudotime Analysis (DPT)**

- Root cluster selected per dataset
- DPT pseudotime computed
- Pseudotime mapped to UMAP

**Outcome:**  
Gradients reflect infection-driven transitions from Mock → 3 dpi.

**8\. Statistical Testing**

Performed for:

- ACE2 expression differences across clusters and conditions
- Cluster identity stability metrics

**Tests used:**

- Kruskal-Wallis
- Mann-Whitney U with FDR correction
- Chi-square ACE2+ enrichment
- Silhouette score
- ARI subsampling (cluster robustness)
- Bootstrap 95% CIs for ACE2 means

**Insights:**

- ACE2 is significantly elevated in infected samples
- Ionocytes show strongest ACE2 enrichment at 3 dpi
- ENO2 shows no meaningful pattern

**📊 Key Biological Insights**

- Infection induces a strong transcriptomic shift by **3 dpi**.
- Ciliated cells and alveolar macrophages show substantial remodeling.
- Clara/club cells diminish or transition during infection.
- ACE2 is **low in magnitude but clearly upregulated** across infected samples.
- ENO2 shows **no detectable biological pattern**.
- Ionocytes consistently express the **highest ACE2 levels**, especially at 3 dpi.

**Areas for Improvement**

- Repeated blocks → should be looped or moved into functions
- Hardcoded paths → should be moved into a config file
- Reproducibility can be better

**🚀 Requirements**

This workflow uses:

- Python 3.8+
- Scanpy
- AnnData
- decoupler
- Pandas / NumPy
- Matplotlib / Seaborn
- PanglaoDB marker sets
- igraph + leidenalg (via Scanpy)

Install dependencies:

pip install scanpy decoupler pandas numpy matplotlib seaborn

**📁 Notebook**

All processing is contained in:

**Task3.ipynb**

**1\. What cell types did you identify at the different stages of infection?**

Across the four datasets (Mock, 1 dpi, 2 dpi, 3 dpi), the following major cell types were consistently identified using Ulm scoring and cluster-level marker profiles:

- **Mock:** Clara cells, pulmonary alveolar type II cells, alveolar macrophages, ionocytes, and ciliated cells.
- **1 dpi:** Pulmonary alveolar type I and type II cells, ionocytes, alveolar macrophages, and airway goblet-like cells.
- **2 dpi:** Ionocytes, airway goblet-like cells, pulmonary alveolar type I and type II cells.
- **3 dpi:** Pulmonary alveolar type II cells, ionocytes, alveolar macrophages, and ciliated cells.

The composition shifts across infection stages, particularly with changes in ionocytes, macrophages, and ciliated cells.

**2\. Why do these cell types correlate with COVID-19 infection?**

- **Ionocytes** consistently show high ACE2 expression scores and become more distinct across time points. Since ACE2 is the viral entry receptor, ionocytes naturally appear infection-associated in this dataset.
- **Alveolar macrophages** show activation signatures and cluster expansion during later infection stages, which aligns with their immune-response role.
- **Pulmonary alveolar type II cells** appear early and remain present, suggesting they are sensitive to infection-driven transcriptional changes.
- **Ciliated cells** show altered distributions at 3 dpi, consistent with infection-induced epithelial remodeling.

These correlations come from expression patterns and shifts visible in cluster scoring, UMAP structure, and pseudotime trajectories.

**3\. Is ACE2 a good marker for tracking COVID-19 infection rate (based on this dataset)?**

ACE2 as a reliable infection marker.

- **Mean ACE2 expression triples from Mock to 1 dpi**, and remains elevated at 2 and 3 dpi.
- The **95% confidence intervals do not overlap** between Mock and infected groups, indicating statistically meaningful increases.
- A **Chi-square test** shows highly significant enrichment of ACE2-positive cells across clusters.
- A **Kruskal-Wallis test** and pairwise **Mann-Whitney tests** confirm large, non-random differences between clusters.
- ACE2 shows **cluster-specific enrichment**, especially in ionocytes.

Although ACE2 appears visually faint on UMAP (because values are low in magnitude), the statistical tests reveal clear infection-dependent changes.

**4\. What is the difference between ENO2 and ACE2 as biomarkers in your study?**

- **ACE2** shows detectable increases across infection conditions and strong enrichment within specific clusters. It captures biologically meaningful variation and correlates with infection progression.
- **ENO2**, in contrast, has almost no measurable expression in any cluster or condition. UMAPs show uniformly near-zero values, and no cluster-specific or condition-specific pattern is present.

Therefore, ENO2 is not informative in this dataset, while ACE2 is statistically and biologically meaningful.

**5\. Which cell cluster has the highest abundance of ACE2 expression after 3 dpi, and what does that mean biologically (visual interpretation)?**

**Ionocytes (Cluster 2)** show the highest ACE2 expression at 3 dpi.

- In the violin plots and scatter overlays, ionocytes contain the densest and highest range of ACE2-positive values.
- Other clusters (macrophages, ciliated cells, type II cells) have sparse ACE2 expression with much lower means.
- The enrichment is strongly supported by statistics: the Chi-square and Kruskal-Wallis tests both produce extremely significant p-values.

Biologically, this means:

- Ionocytes form the primary ACE2-enriched population in your model system.
- They may represent a highly susceptible subpopulation during infection progression.
- Their ACE2 signature distinguishes them from other epithelial and immune cell types at peak infection (3 dpi).
