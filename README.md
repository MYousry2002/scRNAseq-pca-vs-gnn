# scRNAseq
This project attempts to explore and compare the applications of different data science approaches in preprocessing and processing single-cell RNA sequencing datasets.

## Datasets
The main dataset is single-cell RNA sequencing (scRNA-seq) data for human pancreatic islets obtained from the National Institutes of Health's Gene Expression Omnibus (GEO). The dataset has a GEO accession ID of GSE198623. This dataset consists of high-quality single-cell RNA-seq (scRNA-seq) data for the pancreatic islets from healthy adult human donors and is collected and deposited by Tritschler et al. (2022). This dataset is a collection of 5 samples obtained from 5 different donors, named 22, 24, 61, 63, and 74. Single-cell RNA sequencing (scRNA-seq) comprises thousands of individual cells' produced RNA molecules, thus providing the gene expression profiles of these cells across thousands of genes. Since gene expression profile defines cell identity, scRNA-seq is considered the gold standard for defining cell states and phenotypes.

## Motivation
#### The Significance of scRNA-seq in Biological Research
Single-cell RNA sequencing (scRNA-seq) stands at the forefront of modern biological research, providing unprecedented insights into the complex and dynamic nature of cellular processes. The ability to examine the gene expression profiles of individual cells offers a granular view of cellular diversity, function, and interaction, especially in heterogeneous tissues like the human pancreatic islets. This level of detail is crucial for understanding biological systems in health and disease, leading to potential breakthroughs in areas such as developmental biology, immunology, and personalized medicine.

#### Challenges and Opportunities in scRNA-seq Data Analysis
While scRNA-seq is a powerful tool, it also presents significant data analysis challenges. The datasets generated are large, complex, and often noisy, with issues such as batch effects, dropout events, and technical variations. Efficient and accurate preprocessing and processing of scRNA-seq data are vital to ensure the reliability of subsequent biological interpretations. This project is motivated by the need to explore, compare, and potentially improve various data science approaches to address these challenges.

#### Exploring Diverse Approaches
The project focuses on a range of methodologies, each with its unique strengths and potential applications. By comparing and contrasting different data science approaches, we aim to provide a roadmap for researchers dealing with scRNA-seq data, facilitating more accurate and insightful analyses.

## Methods and Models
### Preprocessing Analysis
1. **Alignment and Transcript Quantification**: STARsolo Software
2. **Dimensionality Reduction**: Principal Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP)
3. **Batch Correction**: Systematic correction with ComBat method, aims to reduce the effect of technical differences among samples

### Processing Analysis
1. **Modularity-Optimization Clustering Methods**: Louvain and Leiden Algorithms (Seurat package in R and Scanpy in Python).
2. **Graph Neural Networks Clustering Methods**: scGNN in R
3. **Gene Markers and Cell Typing**: Analysis and Visualization
4. **Active Pathways and Gene Modules Analysis**: Graph Neural Network approaches, scapGNN package in R.

## References
Han, X., Wang, B., Situ, C., Qi, Y., Zhu, H., Li, Y., & Guo, X. (2023). scapGNN: A graph neural network-based framework for active pathway and gene module inference from single-cell multi-omics data. PLoS biology, 21(11), e3002369. https://doi.org/10.1371/journal.pbio.3002369

Scanpy. (2017). Tutorials. Scanpy 0.1.0.dev documentation. https://scanpy.readthedocs.io/en/stable/tutorials.html#integrating-datasets

Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: Guaranteeing well-connected communities. Scientific Reports, 9(1). https://doi.org/10.1038/s41598-019-41695-z

Tritschler, S., Thomas, M., Böttcher, A., Ludwig, B., Schmid, J., Schubert, U., Kemter, E., Wolf, E., Lickert, H., & Theis, F. J. (2022). A transcriptional cross species map of pancreatic islet cells. Molecular Metabolism, 66, 101595. https://doi.org/10.1016/j.molmet.2022.101595

Wang, J., Ma, A., Chang, Y. et al. scGNN is a novel graph neural network framework for single-cell RNA-Seq analyses. Nat Commun 12, 1882 (2021). https://doi.org/10.1038/s41467-021-22197-x

Wikimedia Foundation. (2023, September 10). Louvain method. Wikipedia. https://en.wikipedia.org/wiki/Louvain_method

Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). Combat-seq: Batch effect adjustment for RNA-seq count data. NAR Genomics and Bioinformatics, 2(3). https://doi.org/10.1093/nargab/lqaa078
