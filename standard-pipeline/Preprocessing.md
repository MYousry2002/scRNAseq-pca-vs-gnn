# scRNA Preprocessing Guide

## Batch Correction: ComBat

From the UMAP plot above, we can clearly notice batch effects between our five samples in the data as points from different samples form visually separable clusters. Batch effects are systematic differences in data that arise not from biological variations but from technical or experimental differences, such as different sample preparations, sequencing runs, or lab conditions. If there were no batch effects, the points from the different five samples would ideally mix together well and not form visually seperable clusters. That is because of the assumption that different samples come from the same underlying biological process.

We need to solve this problem. Batch correction in scRNA-seq is a critical step for ensuring that subsequent analyses reflect true biological differences rather than technical artifacts. However, it is important to avoid overcorrection as this might remove genuine biological variability. Thus, we should carefully select an appropriate batch correction technique.

A common statistical method in the field is Combating Batch Effects (ComBat), which uses an empirical Bayes framework to model and adjust for batch effects in gene expression data. This method leverages the strengths of both frequentist (data-driven estimation) and Bayesian (prior and posterior updating) methods, making it particularly effective for high-dimensional genomic data where direct parameter estimation might be challenging (Zhang et al., 2020). ComBat works as follows:

1- Model Formulation: It formulates a model that includes parameters for both the biological covariates and the batch effects, as follows: 
$$ Y_{ij} = μ + α_j + β_i + ε_{iJ} $$  
Where,
- $Y_{ij}$ is the logged expression level of gene  $j$ in sample $i$
- $μ$ is the overall mean expression of the gene across all samples.
- $α_j$ is the biological effect of the j-th gene.
- $β_i$ is the batch effect for the i-th sample.
- $ϵ_{ij}$ is the residual error, assumed to be normally distributed.

2- Empirical Bayes Estimation: In a fully Bayesian approach, priors for these parameters are specified based on prior knowledge. Empirical Bayes simplifies this by estimating priors directly from the observed data. ComBat uses the Empirical Bayes approach as it involves estimating the prior distributions of parameters based on the observed data.

- Batch effects: It assumes that batch effects $β_j$ come from a shared normal distribution with a mean of zero $N(0, σ_β^2)$. The variance $σ_β^2$ is estimated from the data.

- Gene effects: Similarly, gene effects $α_i$ are also modeled with a shared normal distribution, with its own variance $σ_α^2$  estimated from the data.

The concept of "borrowing strength" comes into play here. ComBat assumes that the batch effects for all genes come from a shared distribution. This means that the information from all genes is pooled to estimate the distribution of batch effects. That allows for individual differences in batch effects for each gene but informs the estimation for each batch effect from the collective behavior of all batches. The "borrowing strength" way is especially useful when some genes are not highly expressed or when sample sizes are small, especially when dealing with high-dimensional data, like gene expression, where the number of features (genes) far exceeds the number of samples (Zhang et al., 2020). 

3- Updating with Bayes' Theorem: ComBat then uses Bayes' theorem to update these priors into posterior estimates for each gene’s batch effects.

$$ P(β_i | Y) ∝ P(Y | β_i) * P(β_i)$$ 

Where, $P(β_i ∣ Y)$ is the posterior probability of the batch effect given the observed data, $P(Y ∣ β_i)$ is the likelihood of the observed data given the batch effect, and $P(β_i)$ is the prior probability of the batch effect (estimated empirically).

The posterior estimates for each gene's batch effects are "shrunk" towards the mean of the shared distribution. This shrinkage is more pronounced for genes with less clear batch effects, preventing over-adjustment based on noisy or limited data. For genes with strong and clear batch effects, the posterior estimates will be notably different from the priors, whereas for genes with ambiguous batch effects, the posteriors will align more closely with the empirical priors. 


5- Adjusting the Data: Finally, the original gene expression data is adjusted by removing the estimated batch effects, resulting in batch-corrected data.
$$ Y_{ij}^{adjusted} = Y_{ij} - \hat β_i $$

This adjustment aims to retain biological variability while reducing technical variability, enhancing the validity of downstream analyses. ComBat's empirical Bayes approach, with its balance of data-driven estimation and Bayesian updating, is adept at mitigating overfitting and ensuring that adjustments are reflective of broader data trends.