# 🧫 Exploring Amilioride's Mechanism of Action and Effects

In this project repository I'll be analyzing RNA-sequencing data from a study titled, Amiloride, [An Old Diuretic Drug, Is a Potential Therapeutic Agent for Multiple Myeloma](https://aacrjournals.org/clincancerres/article/23/21/6602/259285/Amiloride-An-Old-Diuretic-Drug-Is-a-Potential). 

In this study, myeloma cell lines and a xenograft mouse model (i.e., a mouse with human tumor cells implanted in it) were used to evaluate the drug amiloride's toxicity (i.e., cell-killing effects). Additionally, amilioride's mechanism of action was investigated using RNA-Seq experiments, qRT-PCR, immunoblotting, and immunofluorescence assays.

The investigators in this study found amiloride-induced apoptosis in a broad panel of multiple myeloma cell lines and in the xenograft mouse model. Additionally, they found that amiloride has a synergistic effect when combined with other drugs, such as dexamethasone and melphalan.

In this project repository, I will analyze the RNA-sequencing data from this study, made available via the Gene Expression Omnibus (GEO), to better understand amiloride's mechanism of action and effects, with special attention dedicated to analyzing differentially expressed genes and their physiological functions. If you're unfamiliar with differential expression analysis and how it works, I recommend reading the sub-section below titled "Scientific Background" before proceeding. Otherwise, you can skip the section titled "Loading Data & Exploratory Data Analysis."

### Scientific Background
**What Is Differential Gene Expression?**

Gene expression is the process by which information encoded in a gene is used to create a functional gene product, typically a protein or RNA molecule. The first step of gene expression is transcription, which is when the DNA sequence coding for a gene is transcribed into messenger RNA (mRNA) by the enzyme RNA polymerase, as demonstrated in the figure below ([image source](https://www.technologynetworks.com/genomics/articles/rna-polymerase-function-and-definition-346823)). 

<img src="images/DEA1.jpg" alt="Description" width="600" height="300">

The mRNA is then translated into a functional protein with the help of ribosomes and transfer RNA (tRNA), as shown in the image below ([image source](https://en.wikipedia.org/wiki/Transfer_RNA)). Gene expression is a tightly regulated process that allows cells to respond to environmental cues, perform specific functions, and adapt to changing conditions. Not all genes are expressed at all times, and the expression level can vary between different cell types, tissues, and developmental stages. In other words, genes can be differentially expressed, which refers to the ability of genes to be "turned on" or "turned off" in response to specific factors or changes in the cellular environment.

<img src="images/DEA2.jpg" alt="Description" width="600" height="300">

Differential expression analysis is a fundamental technique in bioinformatics used to identify differentially expressed genes between two or more biological conditions, such as healthy and diseased tissues, or before and after a treatment. To perform differential expression analysis, we must perform statistical analysis to discover quantitative changes in expression levels between experimental groups. For example, we use statistical testing to decide whether, for a given gene, an observed difference in read counts is significant, that is, whether it is greater than what would be expected just due to natural random variation.

**How Do We Know Genes Are Differentially Expressed?**

Differential expression analysis aims to understand how gene expression levels change under different conditions, providing insights into the molecular mechanisms underlying biological processes. Before performing differential expression analysis, we need to quantify gene expression levels. This can be done with RNA sequencing (RNA-seq) and microarrays, which are used to identify genes being actively transcribed. RNA-Seq, a cutting-edge technique, reads the entire transcriptome, providing a nuanced and dynamic understanding of gene activity. In contrast, microarrays, akin to genetic snapshots, capture a snapshot of gene expression levels at a given moment. 

After obtaining gene expression data, we can assess differential expression by comparing the expression of genes under different experimental conditions. For example, researchers might compare gene expression in healthy and diseased tissues or in the presence and absence of a specific treatment, or they may investigate expression under different environmental conditions.

A gene is considered upregulated if its expression increases by a statistically significant degree in a particular condition and downregulated if its expression decreases by a statistically significant degree. In the image above, genes that are upregulated in the control group as compared to the zinc-exposed group are indicated in orange, whereas genes that are down-regulated in the control group as compared to the zinc-exposed group are indicated in blue. The genes indicated in brown weren't differentially expressed between the two groups. In the next section, we'll explore how we can determine this ourselves, using the raw data from a study titled, Amiloride, An Old Diuretic Drug, Is a Potential Therapeutic Agent for Multiple Myeloma. 

## 🧫 Loading Data & Exploratory Data Analysis

The data from this study were made available via the Gene Expression Omnibus under the accession number [GSE95077](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95077).  To access this data, I used Bash's ```wget``` command, which allows you to download files from the internet using the HTTP, HTTPS, or FTP protocols. 

In the code block below, I'll demonstrate how to access the count matrix data from this study, which contains a measure of gene expression for every gene in each cell line sample, and how to decompress the data using Bash's ```gunzip``` command:
```
# [Bash]
wget -O GSE95077_Normalized_Count_Matrix_JJN3_Amiloride_and_CTRL.txt.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE95077&format=file&file=GSE95077%5FNormalized%5FCount%5FMatrix%5FJJN3%5FAmiloride%5Fand%5FCTRL%2Etxt%2Egz'

gunzip GSE95077_Normalized_Count_Matrix_JJN3_Amiloride_and_CTRL.txt
```
Following that, I'll show you how to import the libraries needed for this project and how to load the count matrix data into a Pandas DataFrame:
```
# [Python]
# import libraries
import matplotlib.pyplot as plt
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
import gseapy as gp
import numpy as np
import seaborn as sns
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

# define the file path, load the data into a DataFrame, and view the first 5 rows
path = 'GSE95077_Normalized_Count_Matrix_JJN3_Amiloride_and_CTRL.txt'
data = pd.read_csv(path, sep='\t', index_col=0)
data.head()
```
Which, produces the following output:

<img src="images/data_head.png" alt="Description" width="1000" height="175">

The data frame above is indexed by gene ID (ENSG), and then there are six columns of RNA-sequencing expression data (i.e., counts). The first three columns contain expression counts for the control group, and the last three columns contain expression counts for the experimental group (i.e., the amiloride treatment group). Note that both groups use the JJN3 cell line, which was established from the bone marrow of a 57-year-old woman with plasma cell leukemia.

Next, we'll want to perform preliminary data analysis to understand the distribution and variability of RNA sequencing counts across the samples and check for missing data before performing any downstream analysis. First, let's check out the sample quality:
```
# [Python]
# check for missing values and get summary statistics 
print(data.isnull().sum())
print(data.describe())
```
Which, produces the following output:

<img src="images/summary_stats.png" alt="Description" width="400" height="450">

Notably, the dataset has no null (missing) values, and all samples' total counts (i.e., gene IDs) are the same. Next, we'll explore the distribution and variability in the dataset, as demonstrated in the code block below:
```
# [Python]
# calcualte total counts per sample and log transform counts
total_counts = data.sum(axis=0)
log_counts = data.apply(lambda x: np.log2(x + 1)) #+1 ensures zero counts are transformed to log2(1) = 0 instead of being undefined

# create subplots
fig, axes = plt.subplots(1, 2, figsize=(18, 6))

# subplot 1: total counts per sample
axes[0].bar(data.columns, total_counts, color='skyblue')
axes[0].set_ylabel('Total Counts')
axes[0].set_title('Total Counts per Sample')
axes[0].tick_params(axis='x', rotation=85)

# subplot 2: log transformed counts per sample
log_counts.boxplot(ax=axes[1])
axes[1].set_ylabel('Log2(Counts + 1)') 
axes[1].set_title('Log Transformed Counts per Sample')
axes[1].tick_params(axis='x', rotation=85)

plt.tight_layout()
plt.show()
```
Which produces the following output:

<img src="images/EDA.png" alt="Description" width="1000" height="450">

The chart on the left, titled "Total Counts," helps visualize the overall sequencing depth across the samples. Ideally, the bars, representing the total counts, should be of similar height, indicating that sequencing depth is consistent across samples, which is the case with this data set. 

Conversely, the rightmost chart is used to assess the variability and distribution of gene expression counts across samples. In this case, we can see that the boxes, representing the interquartile range, and the whiskers, representing variability outside the upper and lower quartiles, are of similar sizes across the samples, indicating a consistent data distribution. 

Now, before moving on to quality control and filtering, we'll use one last visualization to explore the similarity, or dissimilarity, between our six samples:
```
# [Python]
# perform hierarchical clustering and create dendrogram
h_clustering = linkage(log_counts.T, 'ward')
plt.figure(figsize=(8, 6))
dendrogram(h_clustering, labels=data.columns)
plt.xticks(rotation=90)
plt.ylabel('Distance')
plt.show()
```
Which produces the following output:

<img src="images/denodrogram1.png" alt="Description" width="400" height="450">

The image above shows the results of hierarchical clustering, which can be visualized via a dendrogram. When viewing a dendrogram, special attention should be paid to the cluster groupings and branches. Samples clustered together are more similar to each other, and the length of the branches (vertical lines) connecting clusters represents the distance or dissimilarity between clusters. 

The chart above shows that our three control samples are clustered on the left, whereas our three experimental (i.e.,m amiloride-exposed) samples are clustered together on the right. This is a good sign, suggesting that the control and experimental groups are distinct and that there is biological variation between the two groups of samples. Thus, we can feel confident that our downstream differential expression analyses will provide meaningful results. 

## 🧫 Quality Control, Filtering, and Normalization

The next step in our analysis is to filter out genes with low expression levels across all samples, which can introduce noise in the data. By filtering these out, you can make your results more reliable and improve your statistical power, making detecting real biological differences between conditions easier. Additionally, filtering out genes with low expression counts decreases computational load by reducing the number of genes in your dataset, making future downstream analyses faster. In the code block below, I'll show you how perform basic filtering and normalization:
```
# [Python]
# create function to normalize and filter genes 
def filter_genes(data, min_cpm=1, min_samples=3):
    cpm = data.apply(lambda x: (x / x.sum()) * 1e6) 
    mask = (cpm > min_cpm).sum(axis=1) >= min_samples 
    return data[mask]
# filter genes 
data = filter_genes(data)
```
Now, let's break this code down step by step to see how it works:

- In the code block above, the ```filter_genes``` function normalizes the gene expression data to CPM (counts per million), then filters out genes that do not meet the criteria for minimum expression across a sufficient number of samples, thus achieving normalization and filtering in one swoop.
- CPM is a common normalization method used to account for differences in library sizes across samples. It works by scaling the raw counts to a common unit, allowing comparisons across samples. In the code block above, ```cpm = data.apply(lambda x: (x / x.sum()) * 1e6) ``` uses a lambda function. For each row (i.e., gene), the lambda function divides each value (i.e., count) by the sum of counts for that gene across all samples, then multiplies the result by 1*10^6, converting raw counts into CPM. 
- Then, ```mask = (cpm > min_cpm).sum(axis=1) >= min_samples``` identifies genes that have sufficient expression across a minimum number of samples, which helps retain genes that are sufficiently expressed in a significant number of samples. First, ```(cpm > min_cpm)```creates a boolean DataFrame where True indicates that the CPM for a given gene is greater than ```(min_cpm```. Then, ```.sum(axis=1)``` sums these boolean values across samples for each gene, resulting in the count of samplers where the gene's expression exceeds ```min_cpm```. Finally, ```>=min_samples``` checks if this count is greater than or equal to ```min_samples```, creating a boolean series where True indicates the gene meeting the filtering criteria. 
- Finally, ```return data[mask``` returns the subset of the original data that meets the filtering criteria. After calling this function with the code ```data = filter_genes(data)```, the result is that the data set is filtered from 23,044 genes down to 13,335 genes, significantly reducing the noise in the dataset. 

## 🧫 Differential Expression Analysis
Now that we've loaded our data and performed quality control and normalization, we can perform differential expression analysis. In this case, I used a pairwise analysis, which involves comparing gene expression levels between individual pairs of control and experimental samples. For example, I compared control sample 1 to experimental sample 1, control sample 2 to experimental sample 2, etc. 

Pairwise analyses are useful when working with small sample sizes, as we currently are. Additionally, pairwise comparison can be more precise because it compares matched samples, reducing variability caused by biological differences between samples and batch effects. In the code block below, I'll demonstrate how to perform a pairwise analysis, multiple testing corrections, and how to identify differentially expressed genes:
```
# [Python]
results = []
for gene in data.index:
    control = data.loc[gene, ['JJ_AMIL_141050_INTER-Str_counts', 'JJ_AMIL_141056_INTER-Str_counts', 'JJ_AMIL_141062_INTER-Str_counts']]
    treated = data.loc[gene, ['JJ_CTRL_141048_INTER-Str_counts', 'JJ_CTRL_141054_INTER-Str_counts', 'JJ_CTRL_141060_INTER-Str_counts']]
    mean_control = np.mean(control)
    mean_treated = np.mean(treated)
    log2fc = np.log2((mean_treated + 1) / (mean_control + 1))  # Adding 1 to avoid log of 0
    t_stat, p_val = ttest_ind(control, treated)
    results.append({'gene': gene, 'log2fc': log2fc, 't_stat': t_stat, 'p_val': p_val})

results_df = pd.DataFrame(results)

results_df['p_adj'] = multipletests(results_df['p_val'], method='fdr_bh')[1]

results_df['abs_log2fc'] = results_df['log2fc'].abs()

deg = results_df[(results_df['p_adj'] < 0.01) & (results_df['abs_log2fc'] > 1)]
```
The code above performs a differential expression analysis on gene expression data, and the final output, ```deg```, is a DataFrame containing the genes that are significantly differentially expressed between the control and treated samples. Now, let's break this code down step by step to see how it works:

- First, ```results = []``` initializes an empty list to store the results of the t-tests and fold change calculations for each gene.
- Then, ```for gene in data.index:``` initializes a for loop to iterate over each gene in the dataset (recall, ```data.index:``` contains the gene identifiers for this dataset).
    - 










## 🧫 Visualizing and Understanding the Results



## 🧫 Conclusion: What is Amilioride's Mechanism of Action? Is It Effective?

- bring in mouse data for effectivenss 





