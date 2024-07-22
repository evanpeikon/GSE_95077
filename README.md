# 🧫 Exploring Amilioride's Mechanism of Action and Effects

In this project repository I'll be analyzing RNA-sequencing data from a study titled, Amiloride, [An Old Diuretic Drug, Is a Potential Therapeutic Agent for Multiple Myeloma](https://aacrjournals.org/clincancerres/article/23/21/6602/259285/Amiloride-An-Old-Diuretic-Drug-Is-a-Potential). 

In this study, myeloma cell lines and a xenograft mouse model (i.e., a mouse with human tumor cells implanted in it) were used to evaluate the drug amiloride's toxicity (i.e., cell-killing effects). Additionally, amilioride's mechanism of action was investigated using RNA-Seq experiments, qRT-PCR, immunoblotting, and immunofluorescence assays.

The investigators in this study found amiloride-induced apoptosis in a broad panel of multiple myeloma cell lines and in the xenograft mouse model. Additionally, they found that amiloride has a synergistic effect when combined with other drugs, such as dexamethasone and melphalan. Thus, in this project repository, I will analyze the RNA-sequencing data from this study, made available via the Gene Expression Omnibus (GEO), to better understand amiloride's mechanism of action and effects.

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

The next step in our analysis is to filter out genes with low expression levels across all samples, which can introduce noise in the data. By filtering these out, you can make your results more reliable and improve your statistical power, making detecting real biological differences between conditions easier. Additionally, filtering out genes with low expression counts decreases computational load by reducing the number of genes in your dataset, making future downstream analyses faster. In the code block below, I'll show you how perform basic filtering and normalization, then I'll break the code down step by step, explaining how it works:
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
In the code block above, the ```filter_genes``` function normalizes the gene expression data to CPM (counts per million), then filters out genes that do not meet the criteria for minimum expression across a sufficient number of samples, thus achieving normalization and filtering in one swoop.

CPM is a common normalization method used to account for differences in library sizes across samples. It works by scaling the raw counts to a common unit, allowing comparisons across samples. In the code block above, ```cpm = data.apply(lambda x: (x / x.sum()) * 1e6) ``` uses a lambda function. For each row (i.e., gene), the lambda function divides each value (i.e., count) by the sum of counts for that gene across all samples, then multiplies the result by 1*10^6, converting raw counts into CPM. 

Then, ```mask = (cpm > min_cpm).sum(axis=1) >= min_samples``` identifies genes that have sufficient expression across a minimum number of samples, which helps retain genes that are sufficiently expressed in a significant number of samples. First, ```(cpm > min_cpm)```creates a boolean DataFrame where True indicates that the CPM for a given gene is greater than ```(min_cpm```. Then, ```.sum(axis=1)``` sums these boolean values across samples for each gene, resulting in the count of samplers where the gene's expression exceeds ```min_cpm```. Finally, ```>=min_samples``` checks if this count is greater than or equal to ```min_samples```, creating a boolean series where True indicates the gene meeting the filtering criteria. 

Finally, ```return data[mask``` returns the subset of the original data that meets the filtering criteria. After calling this function with the code ```data = filter_genes(data)```, the result is that the data set is filtered from 23,044 genes down to 13,335 genes, significantly reducing the noise in the dataset. 

## 🧫 Differential Expression Analysis



## 🧫 Visualizing and Understanding the Results



## 🧫 Conclusion: What is Amilioride's Mechanism of Action? Is It Effective?

- bring in mouse data for effectivenss 





