# 🧫 Project Introduction 

In this project repository I'll be analyzing RNA-sequencing data from a study titled, Amiloride, [An Old Diuretic Drug, Is a Potential Therapeutic Agent for Multiple Myeloma](https://aacrjournals.org/clincancerres/article/23/21/6602/259285/Amiloride-An-Old-Diuretic-Drug-Is-a-Potential). 

In this study, myeloma cell lines and a xenograft mouse model (i.e, a mouse with human tumor cells implanted in it) were used to evaluate the drug amilioride's cytoxicity (i.e. cell killing effects). Additionally, amilioride's mechanism of action was investigated using RNA-Seq experiments, qRT-PCR, immunoblotting, and immunofluorescence assays.

The investigators in this study found that amilioride induced apoptosis in a broad panel of multiple myeloma cell lines and in the xenograft mouse model. Additionally, they found that amilioride has a synergistic effect when combined with other drugs such as dexamethasone and melphalan. Thus, in this project repository I will analyze the RNA-sequencing data from this study, made available via the Gene Expression Omnibus (GEO), to better understand amilioride's mechanism of action and effects. 

# 🧫 Loading Data & Exploratory Data Analayis

The data from this study was made available via the Gene Expression Omnibus, under the accension number [GSE95077](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95077).  
To access this data, I used Bash's ```wget``` command, which allows you to download files from the internet using the HTTP, HTTPS, or FTP protocols. In the code block below, I'll demonstrate how to access the count matrix data from this study, which contains a measure of gene expression for every gene in each cell line sample, and how to decompress the data using Bash's ```gunzip``` command:
```
# [Bash]
wget -O GSE95077_Normalized_Count_Matrix_JJN3_Amiloride_and_CTRL.txt.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE95077&format=file&file=GSE95077%5FNormalized%5FCount%5FMatrix%5FJJN3%5FAmiloride%5Fand%5FCTRL%2Etxt%2Egz'

gunzip GSE95077_Normalized_Count_Matrix_JJN3_Amiloride_and_CTRL.txt
```
Following that, I'll show you how to import the libraries needed for this project, and how to load the count matrix data into a Panda's DataFrame:
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
data = pd.read_csv(file_path, sep='\t', index_col=0)
data.head()
```
Which, produces the following output:




Data frame where index is gene ID (ENSG), then there are six columns with RNA sequencing expression data. The first three columns are a JJ cell line that has not been exposed to a drug, then the last three columns are a JJ cell line that was exposed to a drug.
