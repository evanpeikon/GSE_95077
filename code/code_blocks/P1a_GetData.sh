wget -O GSE95077_Normalized_Count_Matrix_JJN3_Amiloride_and_CTRL.txt.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE95077&format=file&file=GSE95077%5FNormalized%5FCount%5FMatrix%5FJJN3%5FAmiloride%5Fand%5FCTRL%2Etxt%2Egz'
gunzip GSE95077_Normalized_Count_Matrix_JJN3_Amiloride_and_CTRL.txt

python3 #switch to python environment 
path = 'GSE95077_Normalized_Count_Matrix_JJN3_Amiloride_and_CTRL.txt'
data = pd.read_csv(path, sep='\t', index_col=0)
data.head()
