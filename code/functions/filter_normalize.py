def filter_normalize(data, min_cpm=0.75, min_samples=3):
    """
    Filter and normalize gene expression data.

    Parameters:
    - data: raw gene expression data where rows represent genes and columns represent samples.
    - min_cpm: min counts per million (CPM) threshold for filtering genes.
    - min_sample: min number of samples with CPM above min_cpm to keep a gene.

    Returns:
    - pd.DataFrame: filtered and normalized gene expression data.
    """
    
    # convert raw counts to CPM
    cpm = data.apply(lambda x: (x / x.sum()) * 1e6)  
    
    # filter genes based on CPM thresholds
    mask = (cpm > min_cpm).sum(axis=1) >= min_samples # keep genes with CPM > min_cpm in at least min_samples
    filtered_data = data[mask]  # apply filter 

    # Compute geometric mean of non-zero values for each gene
    geometric_means = filtered_data.apply(lambda row: np.exp(np.log(row[row > 0]).mean()), axis=1)
    
    # calculate size factors by dividing each gene expression by its geometric mean
    size_factors = filtered_data.div(geometric_means, axis=0).median(axis=0)
    
    # normalize data by dividing each gene expression by the size factors
    normalized_data = filtered_data.div(size_factors, axis=1)

    # return normalized data as DF w/ the same index and columns as the filtered data
    return pd.DataFrame(normalized_data, index=filtered_data.index, columns=filtered_data.columns)

# overwrite DF with filtered and normalized data
data_frame = filter_normalize(data_frame)
