def plot_genes_retained_by_cpm(data, min_samples=3):
    """
    Plots the number of genes retained as a function of different CPM thresholds
    
    Parameters:
    - data: DF with raw count data where rows are genes and columns are samples.
    - min_samples: min number of samples in which a gene must have a CPM greater than the threshold to be retained.
    """
    # convert raw counts to CPM to normalize the data
    cpm = data.apply(lambda x: (x / x.sum()) * 1e6)

    # define a range of CPM thresholds to test, from 0 to 5 with increments of 0.1
    thresholds = np.arange(0, 5, 0.1)
    
    # initialize list to store the # of genes retained for ea/ threshold
    genes_retained = []

    # loop through ea/ threshold value to determine the # of genes retained
    for min_cpm in thresholds:
        # create mask where CPM > min_cpm in at least min_samples samples
        mask = (cpm > min_cpm).sum(axis=1) >= min_samples
        
        # count # of genes that meet the criteria and append to the list
        genes_retained.append(mask.sum())

    # plot # of genes retained as a function of CPM threshold
    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, genes_retained, marker='o', color='green')
    
    # add vertical line at CPM = 1 as a rough heuristic for comparison
    plt.axvline(x=1.0, color='red', linestyle='--', label='CPM = 1')

    # add labels
    plt.xlabel('Threshold (CPM)')
    plt.ylabel('Num Genes Retained')
    plt.legend()
    plt.show()

# call function
plot_genes_retained_by_cpm(data_frame)
