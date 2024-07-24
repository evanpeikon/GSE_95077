# create function to normalize and filter genes 
def filter_genes(data, min_cpm=1, min_samples=3):
    cpm = data.apply(lambda x: (x / x.sum()) * 1e6) 
    mask = (cpm > min_cpm).sum(axis=1) >= min_samples 
    return data[mask]
# filter genes 
data = filter_genes(data)
