# check for missing values and get summary statistics 
print(data.isnull().sum())
print(data.describe())

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

# perform hierarchical clustering and create dendrogram
h_clustering = linkage(log_counts.T, 'ward')
plt.figure(figsize=(8, 6))
dendrogram(h_clustering, labels=data.columns)
plt.xticks(rotation=90)
plt.ylabel('Distance')
plt.show()
