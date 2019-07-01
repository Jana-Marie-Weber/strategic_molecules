# Isolation forest 

The script paper_full_routine.py performs the classification of the network with different contamination rates using the isolation forest classifier. The steps are:
1. We load the network with the computed features.
2. We iterate once through the network and save all feature values and the nodes identification to python lists. 
3. We use a min/max normalisation from sklearn to normalise the values. 
4. To prove normalisation we print the new min and max values of all features to a separate file. 
5. We perform PCA on the data. First considering all 5 principal components and looking at the correlations of the features with them and second reducing the data to two dimensions. 
6. On the reduced and normalised data we perform an isolation forest outlier detection. A set of contamination rates was decided on previously, based on rough estimates about the number of strategic mlecules. All samples are used to build the decision trees. For further parameterisation please see the code.
7. The prediction of out or inlier are saved in csv format for plotting. 
8. Rough plots are generated in python for fast evaluation of results.
9. All hubs found in that specific rate are saved to a dictionary and printed at the end. 

All results are saved.

One Warning I came across during the indexing. The version now is running without any errors, but for future application please check the package versions and make sure the indexing runs correctly. 

/home/jmw254/.local/lib/python2.7/site-packages/scipy/stats/stats.py:1713: FutureWarning: Using a non-tuple sequence for multidimensional indexing is deprecated; use `arr[tuple(seq)]` instead of `arr[seq]`. In the future this will be interpreted as an array index, `arr[np.array(seq)]`, which will result either in an error or a different result.return np.add.reduce(sorted[indexer] * weights, axis=axis) / sumval
