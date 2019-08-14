# GeodesicDistance

Created by Didong Li, June 2018

The function `GeoDist` is designed to estimate the geodesic distance between samples assuming the data lie in some d-dimensional Riemannian manifold embedded in R^p
See more details about the spherical distance and spherelets here:
https://arxiv.org/abs/1706.08263

Input: 
- X = data set, represented by n by p matrix where n is the sample size, p is the number of features
- k = number of neighbors used in k-nearest neighbors
- d = intrinsic dimension

Output: 
- D = Euclidean pairwise distance matrix
- SD = spherical pairwise distance matrix

functions `pca`, `SPCA` and `LocalDist` are all necessary to run `GeoDist`.
