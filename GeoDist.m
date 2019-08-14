function [D,SD] = GeoDist(X,k,d)
% find the pairwise distance matrix
% input: X = data set
%           k = number of neighbors in k-nn
%           d = intrinsic dimension
% output: D = Euclidean pairwise distance matrix
%              SD = spherical pairwise distrance matrix
%
% History:
%   Didong Li       June 1, 2018, created



n = size(X,1); % n is sample size
% pairwise distance matrix between neighbors
% D is the locally Euclidean distance
% SD is the locally spherical distance
[D,SD] = LocalDist(X,k,d);  % first calculate distance between k-nearest neighbors


% find the global distance from shortest path search on graphs by Floyd's algorithm 
for i = 1:n
     D = min(D , repmat(D(:,i), [1 n]) + repmat(D(i,:) , [n 1]));  
     SD = min(SD , repmat(SD(:,i) , [1 n]) + repmat(SD(i,:) , [n 1]));  
end

return
