function [V,mu] = pca(X,d)
% principal component analysis       
% Input: X = data set with size n by p, n is the sample size, p is the
%                   ambient dimension
%           d = dimension of the subspace
% 
% Output: V = a p by d matrix which stores principle components in columns.
%              mu = mean of samples
% 
% History:
%   Didong Li       June 1, 2018, created

[n,p] = size(X);          % n is the sample size, p is the ambient dimension

% default value of d is min(n,p)
if ~exist('d','var') || isempty(d)
    d = min(n,p); 
end; 


mu = mean(X,1);           % mu is the mean
X = X-ones(n,1)*mu;     % centralize the data with zero mean             

if p<=n
    S = X'*X;           % p by p covariance matrix.
else
    S = X*X';           % n by n Gram matrix.
end

[V,lambda] = eig(S);      % eigen decomposition of S
lambda = diag(lambda);    % diagonal elements of lamba are eigen values of S
[lambda,I]=sort(lambda,'descend');  % sort the eigenvalues in descending order
V = V(:,I(1:d));          % take the first d eigenvectors

return



