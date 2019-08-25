function [local_D,local_SD] = LocalDist(X,k,d)
% calculate pairwise distances between neighbors
% input: X = data set
%           k = number of neighbors in k-nn
%           d = intrinsic dimenison
% output: local_D: locally Euclidean pairwise distance matrix
%              local_SD: locally spherical pairwise distance matrix
%
% History:
%   Didong Li       June 1, 2018, created



% calculating global pairwise distance matrix to find neighbors
D = pdist(X);
D = squareform(D);

[n,m] = size(X); % n is the sample size



% initialize distance matrices
L = (max(D(:))+1)*100; % can be replaced by larger numbers depending on the scale of the data. 
local_D = L*ones(n,n);  % for non-neighbors, the distance is assigned to be L representing infinity.
local_SD = L*ones(n,n);


for i = 1:n
    [out,ind] = sort(D(i,:));
    kind = ind(1:1+k);
    Xcls = X(kind,:); % neighbors of the i-th sample
    [c,V,r] = SPCA(Xcls,d); % find spherelet by SPCA
    Ycls = zeros(k+1,m); 
    curY=c.'+r*(X(i,:)-c.')*V*V.'/norm(V*V.'*(X(i,:).'-c));
    for j = 1:k+1
        Ycls(j,:) = c.'+r*(Xcls(j,:)-c.')*V*V.'/norm(V*V.'*(Xcls(j,:).'-c)); % project the points to the spherelets
        local_SD(i,ind(j)) = r*acos((curY-c.')*(Ycls(j,:).'-c)/r^2); % local_SD is the locally spherical distance
        local_D(i,ind(j)) = D(i,ind(j)); % local_D is the locally Euclidean distance
    end
    local_SD(i,i) = 0; % set the diagonal element to be zero
    local_D(i,i) = 0;
end
%local_SD = (local_SD + local_SD.')/2;  % symmetrize both matrices.
%local_D = (local_D + local_D.')/2;

local_SD = min(local_SD, local_SD.');  % symmetrize both matrices.
local_D = min(local_D, local_D.');

return
