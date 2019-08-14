function [c,V,r]=SPCA(X,d)
% Find the best  d-dimensional sphere to fit data X
% input:  X = data matrix
%            d = intrinsic dimension
% output: c = center of the spherelet
%              r = radius of the spherelet
%              V = the subspace where the sphere lies in the affine space c+V
%
% History:
%   Didong Li       June 1, 2018, created


[n,m]=size(X); % n=sample size, m = ambient dimension

if n>d+1   % if there are enough samples, fit the data by a sphere or a hyperplane
    
 % do d+1 dimensional PCA first
[V,mu] = pca(X,d+1); 
Y=ones(n,1)*mu+(X-ones(n,1)*mu)*V*V.'; % projection of X onto the d+1 dimensional affine space mu+V
l=zeros(n,1);
for i=1:n
    l(i)=norm(Y(i,:))^2;
end
lbar=mean(l);
H=zeros(m,m);
f=zeros(m,1);
for i=1:n
    H=H+(mu-Y(i,:)).'*(mu-Y(i,:));
    f=f+(l(i)-lbar)*((mu-Y(i,:)).');
end
c=mu.'+V*V.'*(-0.5*pinv(H)*f-mu.');  % center of the sphere

Riemd=zeros(n,1); % distance between sample and center
for i=1:n
    Riemd(i)=norm((c.'-Y(i,:))); 
end
r=mean(Riemd); % radius

%SS=n*var(Riemd);
% SS=0;
% for i=1:n
%    % SS=SS-norm((c.'-X(i,:))*V)^2+norm(c.'-X(i,:))^2;
%    SS=SS+norm(X(i,:)-Y(i,:))^2; %first part of the error (linear projection)
% end
% SS=SS+n*var(Riemd); % second part of the error(spherical error)
%-------this is the sum of square error of the sphere



%-------then we do d dimensional PCA to find the error of a hyperplane and
%compare them
% [coeff_new,score_new,latent_new,tsquared_new,explained_new,mu_new] = pca(X); 
% V_new=coeff(:,1:d);
%Sum of square error of d dimensional PCA is the sum of the last m-d singular values 
% SS_new=0;
% for i=d+1:m
%     SS_new=SS_new+explained_new(i,1);
% end
% compare the two errors
% if SS>SS_new % hyperplane is better
%     c=mu;
%     V=V_new;
%     r=inf;
%     SS=SS_new;
% end 

% if there is no enough sample
else 
    c=mean(X,1);
    r=norm(X(1,:)-c);
    V=ones(m,d+1);
    display('no enough samples')
end

return

