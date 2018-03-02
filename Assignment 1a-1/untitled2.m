clearvars
clear all
indx = 4000;
nd   = 16;
K    = 40;


% To return indexing of DV in math model from CPLEX index
i = floor((indx-1)/(nd*K))+1;
j = floor(((indx-(i-1)*nd*K)-1)/K)+1;
k = mod(indx,K);
if k == 0
    k = K;
end


Xindex(i,j,k,nd,K) 

function out = Xindex(m,n,p,nd,k)
    out =  (m - 1)*nd*k + (n-1)*k + p;  % Function given the variable index for each X(i,j,k) [=(m,n,p)]  
end