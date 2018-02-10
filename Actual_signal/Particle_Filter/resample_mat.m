function [ xk_resampled ] = resample_mat( xk,wk )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N_s=length(wk);

c=zeros(N_s,1);
u=zeros(N_s,1);
xk_resampled=zeros(size(xk));

%%% CDF of the current weights vector (=distribution)

for i=2:N_s
    c(i)=c(i-1)+wk(i); % Cumulative distribution function (CDF)
end

%%%% Comparing the input vector CDF to a uniform CDF and taking
%%%% N_s significant terms. we may very well take few times the same term
%%%% (if its weight is very high)
u(1)=1/N_s; 
i=1;
for j=1:N_s
    
    u(j)=u(1)+1/N_s*(j-1); % new CDF of uniform distribution. 
    
    while c(i) < u(j) && i<N_s % find the index of the first significant term in the CDF wrt to u(j) 
        i=i+1; % we enter here if c(i) is big enough; if the probability of xk(i) is significant
    end
    
    xk_resampled(:,j)=xk(:,i); % assign the term with that significant weight to the new vector
    
end



end

