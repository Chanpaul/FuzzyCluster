function m = ComputeM( data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[N,n]=size(data);

%***************standardize the dataset********************
xbar=sum(data)/N;
diff=data-ones(N,1)*xbar;
s_dataset=diff./sqrt((ones(N,1)*sum((diff.^2)./(N-1))));
%*****************************End**************************
% dataset=data;
dataset=s_dataset;
xbar=sum(dataset)/N;
diff=dataset-ones(N,1)*xbar;
d=sqrt(sum(diff.^2,2));
% F0=dataset-(1./d)*xbar;
F0=diff.*((1./d)*ones(1,n));
F1=(transpose(F0)*F0)./N;
lambda=eig(F1);
lambda_max=max(lambda);
if lambda_max<0.5,
    t_m=1/(1-2*lambda_max);
else
    t_m=Inf;
end;
m=t_m;
end

