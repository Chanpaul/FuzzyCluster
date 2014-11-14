function [ output_args ] = sammon_mapping(dataset,dataLabels,clustRes,sammonMapName,m)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
colors={'r.' 'b*' 'gx' 'b+' 'ys' 'm.' 'c.' 'k.' 'r*' 'g*' 'y*' 'm*' 'c*' 'k*' };
data.X=dataset; %only features.
C=dataLabels;
result=clustRes;
param.m=m;
% C=zeros(length(data),1);
% for i=1:3
%     C(find(iris(:,end)==i))=i;
% end
% % %normalization of the data
% data.X=data;
% data=clust_normalize(data,'range');
%fuzzy c-means clustering 
% param.m=2;
% param.c=3;  %2
% param.val=1;
% param.vis=0;
% result=Kmeans(data,param);
% result=validity(result,data,param);
%Assignment for classification
[d1,d2]=max(result.data.f');
Cc=[];
param.c=size(clustRes.cluster.v,1);
for i=1:param.c
    Ci=C(find(d2==i));
    dum1=hist(Ci,1:param.c);
    [dd1,dd2]=max(dum1);
    Cc(i)=dd2;
end
%Principal Component Projection of the data and the cluster centers
param.q=2;
result = PCA(data,param,result); 
%SAMMON mapping

proj.P=result.proj.P;   %Sammon uses the output of PCA for initializing
param.alpha = 0.4;
param.max=100;

%figure(2) Sammonv1
result = Sammon(proj,data,result,param);
%Modified fuzzy SAMMON mapping
  
proj.P=result.proj.P; %FuzSam uses the output of Sammon for initializing
param.alpha = 0.4;
param.max=100;

% figure(3)
simmonMapHandle=figure;
result=FuzSam(proj,result,param);

clf
for i=1:max(C)+1
    index=find(C==i-1);
    err=(Cc(d2(index))~=i-1);
    eindex=find(err);
    %misclass(i-1)=sum(err);
    plot(result.proj.P(index,1),result.proj.P(index,2),[colors{i}] )
    hold on
    plot(result.proj.P(index(eindex),1),result.proj.P(index(eindex),2),'o','MarkerSize',10);
    hold on
end    
    xlabel('y_1')
    ylabel('y_2')
    title('Fuzzy Sammon mapping')
    
plot(result.proj.vp(:,1),result.proj.vp(:,2),'r*')
%calculating realtion-indexes
result = samstr(data,result);
perff = [projeval(result,param) result.proj.e];
saveas(simmonMapHandle,sammonMapName,'png');
close(simmonMapHandle);
end

