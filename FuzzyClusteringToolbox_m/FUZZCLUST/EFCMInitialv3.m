function result = EFCMInitialv3(dist,dc)
%determine the parameters in inilization with optimizational methods.
ND=size(dist,1);
% fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
rho=zeros(1,ND);
delta=zeros(1,ND);
gamma=zeros(1,ND);
nneigh=zeros(1,ND);
dp=[];
dp(ND).dirNeigh=[];
dp(ND).dcNeigh=[];

% for i=1:ND
%   rho(i)=0.;
% end
%
% Gaussian kernel
%


for i=1:ND-1    
  for j=i+1:ND
     rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
     rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));     
  end 
end
%
% "Cut off" kernel
%
% for i=1:ND-1
%  for j=i+1:ND
%    if (dist(i,j)<dc)
%       rho(i)=rho(i)+1.;
%       rho(j)=rho(j)+1.;
%    end
%  end
% end

maxd=max(max(dist));

[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;
potClustRt=[ordrho(1)];
%dcl=zeros(1,ND);
dp(ordrho(1)).dcNeigh=find(dist(ordrho(1),:)<=dc); 
dclMx=zeros(ND,ND);  %indicating the relationship of data point and the cluster.
dclMx(ordrho(1),ordrho(1))=1;
%dcl(ordrho(1))=ordrho(1); 
for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);          
     end
   end
   dp(nneigh(ordrho(ii))).dirNeigh(end+1)=ordrho(ii);
   dp(ordrho(ii)).dcNeigh=find(dist(ordrho(ii),:)<=dc); 
   if dist(nneigh(ordrho(ii)),ordrho(ii))>dc,
       potClustRt=[potClustRt,ordrho(ii)];
       %dcl(ordrho(ii))=ordrho(ii);
       dclMx(ordrho(ii),ordrho(ii))=1;
   else
       %dcl(ordrho(ii))=dcl(nneigh(ordrho(ii)));
       dclMx(ordrho(ii),nneigh(ordrho(ii)))=1;
   end;
end
delta(ordrho(1))=max(delta(:));
% figure(10);
% plot(rho,delta,'o');
% xlabel('\rho');
% ylabel('\delta');

%*************1ST method to compute initial centroids*****************
rhoth=rho(ordrho(1))/10;   %7
deltath=delta(ordrho(1))/10;   %7
clustSeq=[ordrho(1)];
distDg=[];
for i=1:ND,
  if ( (rho(i)>rhoth) && (delta(i)>deltath))
      if i~=ordrho(1),          
          clustSeq(end+1)=i;
          distDg(end+1)=sqrt((rho(i)-rho(ordrho(1)))^2+(delta(i)-delta(ordrho(1)))^2);
      end;     
  end;
end;

[sorted_clustCent_delta,ordinx]=sort(delta(clustSeq),'descend');
sortedClustSeq1=[];
for j=1:size(clustSeq,2),
    sortedClustSeq1(j)=clustSeq(ordinx(j));
end;
%****************************end*************************************
%***************I. finding the initial clustering centers****************
fClustRt=[];  %final clustering centers;
fdclMx=dclMx;  %final clustering relations;
potClustSz=size(potClustRt,2);
xbv1=-1*ones(1,potClustSz);
minXbv1=100000000;
[sortedRtDelta,origInd]=sort(delta(potClustRt),'descend');
tempDclMx=dclMx;
for i=1:potClustSz-1,
    tempClust=potClustRt(origInd(1:potClustSz-i));
    tempRt=potClustRt(origInd(potClustSz-i+1));
    tempDclMx(:,nneigh(tempRt))=tempDclMx(:,nneigh(tempRt))+tempDclMx(:,tempRt);
    tempDclMx(:,tempRt)=0;
    minDelta=sortedRtDelta(potClustSz-i);
    tempRes=dist.^2.* tempDclMx;   
    tempXBv1=sum(tempRes(:))/(ND*minDelta);
    xbv1(potClustSz-i)=tempXBv1;
    if tempXBv1<minXbv1,
        minXbv1=tempXBv1;
        fClustRt=tempClust;
        fdclMx=tempDclMx;
    end;
    
end;
    
%*********************************The end********************************


% result.cl=cl;
result.rho=rho;
result.delta=delta;
result.dclMx=fdclMx;
result.clustRt=fClustRt;
% result.gama_sorted=gama_sorted;
% result.ordgama=ordgama;
result.nneigh=nneigh;
result.rho_sorted=rho_sorted;
result.ordrho=ordrho;
% result.sortedClustSeq=sortedClustSeq;
% result.finalClustDom=finalClustDom;

end




