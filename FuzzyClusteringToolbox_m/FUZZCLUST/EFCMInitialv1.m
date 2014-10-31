function result = EFCMInitialv1(dist,dc)

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
dp(ordrho(1)).dcNeigh=find(dist(ordrho(1),:)<=dc); 
for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);          
     end
   end
   dp(nneigh(ordrho(ii))).dirNeigh(end+1)=ordrho(ii);
   dp(jj).dcNeigh=find(dist(jj,:)<=dc); 
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
%************************separate into dc neighorhood********************
clustDom=[];
clustDom(1).ele=[ordrho(1)];
clustDom(1).con=[];
clustDom(1).converged=[];
clustDomFeat=[0];
waitforcheck=[];
theleft=ordrho;
used=[];
virtualCl=zeros(1,ND);
index=1;
while ~isempty(theleft),    
    waitforcheck=[waitforcheck,theleft(1)];  
    if theleft(1)~=ordrho(1), 
        index=index+1;
        clustDom(index).ele=[theleft(1)];
        clustDom(index).converged=[];  
        clustDomFeat(index)=0;
    end;    
    theleft(1)=[]; 
    
    while ~isempty(waitforcheck),
        used=[used,waitforcheck(1)];
        temp=dp(waitforcheck(1)).dcNeigh;
        temp1=dp(waitforcheck(1)).dirNeigh;
        if ~isempty(temp1),
            intersectEle=intersect(temp,temp1);
            if ~isempty(intersectEle),
                clustDom(end).ele=unique([clustDom(end).ele,intersectEle]);                
                for j=1:size(intersectEle,2),
                    tt1=find(theleft==intersectEle(j));
                    if ~isempty(tt1),
                        theleft(tt1)=[];
                    end;  
                    tt2=find(used==intersectEle(j));
                    if isempty(tt2),
                        waitforcheck=[waitforcheck,intersectEle(j)];
                    end;                     
                end;
            end;          
            
        end;
        waitforcheck(1)=[];
        waitforcheck=unique(waitforcheck);
    end;
    if ~isempty(clustDom(index).ele),
        clustDomFeat(index)=max(rho(clustDom(index).ele));
    end;
    clustDom(index).converged=clustDom(index).ele;
    virtualCl(clustDom(index).ele)=index;
end;

[trho,ordIndex]=sort(clustDomFeat);
sizeOfClustDom=size(clustDom,2);
borderRho=zeros(1,sizeOfClustDom);
avgBorderRho=0.0;
totalConNum=0;
for jj=1:sizeOfClustDom,
    tsz=size(clustDom(jj).ele,2);
    for ii=1:tsz,
        tid=clustDom(jj).ele(ii);        
        for ss=1:size(dp(tid).dcNeigh,2);
            isCon=find(clustDom(jj).ele==dp(tid).dcNeigh(ss));
            if isempty(isCon),
                clustDom(jj).con(end+1)=dp(tid).dcNeigh(ss); 
                avgRho=(rho(dp(tid).dcNeigh(ss))+rho(tid))/2;
                avgBorderRho=avgBorderRho+avgRho;
                totalConNum=totalConNum+1;
                if avgRho>borderRho(jj),
                    borderRho(jj)=avgRho;
                end;
            end;
        end;        
    end;    
end;
% avgBorderRho=sum(borderRho)/sizeOfClustDom;
% figure(22);
% plot(1:sizeOfClustDom,borderRho,'o');

avgBorderRho1=avgBorderRho/totalConNum;
avgBorderRho=9;
% avgBorderRho1=max(borderRho(find(borderRho<max(rho))));
convergedDom=clustDom;
%******************Converge the Domain**********************
isExist=ones(1,sizeOfClustDom);
for jj=1:sizeOfClustDom,
    id=ordIndex(jj);
    temp1=clustDom(id).con;
    if ~isempty(temp1),        
        comEle=[];        
        for ii=1:sizeOfClustDom, 
            if ii~=jj && isExist(ordIndex(ii))==1,
                temp2=clustDom(ordIndex(ii)).converged;
%                 temp2=clustDom(ordIndex(ii)).ele;
                if ~isempty(temp2),
                    comEle=intersect(temp1,temp2);
                    if ~isempty(comEle),
                        
                        realCon=find(rho(comEle)>=avgBorderRho);
                        if size(realCon,2)>0,
                            clustDom(ordIndex(ii)).converged=unique([clustDom(ordIndex(ii)).converged,clustDom(id).converged]);
                            clustDom(id).ele=[];
                            clustDom(id).con=[];
                            clustDom(id).converged=[];
                            isExist(id)=0;    
                        end;                        
                    end;
                end;
            end;
        end;
        
    end;
end;
ss=find(isExist==0);
if ~isempty(ss),
    clustDom(ss)=[];
end;
%**************************End of Converge*******************************

sortedClustSeq=[];

temprho=[];
cordp=[];
for j=1:size(clustDom,2),
    [ttrho,ttordind]=max(rho(clustDom(j).converged)); 
    temprho(j)=ttrho;
    cordp(j)=clustDom(j).converged(ttordind);    
end;
[sorted_temprho,ordindex]=sort(temprho,'descend');
% figure(12);
% plot(1:size(clustDom,2),sorted_temprho,'o');

if size(sorted_temprho,2)==1,
    sorted_temprho(end+1)=2;
end;
sortedClustSeq=cordp(ordindex);
finalClustDom=clustDom(ordindex);
sortedClustSeq=cordp(ordindex);  
diff=zeros(1,size(sortedClustSeq,2));
% for jj=2:size(sortedClustSeq,2)-1,
%     diff(jj)=rho(sortedClustSeq(jj))-rho(sortedClustSeq(jj+1));
% end;
% [d,id]=max(diff);
% sortedClustSeq=sortedClustSeq(1:id);


%******************Revise the cluster domain*************************


%**********************end of revise*************************

%****************************End***********************************
    
%*********************************The end********************************


% result.cl=cl;
result.rho=rho;
result.delta=delta;
% result.gama_sorted=gama_sorted;
% result.ordgama=ordgama;
result.nneigh=nneigh;
result.rho_sorted=rho_sorted;
result.ordrho=ordrho;
result.sortedClustSeq=sortedClustSeq;
result.finalClustDom=finalClustDom;

end




