function result = EFCMclustNewv1(dataset,param,iniParam)
%data normalization

X=dataset;
clusterNum=param.c;    
dist=param.dist;

m = param.m;
e = param.e;
dc=param.dc;

[N,n] = size(X);

nneigh=iniParam.nneigh;
rho_sorted=iniParam.rho_sorted;
ordrho=iniParam.ordrho;
dclMx=iniParam.dclMx;   %dclMx(i,j)=1 means that the ith data belongs to the cluster centered at the jth data.
clustRt=iniParam.pclustRt;  %  iniParam.clustRt

delta=iniParam.delta;
icl=zeros(clusterNum,n);
rho=iniParam.rho;

d=zeros(N,clusterNum);
clustRtSz=size(clustRt,2);
if clusterNum<clustRtSz,
    [sortedRtDelta,origInd]=sort(delta(clustRt),'descend');
    for j=1:clustRtSz-clusterNum,
        tempRt=clustRt(clustRtSz-j+1);
        dclMx(:,nneigh(tempRt))=dclMx(:,nneigh(tempRt))+dclMx(:,tempRt);
        dclMx(:,tempRt)=0;
    end;
    clustRt=clustRt(origInd(1:clusterNum));
elseif clusterNum>clustRtSz,
    disp('the input cluster number exceeds the optimized cluster number');
    clusterNum=clustRtSz;
end;
%normRho=rho./max(rho);
for i=1:clusterNum, 
    icl(i,:)=full((transpose(dclMx(:,clustRt(i)))*dataset))./sum(dclMx(:,clustRt(i)));
    for j=1:N,
        d(j,i)=sum((dataset(j,:)-icl(i)).^2,2);        
    end;    
end;
%Introducing null cluster
%d(:,end+1)=transpose(max(dist).^2);
%
% rhoCoef=clustDomMx.*(normRho'*ones(1,clusterNum));
%********************Computing the initial membership********************

d0=d;
d = (d+1e-10).^(-1/(m-1));
f0 = (d ./ (sum(d,2)*ones(1,clusterNum)));

%J0=sum(diag(transpose(rhoCoef.*f0.^m)*d0));
J0=sum(diag(transpose(f0.^m)*d0));
%f0=f0(:,1:clusterNum); %exclude the null cluster.
%********************End of computing initial membership******************

% f = zeros(N,clusterNum);                % partition matrix
iter = 1;                       % iteration counter
initialJ=10000000;
% J0=sum(diag(transpose(f0.^m)*d));
X1 = ones(N,1);
distout=d;
v=zeros(clusterNum,n);
J=[];

% figure(clusterNum);
% cmap=colormap;
% plot(dataset(:,1),dataset(:,2),'.','MarkerSize',3,'MarkerFaceColor',[0 0 0]);
% hold on
% marker=['o','*','s','<','d'];
% for k=1:clusterNum,
%     ic=int8((k*64.)/(clusterNum*1.));
%     color=cmap(ic,:);
%     plot(icl(k,1),icl(k,2),marker(k),'MarkerSize',8,'MarkerFaceColor',color);
%     hold on
% end;

while abs(J0-initialJ)>e
  initialJ=J0;
  iter = iter + 1;
  % Calculate centers
  
  fm = f0.^m;
  sumf = sum(fm);
  v = (fm'*X)./(sumf'*ones(1,n));
%   v = (fm'*X)./(sumf'*ones(1,n));

  for j = 1 : clusterNum,
%       sumf=transpose(dclMx(:,clustRt(j)))*fm(:,j);
% %       clustElem=find(cl==j);      
% %     
% %       sumf=sum(fm(clustElem,j));
% %       v(j,:)= fm(clustElem,j)'*X(clustElem,:)./sumf;
%       v(j,:)= full(transpose(dclMx(:,clustRt(j)).*fm(:,j))*X./sumf);
      xv = X - X1*v(j,:);
      d(:,j) = sum((xv*eye(n).*xv),2);      
  end;
  distout=sqrt(d);  
  d0=d;
  d = (d+1e-10).^(-1/(m-1));
  f0 = (d ./ (sum(d,2)*ones(1,clusterNum)));
%   J(iter) = sum(sum(f0.*d));
  J(iter)=sum(diag(transpose(f0.^m)*d0));
  
%   J(iter)=sum(diag(transpose(rhoCoef.*f0.^m)*d0));
  J0=J(iter);
  %f0=f0(:,1:clusterNum); %exclude the null cluster.
%   for k=1:clusterNum,
%       ic=int8((k*64.)/(clusterNum*1.));
%       color=cmap(ic,:);
%       plot(v(k,1),v(k,2),marker(k),'MarkerSize',8,'MarkerFaceColor',color);
%       hold on
%   end;
end

result.data.f=f0;
result.data.d=distout;
result.cluster.v=v;
result.iter = iter-1;
result.cost = J;
