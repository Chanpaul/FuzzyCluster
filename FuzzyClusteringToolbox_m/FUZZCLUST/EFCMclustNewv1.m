function result = EFCMclustNewv1(dataset,param,iniParam)
%data normalization

X=dataset;
clusterNum=param.c;    
dist=param.dist;
% dc=param.dc;
m = param.m;
e = param.e;
dc=param.dc;
% f0=param.c;
%checking the parameters given
%default parameters
% if exist('param.m')==1, m = param.m;else m = 2;end;
% if exist('param.e')==1, e = param.e;else e = 1e-4;end;

[N,n] = size(X);
% [Nf0,nf0] = size(f0); 
% X0 = ones(N,1);


% ordgama=iniParam.ordgama;
nneigh=iniParam.nneigh;
rho_sorted=iniParam.rho_sorted;
ordrho=iniParam.ordrho;
icl=zeros(clusterNum,n);
% cl=iniParam.cl;
rho=iniParam.rho;
clustSeq=iniParam.sortedClustSeq;
clustDom=iniParam.finalClustDom;
delta=iniParam.delta;
% [delta_sorted,orddelta]=sort(delta,'descend');
d=zeros(N,clusterNum);

normRho=rho./max(rho);
clustDomMx=zeros(N,clusterNum);

%**********************************************
cl=-1*ones(1,N);
for t=1:clusterNum,
    cl(clustSeq(t))=t;
end;
for i=1:N,
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));
  end
end

%************************************************

for i=1:clusterNum,  
    domEle=find(cl==i);
%     clustDomMx(domEle,i)=1;
    icl(i,:)=mean(dataset(domEle,:));
    clustDomMx(clustDom(i).converged,i)=1;
%     icl(i,:)=mean(dataset(clustDom(i).converged,:));
    for j=1:N,
        d(j,i)=sqrt(sum((dataset(j,:)-icl(i)).^2,2));        
    end;    
end;

rhoCoef=clustDomMx.*(normRho'*ones(1,clusterNum));
%********************Computing the initial membership********************

d0=d;
d = (d+1e-10).^(-1/(m-1));
f0 = (d ./ (sum(d,2)*ones(1,clusterNum)));
% f0=0.8*f0+0.2*clustDomMx;
% for j=1:N,
%     tt=find(clustDomMx(j,:)==1);
%     if ~isempty(tt),
%         f0(j,:)=0.7.*f0(j,:)+0.3*clustDomMx(j,:);
%     end;
% end;
J0=sum(diag(transpose(rhoCoef.*f0.^m)*d0));
% J0=sum(diag(transpose(f0.^m)*d0));
% J0=sum(diag(transpose(rhoCoef.*f0.^m)*d0));
%********************End of computing initial membership******************

% f = zeros(N,clusterNum);                % partition matrix
iter = 1;                       % iteration counter
initialJ=10000000;
% J0=sum(diag(transpose(f0.^m)*d));
X1 = ones(N,1);
distout=d;
v=zeros(clusterNum,n);
J=[];

figure(clusterNum);
cmap=colormap;
plot(dataset(:,1),dataset(:,2),'.','MarkerSize',3,'MarkerFaceColor',[0 0 0]);
hold on
marker=['o','*','s','<','d'];
for k=1:clusterNum,
    ic=int8((k*64.)/(clusterNum*1.));
    color=cmap(ic,:);
    plot(icl(k,1),icl(k,2),marker(k),'MarkerSize',8,'MarkerFaceColor',color);
    hold on
end;

while abs(J0-initialJ)>e
  initialJ=J0;
  iter = iter + 1;
  % Calculate centers
  
  fm = f0.^m;
%   sumf = sum(fm);
%   v = (fm'*X)./(sumf'*ones(1,n));
%   v = (fm'*X)./(sumf'*ones(1,n));

  for j = 1 : clusterNum,
      clustElem=find(cl==j);      
% %       clustElem=clustDom(j).converged;      
      sumf=sum(fm(clustElem,j));
      v(j,:)= fm(clustElem,j)'*X(clustElem,:)./sumf;
      xv = X - X1*v(j,:);
      d(:,j) = sum((xv*eye(n).*xv),2);      
  end;
  distout=sqrt(d);
  d0=d;
  d = (d+1e-10).^(-1/(m-1));
  f0 = (d ./ (sum(d,2)*ones(1,clusterNum)));
%   f0=0.8*f0+0.2*clustDomMx;
  %*******************revise the fuzzy membership*************************
%   for j=1:N,
%       tt=find(clustDomMx(j,:)==1);
%       if ~isempty(tt),
%           f0(j,:)=0.7.*f0(j,:)+0.3*clustDomMx(j,:);
%       end;
%   end;
  %**********************end of revise************************
%   f0=0.4.*f0+0.6.*f_density;
%   J(iter) = sum(sum(f0.*d));
%   J(iter)=sum(diag(transpose(f0.^m)*d0));
  
  J(iter)=sum(diag(transpose(rhoCoef.*f0.^m)*d0));
  J0=J(iter);
  
  for k=1:clusterNum,
      ic=int8((k*64.)/(clusterNum*1.));
      color=cmap(ic,:);
      plot(v(k,1),v(k,2),marker(k),'MarkerSize',8,'MarkerFaceColor',color);
      hold on
  end;
end
% hold off



% fm = f.^m; 
% sumf = sum(fm);

%results
result.data.f=f0;
result.data.d=distout;
result.cluster.v=v;
result.iter = iter-1;
result.cost = J;
