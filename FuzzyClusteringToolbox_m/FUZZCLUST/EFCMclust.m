function result = EFCMclust(dataset,param,iniParam)
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


ordgama=iniParam.ordgama;
nneigh=iniParam.nneigh;
rho_sorted=iniParam.rho_sorted;
ordrho=iniParam.ordrho;
icl=zeros(1,clusterNum);
cl=iniParam.cl;
d=zeros(N,clusterNum);
for i=1:clusterNum,
    cl(ordgama(i))=i;
    icl(i)=ordgama(i);
    d(:,i)= dist(:,icl(i));
end;



disp('Performing assignation')

%assignation
for i=1:N
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));
  end
end
%halo
for i=1:N
  halo(i)=cl(i);
end
if (clusterNum>1)
  for i=1:clusterNum
    bord_rho(i)=0.;
  end
  for i=1:N-1
    for j=i+1:N
      if ((cl(i)~=cl(j))&& (dist(i,j)<=dc))
        rho_aver=(rho(i)+rho(j))/2.;
        if (rho_aver>bord_rho(cl(i))) 
          bord_rho(cl(i))=rho_aver;
        end
        if (rho_aver>bord_rho(cl(j))) 
          bord_rho(cl(j))=rho_aver;
        end
      end
    end
  end
  for i=1:N
    if (rho(i)<bord_rho(cl(i)))
      halo(i)=0;
    end
  end
end

%********************Computing the initial membership********************
% cl=iniParam.cl;
% rho=iniParam.rho;
% delta=iniParam.delta;
% gama_sorted=iniParam.gama_sorted;
% ordgama=iniParam.ordgama;
% nneigh=iniParam.nneigh;
% rho_sorted=iniParam.rho_sorted;
% ordrho=iniParam.ordrho;
% icl=zeros(1,clusterNum);
% d=zeros(N,clusterNum);
% for i=1:clusterNum,
% %     cl(ordgama(i))=i;
%     icl(i)=ordgama(i);
%     d(:,i)= dist(:,icl(i));
% end;
figure(5);
plot(dataset(icl,1),dataset(icl,2),'ro');
hold on
plot(dataset(:,1),dataset(:,2),'b.');
hold off;
d0=d;
d = (d+1e-10).^(-1/(m-1));
f0 = (d ./ (sum(d,2)*ones(1,clusterNum)));
J0=sum(diag(transpose(f0.^m)*d0));
%********************End of computing initial membership******************

% f = zeros(N,clusterNum);                % partition matrix
iter = 0;                       % iteration counter
initialJ=10000000;
% J0=sum(diag(transpose(f0.^m)*d));
X1 = ones(N,1);
distout=d;
v=zeros(clusterNum,n);
J=[];
% Iterate
% while abs(J0-initialJ)>e
%   initialJ=J0;
%   iter = iter + 1;
% %   f = f0;
%   fm=f0.^m;
% %   fm=rho'*ones(1,clusterNum).*(f.^m);
%   % Calculate centers
% %   fm = f.^m;
% 
%   sumf = sum(fm);
%   v = (fm'*X)./(sumf'*ones(1,n));
%   for j = 1 : clusterNum,
%     xv = X - X1*v(j,:);
%     d(:,j) = sum((xv*eye(n).*xv),2);
%   end;
%   distout=sqrt(d);
%   d = (d+1e-10).^(-1/(m-1));
%   f0 = (d ./ (sum(d,2)*ones(1,clusterNum)));
%   fm=f0.^m;
% %   J(iter) = sum(sum(f0.*d));
% %   J(iter)=sum(diag(transpose(f0.^m)*d));
%   J(iter)=sum(diag(transpose(fm)*d));
%   J0=J(iter);
%   % Update f0
%   
%    
% end
figure(6);
plot(dataset(:,1),dataset(:,2),'b.');
hold on
while abs(J0-initialJ)>e
  initialJ=J0;
  iter = iter + 1;
  % Calculate centers
  fm = f0.^m;
  sumf = sum(fm);
  v = (fm'*X)./(sumf'*ones(1,n));
  for j = 1 : clusterNum,
    xv = X - X1*v(j,:);
    d(:,j) = sum((xv*eye(n).*xv),2);
  end;
  distout=sqrt(d);
  d0=d;
  d = (d+1e-10).^(-1/(m-1));
  f0 = (d ./ (sum(d,2)*ones(1,clusterNum)));
%   J(iter) = sum(sum(f0.*d));
  J(iter)=sum(diag(transpose(f0.^m)*d0));
  J0=J(iter);
  plot(v(:,1),v(:,2),'ro');
end




% fm = f.^m; 
% sumf = sum(fm);

%results
result.data.f=f0;
result.data.d=distout;
result.cluster.v=v;
result.iter = iter;
result.cost = J;
