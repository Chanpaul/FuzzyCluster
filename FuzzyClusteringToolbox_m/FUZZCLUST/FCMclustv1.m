function result = FCMclustv1(dataset,param)


%data normalization
X=dataset;
initValidityXB=10000000;
c=param.c;
clutsType=[];
%checking the parameters given
%default parameters
% if exist('param.m')==1, m = param.m;else m = 2;end;
% if exist('param.e')==1, e = param.e;else e = 1e-4;end;
m = param.m;
e = param.e;
[N,n] = size(X);
% [Nf0,nf0] = size(f0); 
%**********************Initialization I*********************
X1 = ones(N,1);
% Initialize fuzzy partition matrix
rand('state',0)
mm = mean(X);             %mean of the data (1,n)
aa = max(abs(X - ones(N,1)*mm)); %
v = 2*(ones(c,1)*aa).*(rand(c,n)-0.5) + ones(c,1)*mm;
for j = 1 : c,
    xv = X - X1*v(j,:);
    d(:,j) = sum((xv*eye(n).*xv),2);
end;
d0=d;
d = (d+1e-10).^(-1/(m-1));
f0 = (d ./ (sum(d,2)*ones(1,c)));
iter = 0;                       % iteration counter
initialJ=10000000000000000;
J0=sum(diag(transpose(f0.^m)*d0));
distout=zeros(N,c);
J=[];
%*************************End***********************
%******************Initialization II**************************
% X1 = ones(N,1);
% f0=initfcm(c,N);
% f0=f0';
% iter = 0;                       % iteration counter
% initialJ=10000000000000000;
% J0=0;
% distout=zeros(N,c);
% J=[];
%******************end****************************
%****************Initialization III --randomly select centroids********
% X1 = ones(N,1);
% rng('default');
% % centroidIndex1=randi(N,1,c);  %there are some repeatable integers.
% centroidIndex=randperm(N);   %Generate N unrepeatable integers, then we select c of them.
% v=zeros(c,n);
% for j = 1 : c,
%     v(j,:)=X(centroidIndex(j),:);
%     xv = X - X1*v(j,:);
%     d(:,j) = sum((xv*eye(n).*xv),2);
% end;
% d0=d;
% d = (d+1e-10).^(-1/(m-1));
% f0 = (d ./ (sum(d,2)*ones(1,c)));
% iter = 0;                       % iteration counter
% initialJ=10000000000000000;
% J0=sum(diag(transpose(f0.^m)*d0));
% distout=zeros(N,c);
% J=[];
%*****************end*****************************
%*********************test****************************
figure(c);
cmap=colormap;
plot(dataset(:,1),dataset(:,2),'.','MarkerSize',3,'MarkerFaceColor',[0 0 0]);
hold on
markerSz=8;
marker=['o','*','s','<','d'];
for k=1:c,
    ic=int8((k*64.)/(c*1.));
    color=cmap(ic,:);
    plot(v(k,1),v(k,2),marker(k),'MarkerSize',markerSz,'MarkerFaceColor',color);
    hold on
end;
%********************end test******************************

% Iterate

while abs(J0-initialJ)>e
  initialJ=J0;
  iter = iter + 1;
  % Calculate centers
  fm = f0.^m;
  sumf = sum(fm);
  v = (fm'*X)./(sumf'*ones(1,n));
  for j = 1 : c,
    xv = X - X1*v(j,:);
    d(:,j) = sum((xv*eye(n).*xv),2);
  end;
  distout=sqrt(d);
  d0=d;
  d = (d+1e-10).^(-1/(m-1));
  f0 = (d ./ (sum(d,2)*ones(1,c)));
%   J(iter) = sum(sum(f0.*d));
  J(iter)=sum(diag(transpose(f0.^m)*d0));
  J0=J(iter);
  %***********************test for illustrating the variation of cluster center*********************
  markerSz=markerSz+2;
  for k=1:c,
      ic=int8((k*64.)/(c*1.));
      color=cmap(ic,:);
      plot(v(k,1),v(k,2),marker(k),'MarkerSize',8,'MarkerFaceColor',color);
      hold on
  end;
  %*********************end of test****************************
  
end

%results
result.data.f=f0;
result.data.d=distout;
result.cluster.v=v;
result.cluster.type=clutsType;
result.iter = iter;
result.cost = J;