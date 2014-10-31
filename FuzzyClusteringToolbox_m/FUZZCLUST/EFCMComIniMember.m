function result = EFCMComIniMember( iniResult,clusterNum,m,dataNum,dist)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cl=iniResult.cl;
rho=iniResult.rho;
delta=iniResult.delta;
gama_sorted=iniResult.gama_sorted;
ordgama=iniResult.ordgama;
nneigh=iniResult.nneigh;
rho_sorted=iniResult.rho_sorted;
ordrho=iniResult.ordrho;
icl=zeros(1,clusterNum);
ND=dataNum;
for i=1:clusterNum,
    cl(ordgama(i))=i;
    icl(i)=ordgama(i);
end;

% fprintf('NUMBER OF CLUSTERS: %i \n', clusterNum);
% disp('Performing assignation')

%assignation
% for i=1:ND
%   if (cl(ordrho(i))==-1)
%     cl(ordrho(i))=cl(nneigh(ordrho(i)));
%   end
% end
% %halo
% for i=1:ND
%   halo(i)=cl(i);
% end
% for i=1:clusterNum
%     bord_rho(i)=0.;
% end;
% for i=1:ND-1,
%     for j=i+1:ND,
%         if ((cl(i)~=cl(j))&& (dist(i,j)<=dc)),
%             rho_aver=(rho(i)+rho(j))/2.;
%             if (rho_aver>bord_rho(cl(i))) 
%                 bord_rho(cl(i))=rho_aver;
%             end;
%             if (rho_aver>bord_rho(cl(j))) ,
%                 bord_rho(cl(j))=rho_aver;
%             end;
%         end;
%     end;
% end;
% for i=1:ND
%     if (rho(i)<bord_rho(cl(i)))
%         halo(i)=0;
%     end;
% end;

%*****************initialize the membership************************
J=0;
initiald=zeros(ND,clusterNum);
for j = 1 : clusterNum,
    clusterCenter=icl(j);    
    initiald(:,j) = dist(:,clusterCenter);
end;
initiald = (initiald+1e-10).^(-1/(m-1));
f0 = (initiald ./ (sum(initiald,2)*ones(1,clusterNum)));
J=sum(diag(transpose(f0.^m)*initiald));


%*****************************************


% for i=1:ND,  
%   for j=1:clusterNum,
%     mu(i,j)=0.0;
%     if (halo(i)==j),
%         mu(i,j)=unifrnd(0.80,0.99);
%     elseif (halo(i)==0),
%         mu(i,j)=unifrnd(0.01,0.05);
%     else
%         mu(i,j)=unifrnd(0.10,0.20);
%     end;
%     
%   end;
%   
% end;
% mu=mu./(sum(mu,2)*ones(1,clusterNum));
% for i=1:clusterNum
%   nc=0;
%   nh=0;
%   for j=1:ND
%     mu(i,j)=0.0;
%     if (cl(j)==i) 
%       nc=nc+1;
%     end
%     if (halo(j)==i) 
%       nh=nh+1;
%     end;
%   end
%   fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
% end



result.v=icl;
result.f0=f0;
result.cl=cl;
% result.halo=halo;
result.d=initiald;
result.J=J;

end

