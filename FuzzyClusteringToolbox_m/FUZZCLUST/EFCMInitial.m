function result = EFCMInitial(dist,dc)

ND=size(dist,1);
% fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
rho=zeros(1,ND);
delta=zeros(1,ND);
gamma=zeros(1,ND);
nneigh=zeros(1,ND);

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
for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);  
        
   end
end
delta(ordrho(1))=max(delta(:));
figure(10);
plot(rho,delta,'o');
xlabel('\rho');
ylabel('\delta');



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
sortedClustSeq=[];
for j=1:size(clustSeq,2),
    sortedClustSeq(j)=clustSeq(ordinx(j));
end;


% cl=-1*ones(1,ND);
% for i=1:ND
%   cl(i)=-1;
% end;
% result.cl=cl;
result.rho=rho;
result.delta=delta;
% result.gama_sorted=gama_sorted;
% result.ordgama=ordgama;
result.nneigh=nneigh;
result.rho_sorted=rho_sorted;
result.ordrho=ordrho;
result.sortedClustSeq=sortedClustSeq;

end




