function result = myvalidity(result,data,param)
%validation of the clustering

N = size(result.data.f,1);
c = size(result.cluster.v,1);
n = size(result.cluster.v,2);
v = result.cluster.v;

% if exist('param.m')==1, m = param.m;else m = 2;end;
m = param.m;
if strcmp(param.val,'PC')
        %partition coefficient (PC)
        fm = (result.data.f).^m;
        PC = 1/N*sum(sum(fm));
        %classification entropy (CE)
        fm = (result.data.f).*log(result.data.f);
        CE = -1/N*sum(sum(fm));
     
     %results   
     result.validity.PC = PC;
     result.validity.CE = CE;        
        
        
 elseif strcmp(param.val,'xb')
%         %partition index(SC)
%         ni = sum(result.data.f);                        %calculate fuzzy cardinality
%         si = sum(result.data.d.*result.data.f.^(m/2));  %calculate fuzzy variation
%         pii=si./ni;
%         mask = zeros(c,n,c);                            %calculate separation of clusters 
%         for i = 1:c
%             for j =1:c
%                 mask(j,:,i) = v(i,:);
%             end
%             dist(i) = sum(sum((mask(:,:,i) - v).^2));
%         end
%         s = dist;
% 
%         SC = sum(pii./s);
% 
%         %separation index (S)
%         S = sum(pii)./(N*min(dist));

        %Xie and Beni's index (XB)
%         XB = sum((sum(result.data.d.*result.data.f.^2))./(N*min(result.data.d)));
        [clusterNum,colNum]=size(result.cluster.v);
        
        minDist=10000000;
%         for t=1:clusterNum-1,
%             X0 = ones(clusterNum-t,1);
%             tempDiff = result.cluster.v(t+1:clusterNum,:)- X0*result.cluster.v(t,:);
%             tempMin=min(sum((tempDiff*eye(colNum).*tempDiff),2));
%             if  tempMin<minDist,
%                 minDist=tempMin;
%             end
%                
%         end
%         XB = sum(diag(transpose(result.data.f.^m)*result.data.d.^2))/(N*minDist);
        if c==1,
            XB=sum(sum(result.data.f.^m.*result.data.d.^2))/N;
        else
            
            for t=1:clusterNum-1,
                for j=t+1:clusterNum,
                    tempDiff=result.cluster.v(t,:)-result.cluster.v(j,:);
                    tempD=sum(tempDiff.^2);
                    minDist=min(minDist,tempD);
                end;
            end;
            XB = sum(sum(result.data.f.^m.*result.data.d.^2))/(N*minDist);
        end
        
        %results    
%     result.validity.SC = SC;
%     result.validity.S = S;
    result.validity.XB = XB;    
        
        
        
elseif strcmp(param.val,'DI')
        %Dunn's index (DI)
        %%for identification of compact and well separated clusters
        [m,label] = min(result.data.d');%crisp clustering(Kmeans)

        for i = 1:c
           index=find(label == i);
           dat{i}=data.X(index,:);
           meret(i)= size(dat{i},1);
        end
        mindistmatrix =ones(c,c)*inf;
        mindistmatrix2 =ones(c,c)*inf;
        
        for cntrCurrentClust = 1:c
           for cntrOtherClust = (cntrCurrentClust+1):c
               for cntrCurrentPoints = 1:meret(cntrCurrentClust)
                   dd = min(sqrt(sum([(repmat(dat{cntrCurrentClust}(cntrCurrentPoints,:),meret(cntrOtherClust),1) ...
                   -dat{cntrOtherClust}).^2]')));
                                       %calculate distances for an alternative Dunn index 
                   dd2 = min(abs(result.data.d(cntrCurrentClust,:)-result.data.d(cntrOtherClust,:)));
                       
                   if mindistmatrix(cntrCurrentClust,cntrOtherClust) > dd
                      mindistmatrix(cntrCurrentClust,cntrOtherClust) = dd;
                   end
                   if mindistmatrix2(cntrCurrentClust,cntrOtherClust) > dd2
                      mindistmatrix2(cntrCurrentClust,cntrOtherClust) = dd2;
                   end
               end
            end
        end

        minimalDist = min(min(mindistmatrix));
        minimalDist2 = min(min(mindistmatrix2));
        
        maxDispersion = 0;
        for cntrCurrentClust = 1:c
            actualDispersion = 0;
            for cntrCurrentPoints1 = 1:meret(cntrCurrentClust)
              dd = max(sqrt(sum([(repmat(dat{cntrCurrentClust}(cntrCurrentPoints1,:),meret(cntrCurrentClust),1) ...
                           -dat{cntrCurrentClust}).^2]')));
              if actualDispersion < dd
                 actualDispersion = dd;
              end
              if maxDispersion < actualDispersion
               maxDispersion = actualDispersion;
              end
           end
        end

        DI = minimalDist/maxDispersion;
        %alternative Dunn index
        ADI = minimalDist2/maxDispersion;
    %results
    result.validity.DI = DI;
    result.validity.ADI = ADI;   
elseif strcmp(param.val,'vos'),
    %vos index
    %computing the overlap measure
    omega=zeros(1,N);
    for i=1:N,
        temp=max(result.data.f(i,:));
        if temp>=0.8,
            omega(i)=0.1;
        elseif temp>=0.7 && temp<0.8,
            omega(i)=0.4;
        elseif temp>=0.6 && temp<0.7,
            omega(i)=0.7;
        else
            omega(i)=1.0;
        end;
    end;
    overlap=0.0;   %overlap measure;
    separation=10000000.0;       %separation measure;
    muSeq=unique(result.data.f(:));
    muSeqLen=size(muSeq,2);
    for j=1:c-1,
        for k=j+1:c,
            sep1=0.0;
            for t=1:N,
                g=min(result.data.f(t,j),result.data.f(t,k));
                sep1=max(sep1,g);
            end;
            separation=min(separation,sep1);
            for m=1:muSeqLen,
                tempMu=muSeq(m);
                for i=1:N,
                    if result.data.f(i,j)>=tempMu && result.data.f(i,k)>=tempMu,
                        overlap=overlap+omega(i)^2;
                    end;
                end;
            end;
        end;
%         separation=min(separation,sep1);
    end;
    %Computing the separation measure
    separation=1-separation;
    result.validity.vos.separation=separation;
    result.validity.vos.overlap=overlap*2/(c*(c-1));
elseif strcmp(param.val,'MPE-DMFP'),   %MPE-DMFP
    mu=result.data.f;
    MPE=mu.^2.*log2(mu);
    DM=0.0;
    Mk=zeros(N,c);  
    tempMean=zeros(N,1);
    for i=1:c,
        tempMean=tempMean+mu(:,i);
        Mk(:,i)=tempMean/N;
    end;
    for i=1:c-1,
        for j=i+1:c,
            diff=Mk(:,i)-Mk(:,j);
            DM=DM+2*sum(diff.^2);
        end;
    end;
    result.validity.MPE_DMFP=DM+sum(MPE(:));    
elseif strcmp(param.val,'aacc'),   %AACC
    [clusterNum,colNum]=size(result.cluster.v);
    minDist=10000000;
    clustCenter=result.cluster.v;
    punishFact.num=0.0;
    punishFact.Denom=1/clusterNum;
    for t=1:clusterNum-1,
        for j=t+1:clusterNum,
            tempDiff=result.cluster.v(t,:)-result.cluster.v(j,:);
            tempD=sum(tempDiff.^2);
            punishFact.num=punishFact.num+2*sqrt(tempD);
            minDist=min(minDist,tempD);
        end;
    end;
    
    punishFact.num=punishFact.num/(clusterNum*(clusterNum-1));    
    if c==1,
        AACC=sum(sum(result.data.f.^m.*result.data.d.^2))/(N+1);
    else
        AACC = (sum(sum(result.data.f.^m.*result.data.d.^2))+punishFact.num)/(N*minDist+punishFact.Denom);
    end;
        
        %results    
    
    result.validity.AACC = AACC;     
    
else
       %partition coefficient (PC)
        fm = (result.data.f).^m;
        PC = 1/N*sum(sum(fm));
        %classification entropy (CE)
        fm = (result.data.f).*log(result.data.f);
        CE = -1/N*sum(sum(fm)); 

     %results   
     result.validity.PC = PC;
     result.validity.CE = CE; 
end



