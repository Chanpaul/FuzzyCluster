function [ output_args ] = ClustResDemo( dataset,clustRes)
%ClustResDemo Summary of this function goes here
%   demonstrate the clustering result with color map.
%******************Cluster result demonstration**************************

cmap=colormap;
[rowNum,colNum]=size(dataset);
NCLUST=size(clustRes.cluster.v,1);
tempmu=clustRes.data.f;
for i=0:NCLUST,
    nn=0;
    if i==0,
        color=[0 0 0];
        for j=1:rowNum,
            if max(tempmu(j,:))<0.6,
                nn=nn+1;
                A(nn,1)=dataset(j,1);
                A(nn,2)=dataset(j,2);
            end;
        end;
    else
        ic=int8((i*64.)/(NCLUST*1.));
        color=cmap(ic,:);
        for j=1:rowNum,
            if (tempmu(j,i)>0.6),
                nn=nn+1;
                A(nn,1)=dataset(j,1);
                A(nn,2)=dataset(j,2);
            end;
        end;
    end;
  hold on
  plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',3,'MarkerFaceColor',color,'MarkerEdgeColor',color);
%   plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',3,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
  
end;
hold off
end

