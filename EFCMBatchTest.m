%Test longitidal dataset, computing the validation index, output Sammon
%mapping graph;
clear all
close all
path(path,'FuzzyClusteringToolbox_m\FUZZCLUST')

dataDir='D:\UmassMed\Code\Dataset\';
%dataFileName='label_clust_dp_fig2_panelB_no_noise';  %2535 
%dataFileName='Afterimmu_mi10.txt';
dataFileName='beforeimmu_mi10.txt';
dataFile=strcat(dataDir,dataFileName);
%dataFile='iris.dat';
 
dataLabelFg=2; %Indicates which column is the label, 0 for unlabeled.
featureFg=[3:9];  %features used for clustering
otherFg=[1,10];  %others, e.g., the id and sth like.
%dataset=load(dataFile);
rowdata=importdata(dataFile,',');
dataset=rowdata.data;
[rowNum,colNum]=size(dataset);
m=2;
%validIndSt=[2 4 6]; % 2--xb, 4--vos, 6--AACC.
validIndSt={'xb','vos','aacc'};
clustName={'FCM','EFCM'};
simResDir='D:\UmassMed\Code\Results\LongitudinalTest\';
simResFile=strcat(simResDir,'simRes.xls');  %saving clustering error
clustResFile=strcat(simResDir,'clustRes.xls');  %saving fuzzy degrees.
s=0;
for k=1:size(clustName,2),
    for i=1:size(validIndSt,2),
        %validIndName=2;
        validIndName=validIndSt{i};
        clustAccuracy=[];
        clustError=[];
        s=s+1;
        validityIndex=[];
        g=1;
        sheetname=strcat(clustName{k},'-',validIndName,'-',dataFileName(1:5));
        %sheetname=sheetname{1};
        for j=1:10,
            tempDat=find(dataset(:,end)==j);
            res=EFCMTest(clustName{k},validIndName,dataset(tempDat,featureFg),m);
            validityIndex=[validityIndex;res.valRes];
            % %****************** Computing the accuracy of the clustering results******
            dataSz=size(tempDat,1);
            truth=dataset(tempDat,dataLabelFg);
            testResult=zeros(dataSz,1);
            for ss=1:dataSz,
                [va,ind]=max(res.data.f(ss,:));
                testResult(ss)=ind-1;
            end;
            [accuracy, true_labels, CM,corLabel] = calculateAccuracy(testResult, truth);
            
            tXlsRange=sprintf('A%d:C%d',g,g);
            xlswrite(clustResFile,rowdata.textdata([otherFg,dataLabelFg]),sheetname,tXlsRange);
            tXlsRange=sprintf('A%d:C%d',g+1,g+dataSz);
            xlswrite(clustResFile,dataset(tempDat,[otherFg,dataLabelFg]),sheetname,tXlsRange);
            tXlsRange=sprintf('D%d:E%d',g,g);
            xlswrite(clustResFile,corLabel,sheetname,tXlsRange);
            tXlsRange=sprintf('D%d:E%d',g+1,g+dataSz);
            xlswrite(clustResFile,res.data.f,sheetname,tXlsRange);
            g=g+dataSz+1;
            clustAccuracy=[clustAccuracy,accuracy];
            clustError=[clustError,1-accuracy];
            % %******************End of Computing the accuracy********************            
        end;
        sheetname1=dataFileName(1:5);
        xlsRange=sprintf('A%d:C%d',s,s);
        xlswrite(simResFile,{clustName{k},validIndName,mean(clustError)},sheetname1,xlsRange);
        figHandle=figure;
        %plot(2:size(validityIndex,2)+1,mean(validityIndex),'-o');
        plotdata=mean(validityIndex);
        plot(2:8,plotdata(2:8),'-o');
        xlabel('#cluster number');
        ylabel(validIndName);
        figName=strcat(simResDir,clustName(k),'-',validIndName,'-',dataFileName(1:5));
        figName=figName{1};
        saveas(figHandle,figName,'png');
        close(figHandle);
        %***************generate Sammon mapping graph**************
%         sammonMapName=strcat(figName,'-sammon');
%         sammon_mapping(dataset(tempDat,featureFg),dataset(tempDat,dataLabelFg),res,sammonMapName,m);
    end;
end;



