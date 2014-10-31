% clear all
% close all
path(path,'FuzzyClusteringToolbox_m\FUZZCLUST')
% validityIndex=[2,4]; %2-XB, 4-Vos
% validityNum=size(validityIndex,2);
% tic
dataDir='D:\UmassMed\Code\Dataset\';
% dataFileName='clust_dp_fig2_panelB.dat';
% dataFileName='clust_dp_fig2_panelC.dat';
% dataFileName='';
% dataFileName='label_clust_dp_fig2_panelB';
% dataFileName='label_clust_dp_fig2_panelC';
% dataFileName='label_clust_dp_fig2_panelC_no_noise';
dataFileName='label_clust_dp_fig2_panelB_no_noise';  %2535 
% dataFileName='label_clust_dp_fig2_panelD_no_noise';    %3119

% dataFileName='breast-cancer-wisconsin\breast-cancer-wisconsin.data';
% dataFileName='pima-indians-diabetes\pima-indians-diabetes.data';
% dataFileName='haberman.data';
% dataFileName='blood-transfusion.data';
% dataFileName='Aggregation.txt';
%the 3rd column of the data is the label.
% dataFileName='flame.txt';
% dataFileName='s3.txt';
% dataFileName='spiral.txt';
dataFile=strcat(dataDir,dataFileName);
% dataFile='iris.dat';
% 
dataLabelFlag=true; %the dataset is labeled at the last column.
% dataLabelFlag=false; %the dataset is unlabeled at the last column.
dataset=load(dataFile);
m=2;
res=EFCMTest('EFCM',2,dataset,dataLabelFlag,m);
% %************************Computing the running time**********************
% dataset=load(dataFile);
% clusType=unique(dataset(:,end));
% clustNum=zeros(1,size(clusType,2));
% clust=[];
% clust(size(clusType,1)).seq=[];
% for j=1:size(clusType,1),
%     clust(clusType(j)).seq=find(dataset(:,end)==clusType(j));
%     clustNum(clusType(j))=size(clust(clusType(j)).seq,1);
% end;
% clustRatio=clustNum./sum(clustNum);
% fDir='D:\UmassMed\Code\Results\';
% comResFile=strcat(fDir,'NewResult.txt');
% f=fopen(comResFile,'w');
% fprintf(f,'name****dataset******accuracy***running time****iterating times \n');
% elapseTime=zeros(1,1);
% m=2;
% validity=2;  % 2 for XB, 4 for VOS 5 for 
% sampleNumSet=[1500, 2000, 2500];
% %     dataset=load(dataFile);
% figure(11);
% plot(dataset(:,1),dataset(:,2),'o');
% loopNum=1;
% for i=1:size(sampleNumSet,2),
%     %****************generate dataset**********************************
%     sampDataset=[];
%     sampleNum=sampleNumSet(i);
%     for i=1:size(clusType,1),
%         nn= round(sampleNum*clustRatio(clusType(i)));
%         rng=('default');
%         rseq=randperm(clustNum(clusType(i)),nn);
%         dseq=clust(clusType(i)).seq(rseq);
%         sampDataset=[sampDataset;dataset(dseq,:)];
%         
%     end;
%     %**************end********************************  
%     
%     elapseTime=zeros(1,loopNum);
%     accuracy=zeros(1,loopNum);
%     for s=1:loopNum,
%         tic
%         res=EFCMTest('EFCM',validity,sampDataset,dataLabelFlag,m);
%         t=toc();
%         elapseTime(s)=t;
%         accuracy(s)=res.accuracy;        
%     end;
%     et=mean(elapseTime);
%     acc=mean(accuracy);
%     itertimes=res.iter;
%     fprintf(f,'CEF****%d******%f***%f*******%d \n',sampleNum,acc,et,itertimes);
%     elapseTime=zeros(1,loopNum);
%     accuracy=zeros(1,loopNum);
%     for s=1:loopNum,
%         tic
%         res=EFCMTest('FCM',validity,sampDataset,dataLabelFlag,m);
%         t=toc();
%         elapseTime(s)=t;
%         accuracy(s)=res.accuracy;        
%     end;
%     et=mean(elapseTime);
%     acc=mean(accuracy);
%     itertimes=res.iter;
%     fprintf(f,'FCM****%d******%f***%f*****%d \n',sampleNum,acc,et,itertimes);
%     
% end;
% 
% 
% %*************************End of computing running time*******************
%*********************COMputing validity values**************************


% %*******************For test II***************************************
% figDir='D:\UmassMed\Code\Results\20141002\';
% % comResFile=strcat(figDir,'pimaIndiaComResult.txt');
% % dataSetFlag='pimaIndia';
% % comResFile=strcat(figDir,'breastCancerComResult.txt');
% % dataSetFlag='breastCancer';
% comResFile=strcat(figDir,'irisComResult.txt');
% dataSetFlag='iris';
% f=fopen(comResFile,'w');
% m=2;
% validityIndex=[];  %2 for XB, 4 for VOS, 6 for ACAA;
% validityIndex(1).na='XB';
% validityIndex(1).id=2;
% validityIndex(2).na='VOS';
% validityIndex(2).id=4;
% validityIndex(3).na='ACAA';
% validityIndex(3).id=6;
% fcmName=[];
% fcmName(1).na='EFCM';
% fcmName(2).na='FCM';
% fprintf(f,'name****validity*****optimal cluster num.*****accuracy***running time \n');
% for i=1:size(fcmName,2),
%     for j=1:size(validityIndex,2),        
%         h=figure;
%         elapseTime=[];
%         for t=1:5,
%             tstart=tic;
%             res=EFCMTest(fcmName(i).na,validityIndex(j).id,dataset,dataLabelFlag,m);
%             telapsed = toc(tstart);
%             elapseTime(t)=telapsed;
%         end;
%         telapsed=mean(elapseTime);
%         plot(2:size(res.valRes,2),res.valRes(2:end),'ro-');
%         xlim([0,8]);
%         % ylim([0,1]);
%         grid;
%         xlabel('clustering number');
%         ylabel(strcat(fcmName(i).na,' on  ',dataSetFlag,' with ',validityIndex(j).na,' index'));
%         figName=strcat(dataSetFlag,'-',fcmName(i).na,'-',validityIndex(j).na);
%         saveas(h,strcat(figDir,figName),'png');
%         close(h);
%         
%         optiClustNum=size(res.cluster.v,1);
%         fprintf(f,'%s     %s    %d     %f    %f \n',fcmName(i).na,validityIndex(j).na,optiClustNum,res.accuracy,telapsed);
%     end;
% end;
% fclose(f);
% %**************************end of test II*********************************
% figure(2);
% res2=EFCMTest('FCM',4,dataset,dataLabelFlag,m);
% toc
%******************End of computing validity values*********************
%**************computing the cluster number with different m************

% clustNumSet=zeros(1,size(1.1:0.1:3,2));
% elapseTime=[];
% count=0;
% for m=1.1:0.1:3,   
%     count=count+1;
%     res=EFCMTest('EFCM',4,dataset,dataLabelFlag,m);
%     clustNum=size(res.cluster.v,1);
%     clustNumSet(count)=clustNum;
% end;
% %     subplot(2,3,j);
% plot(1.1:0.1:3,clustNumSet,'o-');
% xlabel('parameter m');
% ylabel('cluster num');
% ylim([0 5]);
%********************end of computing********************************

%**********************************************
% elapseTime=[];
% count=0;
% j=0;
% for i=4000:4000,
%     j=j+1;
%     count=0;
%     clustNumSet=[];
%     for m=1:0.1:10,   
%         count=count+1;
%         res=EFCMTest('EFCM',2,dataset(1:i,:),dataLabelFlag,m);
%         clustNum=size(res.cluster.v,1);
%         clustNumSet(count)=clustNum;
%     end;
% %     subplot(2,3,j);
%     plot(1:0.1:10,clustNumSet,'o-');
%     xlabel('parameter m');
%     ylabel('cluster num');
% %     figure(1);
% %     subplot(4,3,count);
% %     plot(dataset(1:i,1),dataset(1:i,2),'o','MarkerSize',2);
% %     titleStr=sprintf('Original %d points',i);
% %     title(titleStr);
% % 
% %     tstart=tic;
% %     res=EFCMTest('FCM',4,dataset(1:i,:),dataLabelFlag,2);
% %     telapsed = toc(tstart);
% %     figure(2);
% %     subplot(4,3,count);
% %     ClustResDemo(dataset(1:i,:),res);
% %     titleStr=sprintf('clustering %d points',i);
% %     title(titleStr);
% %     elapseTime(count)=telapsed;
%     
%     
% end;
% % figure(3);
% % plot(2000:200:4000,elapseTime,'o-');
% % xlabel('#Number of points');
% % ylabel('Computing time (s)');
% % toc
% %************************************************************************
% 
% % for i=1:2,
% %     res=EFCMTest('EFCM',4,dataset,dataLabelFlag,2);
% %     val=res.valRes;
% %     n=size(val,2);
% %     if i==1,
% %         t='XB';
% %     else
% %         t='Vos';
% %     end;
% %     subplot(2,2,2*i-1);    
% %     plot(2:n,val(2:n),'o-');
% %     xlabel('#num of cluster');
% %     ylabel(t);
% %     title(strcat('EFCM with ',t));
% %     clear res;
% %     clear val;
% %     
% %     res=EFCMTest('FCM',validityIndex(1));
% %     val=res.valRes;
% %     n=size(val,2);
% %     subplot(2,2,2*i); 
% % %     figure(i);
% %     plot(2:n,val(2:n),'o-');
% %     xlabel('#num of cluster');
% %     ylabel(t);
% %     title(strcat('FCM with ',t));
% % 
% %     clear res;
% %     clear val;
% % end;
