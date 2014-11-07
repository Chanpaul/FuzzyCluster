clear all
close all
path(path,'FuzzyClusteringToolbox_m\FUZZCLUST')

dataDir='D:\UmassMed\Code\Dataset\';
%dataFileName='label_clust_dp_fig2_panelB_no_noise';  %2535 
dataFileName='waveform.data';
%dataFileName='waveform1.data';
dataFile=strcat(dataDir,dataFileName);
%dataFile='iris.dat';
 
dataLabelFlag=true; %the dataset is labeled at the last column.
% dataLabelFlag=false; %the dataset is unlabeled at the last column.
dataset=load(dataFile);
[rowNum,colNum]=size(dataset);
m=2;
res=EFCMTest('EFCM',4,dataset,dataLabelFlag,m);
sprintf('the accuracy of EFCM is %f',res.accuracy);
figure(1);
semilogy(2:size(res.valRes,2),res.valRes(2:end),'o');
% plot(2:size(res.valRes,2)+1,res.valRes,'o');
res1=EFCMTest('FCM',4,dataset,dataLabelFlag,m);
sprintf('the accuracy of FCM is %f',res1.accuracy);
figure(2);
% plot(2:size(res1.valRes,2)+1,res1.valRes,'*');
semilogy(2:size(res1.valRes,2),res1.valRes(2:end),'*');

