clear all
close all
path(path,'FuzzyClusteringToolbox_m\FUZZCLUST')

dataDir='D:\UmassMed\Code\Dataset\';
%dataFileName='label_clust_dp_fig2_panelB_no_noise';  %2535 
dataFileName='waveform.data';
dataFile=strcat(dataDir,dataFileName);
%dataFile='iris.dat';
 
dataLabelFlag=true; %the dataset is labeled at the last column.
% dataLabelFlag=false; %the dataset is unlabeled at the last column.
dataset=load(dataFile);
m=2;
res=EFCMTest('EFCM',2,dataset,dataLabelFlag,m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%gggggggggggggggggggggggggggggggggggggggggg

