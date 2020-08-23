clc;
clear all;
close all;
%% load trained GA-BPNN 
load ('F:\Experiment (Coding)\RandomSubSamplingValidation\NN GA\YeastData.mat'); 
% clearvars -except BPnet

% load ('F:\Experiment (Coding)\RandomSubSamplingValidation\YeastData.mat');
Data = YeastData;
populationsize = size(Data,1);


for i=1:100
%% Random sampling (50/50,60/40,70/30,80/20,90/10)
s1 = round(populationsize* 90/100);
s2= populationsize-s1;

TrainingData  = datasample(Data,s1);
TestingData  = datasample(Data,s2);


%Training
TrainInput = TrainingData(:,1:8);
TrainLabels = TrainingData(:,9);
TrainClasses = full(ind2vec(TrainLabels'));

%Testing
TestingInput =TestingData(:,1:8);
TestLabels = TestingData(:,9);
TestClasses = full(ind2vec(TestLabels'));

TestingInput =TestingInput';
TestLabels= TestLabels';

%% improved BPNN
tic 
outputs = BPnet(TestingInput);
Predictions = vec2ind(outputs);
% [c,cm,ind,per] = confusion(TestClasses,outputs);
BPConfusionmax = confusionmat(TestLabels,Predictions);

[BPResult,BPRefereceResult]=confusion.getValues(BPConfusionmax);
BPtime= toc;

% ROC 
[BPX,BPY,BPT,BPAUC,BPOPTROCPT] = perfcurve(TestLabels,Predictions,'2');
% % plot(BPX,BPY,'-y','LineWiBPh',1.5);
% % hold on 

BP{i}= [BPRefereceResult.AccuracyInTotal;BPRefereceResult.Sensitivity;BPRefereceResult.Specificity;BPRefereceResult.Precision;BPRefereceResult.F1_score;BPAUC;BPtime];
BPXX{i}= BPX;
BPYY{i} =BPY;

%Accuracy
BPAcc(i,:)  =  BPResult.Accuracy; 
BPAcc(isnan(BPAcc)) = 0;

%Sensitivity
BPSen(i,:) =  BPResult.Sensitivity; 
BPSen(isnan(BPSen)) = 0;

%Specificity
BPSpe(i,:) =  BPResult.Specificity;
BPSpe(isnan(BPSpe)) = 0;


%Precision
BPPre(i,:)  =  BPResult.Precision; 
BPPre(isnan(BPPre)) = 0;

%F1_score
BPF1(i,:)  = BPResult.F1_score; 
BPF1(isnan(BPF1)) = 0;

%AUC
BPAUCTotal(i,:)  =BPAUC; 
BPAUCTotal(isnan(BPAUCTotal)) = 0;

%Time
BPTIME(i,:) = BPtime; 



end
% mean and std of each classifier 
finalmatrix  = [BPAcc,BPSen,BPSpe,BPPre,BPF1,BPAUCTotal,BPTIME];
MeanStdBP = [mean(BPAcc), mean(BPSen),mean(BPSpe),mean(BPPre),mean(BPF1),mean(BPAUCTotal),mean(BPTIME);
             std(BPAcc), std(BPSen),std(BPSpe),std(BPPre),std(BPF1),std(BPAUCTotal),std(BPTIME)];
% 
% 
%          
% MeanStd = [MeanStdBP';MeanStdKNN';MeanStdSVM';MeanStBPB';MeanStdBPNN'];