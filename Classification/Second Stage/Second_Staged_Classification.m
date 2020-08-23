clc;
clear all;
close all;
%% 
load ('F:\Experiment (Coding)\Hybrid Feature Selection and Classification Comparison1\WithoutFeatureSelection.mat');
Data = WithoutFeatureSelection;
populationsize = size(Data,1);

for i=1 :100
% Random sampling (50/50,60/40,70/30,80/20,90/10)
s1 = round(populationsize* 70/100);
s2= populationsize-s1;
Data  = datasample(Data,populationsize);
TrainingData  = datasample(Data,s1);
TestingData  = datasample(Data,s2);



%Training & Testing
TrainInput = TrainingData(:,2:32);
TrainLabels = TrainingData(:,1);

TestingInput =TestingData(:,2:32); %*****
TestLabels = TestingData(:,1);


%% Decision Tree
tic 
[trainedDT, DTAccuracy,DTPredictions,DTMScores] = trainDecisionTreeSecondStage(TrainingData);
[DTTestPredictions,DTMTestScores] = trainedDT.predictFcn(TestingInput);
Confusionmax_DT = confusionmat(TestLabels,DTTestPredictions);
[DTResult,DTRefereceResult]=confusion.getValues(Confusionmax_DT);
timeDT= toc;
% DTAUC = 2* ((DTResult.Sensitivity*DTResult.Precision)/(DTResult.Sensitivity + DTResult.Precision));
[DTX,DTY,DTT,DTAUC,DTOPTROCPT] = perfcurve(TestLabels,DTMTestScores(:,2),'1');
DT{i} = [DTResult.Accuracy;DTResult.Sensitivity;DTResult.Specificity;DTResult.Precision;DTResult.F1_score;DTAUC;timeDT];
DTXX{i}= DTX;
DTYY{i} =DTY;

%% KNN
tic
[trainedKNN, KNNAccuracy,KNNPredictions,KNNScores] = trainKNNSecondStage(TrainingData);
[KNNTestPredictions,KNNTestScores] = trainedKNN.predictFcn(TestingInput);
Confusionmax_KNN = confusionmat(TestLabels,KNNTestPredictions);
[KNNResult,KNNRefereceResult]=confusion.getValues(Confusionmax_KNN);
timeKNN= toc;
% KNNAUC = 2* ((KNNResult.Sensitivity*KNNResult.Precision)/(KNNResult.Sensitivity + KNNResult.Precision));
[KNNX,KNNY,KNNT,KNNAUC,KNNOPTROCPT] = perfcurve(TestLabels,KNNTestScores(:,2),'1');
KNN{i} = [KNNResult.Accuracy;KNNResult.Sensitivity;KNNResult.Specificity;KNNResult.Precision;KNNResult.F1_score;KNNAUC;timeKNN];
KNNXX{i}= KNNX;
KNNYY{i} =KNNY;

%% SVM
tic
[trainedSVM, SVMAccuracy,SVMPredictions,SVMScores] = trainSVMSecondStage(TrainingData);
[SVMTestPredictions,SVMTestScores] = trainedSVM.predictFcn(TestingInput);
Confusionmax_SVM = confusionmat(TestLabels,SVMTestPredictions);
[SVMResult,SVMRefereceResult]=confusion.getValues(Confusionmax_SVM);
timeSVM= toc;
% SVMAUC =2* ((SVMResult.Sensitivity*SVMResult.Precision)/(SVMResult.Sensitivity + SVMResult.Precision));
[SVMX,SVMY,SVMT,SVMAUC,SVMOPTROCPT] = perfcurve(TestLabels,SVMTestScores(:,2),'1');
SVM{i} = [SVMResult.Accuracy;SVMResult.Sensitivity;SVMResult.Specificity;SVMResult.Precision;SVMResult.F1_score;SVMAUC;timeSVM];
SVMXX{i} = SVMX;
SVMYY{i} = SVMY;

%% Tree Bagger 
tic
[trainedEnsembleTrees, EnsembleTreesAccuracy,EnsembleTreesPredictions,EnsembleTreesScores] = trainEnsembleTreesSecondStage(TrainingData);
[ETTestPredictions,EnsembleTestTreesScores] = trainedEnsembleTrees.predictFcn(TestingInput);
Confusionmax_EnsembleTrees = confusionmat(TestLabels,ETTestPredictions);
[TBResult,TBRefereceResult]=confusion.getValues(Confusionmax_EnsembleTrees);
timeTB= toc;
% TBAUC =2* ((TBResult.Sensitivity*TBResult.Precision)/(TBResult.Sensitivity + TBResult.Precision));
[TBX,TBY,TBT,TBAUC,TBOPTROCPT] = perfcurve(TestLabels,EnsembleTestTreesScores(:,2),'1');
TB{i}= [TBResult.Accuracy;TBResult.Sensitivity;TBResult.Specificity;TBResult.Precision;TBResult.F1_score;TBAUC;timeTB];
TBXX {i} = TBX;
TBYY{i}=  TBY;

%% RF
tic
[trainedRF, RFAccuracy,RFPredictions,RFScores] = trainRFSecondStage(TrainingData);
[RFTestPredictions,RFTestScores] = trainedRF.predictFcn(TestingInput);
Confusionmax_RF = confusionmat(TestLabels,RFTestPredictions);
[RFResult,RFRefereceResult]=confusion.getValues(Confusionmax_RF);
timeRF= toc;
% RFAUC =2* ((RFResult.Sensitivity*RFResult.Precision)/(RFResult.Sensitivity + RFResult.Precision));
[RFX,RFY,RFT,RFAUC,RFOPTROCPT] = perfcurve(TestLabels,RFTestScores(:,2),'1');
RF{i}= [RFResult.Accuracy;RFResult.Sensitivity;RFResult.Specificity;RFResult.Precision;RFResult.F1_score;RFAUC;timeRF];
RFXX{i} = RFX;
RFYY{i}=  RFY;


%Accuracy
DTAcc(i,:)  =  DTResult.Accuracy; 
KNNAcc(i,:) = KNNResult.Accuracy;
SVMAcc(i,:) = SVMResult.Accuracy;
TBAcc(i,:) = TBResult.Accuracy;
RFAcc(i,:)= RFResult.Accuracy;
 
%Sensitivity

DTSen(i,:) =  DTResult.Sensitivity; 
KNNSen(i,:) = KNNResult.Sensitivity;
SVMSen(i,:) = SVMResult.Sensitivity;
TBSen(i,:) = TBResult.Sensitivity;
RFSen(i,:)= RFResult.Sensitivity;
DTSen(isnan(DTSen)) = 0;
KNNSen(isnan(KNNSen)) = 0;
SVMSen(isnan(SVMSen)) = 0;
TBSen(isnan(TBSen)) = 0;
RFSen(isnan(RFSen)) = 0;

 
%Specificity
DTSpe(i,:) =  DTResult.Specificity; 
KNNSpe(i,:) = KNNResult.Specificity;
SVMSpe(i,:) = SVMResult.Specificity;
TBSpe(i,:) = TBResult.Specificity;
RFSpe(i,:)= RFResult.Specificity;
 
%Precision
DTPre(i,:)  =  DTResult.Precision; 
KNNPre(i,:) = KNNResult.Precision;
SVMPre(i,:) = SVMResult.Precision;
TBPre(i,:) = TBResult.Precision;
RFPre(i,:)= RFResult.Precision;
 
%F1_score
DTF1(i,:)  = DTResult.F1_score; 
KNNF1(i,:)= KNNResult.F1_score;
SVMF1(i,:) = SVMResult.F1_score;
TBF1(i,:) = TBResult.F1_score;
RFF1(i,:)= RFResult.F1_score;
DTF1(isnan(DTF1)) = 0;
KNNF1(isnan(KNNF1)) = 0;
SVMF1(isnan(SVMF1)) = 0;
TBF1(isnan(TBF1)) = 0;
RFF1(isnan(RFF1)) = 0;


%AUC
DTAUCTotal(i,:)  =DTAUC; 
KNNAUCTotal(i,:) = KNNAUC;
SVMAUCTotal(i,:) =SVMAUC;
TBAUCTotal(i,:) = TBAUC;
RFAUCTotal(i,:)= RFAUC;
DTAUCTotal(isnan(DTAUCTotal)) = 0;
KNNAUCTotal(isnan(KNNAUCTotal)) = 0;
SVMAUCTotal(isnan(SVMAUCTotal)) = 0;
TBAUCTotal(isnan(TBAUCTotal)) = 0;
RFAUCTotal(isnan(RFAUCTotal)) = 0;


%Time
DTTIME(i,:) = timeDT; 
KNNTIME(i,:) =timeKNN;
SVMTIME(i,:) =timeSVM;
TBTIME(i,:) =timeTB;
RFTIME(i,:)=timeRF;


end

% mean and std of each classifier 
MeanStdDT = [mean(DTAcc), mean(DTSen),mean(DTSpe),mean(DTPre),mean(DTF1),mean(DTAUCTotal),mean(DTTIME);
             std(DTAcc), std(DTSen),std(DTSpe),std(DTPre),std(DTF1),std(DTAUCTotal),std(DTTIME)];

MeanStdKNN = [mean(KNNAcc), mean(KNNSen),mean(KNNSpe),mean(KNNPre),mean(KNNF1),mean(KNNAUCTotal),mean(KNNTIME);
             std(KNNAcc), std(KNNSen),std(KNNSpe),std(KNNPre),std(KNNF1),std(KNNAUCTotal),std(KNNTIME)];
         
MeanStdSVM = [mean(SVMAcc), mean(SVMSen),mean(SVMSpe),mean(SVMPre),mean(SVMF1),mean(SVMAUCTotal),mean(SVMTIME);
             std(SVMAcc), std(SVMSen),std(SVMSpe),std(SVMPre),std(SVMF1),std(SVMAUCTotal),std(SVMTIME)];
         
MeanStdTB = [mean(TBAcc), mean(TBSen),mean(TBSpe),mean(TBPre),mean(TBF1),mean(TBAUCTotal),mean(TBTIME);
             std(TBAcc), std(TBSen),std(TBSpe),std(TBPre),std(TBF1),std(TBAUCTotal),std(TBTIME)];
         
MeanStdRF = [mean(RFAcc), mean(RFSen),mean(RFSpe),mean(RFPre),mean(RFF1),mean(RFAUCTotal),mean(RFTIME);
             std(RFAcc), std(RFSen),std(RFSpe),std(RFPre),std(RFF1),std(RFAUCTotal),std(RFTIME)];
         
         
MeanStd = [MeanStdDT';MeanStdKNN';MeanStdSVM';MeanStdTB';MeanStdRF'];