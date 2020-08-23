clc;
clear all;
close all;
%%
load ('F:\Experiment (Coding)\Hybrid Feature Selection and Classification Comparison1\LIDC-IDRI.mat');
Data = D;
populationsize = size(Data,1);
for i=1:100
%% Random sampling (50/50,60/40,70/30,80/20,90/10)
s1 = round(populationsize* 70/100);
s2= populationsize-s1;

TrainingData  = datasample(Data,s1);
TestingData  = datasample(Data,s2);


%Training
TrainInput = TrainingData(:,1:14);
TrainLabels = TrainingData(:,15);

%Testing
TestingInput =TestingData(:,1:14);
TestLabels = TestingData(:,15);


%% Decision Tree
tic 
[trainedDT, DTAccuracy,DTPredictions,DTMScores] = trainDecisionTreeFirstStage(TrainingData);
[DTTestPredictions,DTMTestScores] = trainedDT.predictFcn(TestingInput);
Confusionmax_DT = confusionmat(TestLabels,DTTestPredictions);
[DTResult,DTRefereceResult]=confusion.getValues(Confusionmax_DT);
timeDT= toc;

% ROC (Decision Tree)
[DTX,DTY,DTT,DTAUC,DTOPTROCPT] = perfcurve(TestLabels,DTMTestScores(:,2),'1');
% plot(DTX,DTY,'-y','LineWidth',1.5);
% hold on 
DT{i}= [DTRefereceResult.AccuracyInTotal;DTRefereceResult.Sensitivity;DTRefereceResult.Specificity;DTRefereceResult.Precision;DTRefereceResult.F1_score;DTAUC;timeDT];
DTXX{i}= DTX;
DTYY{i} =DTY;
%% KNN
tic
[trainedKNN, KNNAccuracy,KNNPredictions,KNNScores] = trainKNNFirstStage(TrainingData);
[KNNTestPredictions,KNNTestScores] = trainedKNN.predictFcn(TestingInput);
Confusionmax_KNN = confusionmat(TestLabels,KNNTestPredictions);
[KNNResult,KNNRefereceResult]=confusion.getValues(Confusionmax_KNN);
timeKNN= toc;

% ROC (KNN)
[KNNX,KNNY,KNNT,KNNAUC,KNNOPTROCPT] = perfcurve(TestLabels,KNNTestScores(:,2),'1');
% plot(KNNX,KNNY,'-g','LineWidth',1.5);
% hold on 

KNN = [KNNRefereceResult.AccuracyInTotal;KNNRefereceResult.Sensitivity;KNNRefereceResult.Specificity;KNNRefereceResult.Precision;KNNRefereceResult.F1_score;KNNAUC;timeKNN];
KNNXX{i}= KNNX;
KNNYY{i} =KNNY;

%% SVM
tic
[trainedSVM, SVMAccuracy,SVMPredictions,SVMScores] = trainSVMFirstStage(TrainingData);
[SVMTestPredictions,SVMTestScores] = trainedSVM.predictFcn(TestingInput);
Confusionmax_SVM = confusionmat(TestLabels,SVMTestPredictions);
[SVMResult,SVMRefereceResult]=confusion.getValues(Confusionmax_SVM);
timeSVM= toc;

% ROC (SVM)
[SVMX,SVMY,SVMT,SVMAUC,SVMOPTROCPT] = perfcurve(TestLabels,SVMTestScores(:,2),'1');
% plot(SVMX,SVMY,'-b','LineWidth',1.5);
% hold on 

SVM {i} = [SVMRefereceResult.AccuracyInTotal;SVMRefereceResult.Sensitivity;SVMRefereceResult.Specificity;SVMRefereceResult.Precision;SVMRefereceResult.F1_score;SVMAUC;timeSVM];
SVMXX{i} = SVMX;
SVMYY{i} = SVMY;

%% Tree Bagger 
tic
[trainedEnsembleTrees, EnsembleTreesAccuracy,EnsembleTreesPredictions,EnsembleTreesScores] = trainEnsembleTreesFirstStage(TrainingData);
[ETTestPredictions,EnsembleTestTreesScores] = trainedEnsembleTrees.predictFcn(TestingInput);
Confusionmax_EnsembleTrees = confusionmat(TestLabels,ETTestPredictions);
[ETResult,TBRefereceResult]=confusion.getValues(Confusionmax_EnsembleTrees);
timeTB= toc;

% ROC (TB)
[TBX,TBY,TBT,TBAUC,TBOPTROCPT] = perfcurve(TestLabels,EnsembleTestTreesScores(:,2),'1');
% plot(TBX,TBY,'-c','LineWidth',1.5);
% hold on 

TB {i} = [TBRefereceResult.AccuracyInTotal;TBRefereceResult.Sensitivity;TBRefereceResult.Specificity;TBRefereceResult.Precision;TBRefereceResult.F1_score;TBAUC;timeTB];
TBXX {i} = TBX;
TBYY{i}=  TBY;

%% Random Forest
tic
[trainedRF, RFAccuracy,RFPredictions,RFScores] = trainDecisionRFFirstStage(TrainingData);
[RFTestPredictions,RFTestScores] = trainedRF.predictFcn(TestingInput);
Confusionmax_RF = confusionmat(TestLabels,RFTestPredictions);
[RFResult,RFRefereceResult]=confusion.getValues(Confusionmax_RF);
timeRF= toc;

% ROC (RF)
[RFX,RFY,RFT,RFAUC,RFOPTROCPT] = perfcurve(TestLabels,RFTestScores(:,2),'1');
% plot(RFX,RFY,'-r','LineWidth',1.5);
% hold on 

RF{i}= [RFRefereceResult.AccuracyInTotal;RFRefereceResult.Sensitivity;RFRefereceResult.Specificity;RFRefereceResult.Precision;RFRefereceResult.F1_score;RFAUC;timeRF];
RFXX{i} = RFX;
RFYY{i}=  RFY;
%
%Accuracy
DTAcc(i,:)  =  DTRefereceResult.AccuracyInTotal; 
KNNAcc(i,:) = KNNRefereceResult.AccuracyInTotal;
SVMAcc(i,:) = SVMRefereceResult.AccuracyInTotal;
TBAcc(i,:) = TBRefereceResult.AccuracyInTotal;
RFAcc(i,:)= RFRefereceResult.AccuracyInTotal;

%Sensitivity
DTSen(i,:) =  DTRefereceResult.Sensitivity; 
KNNSen(i,:) = KNNRefereceResult.Sensitivity;
SVMSen(i,:) = SVMRefereceResult.Sensitivity;
TBSen(i,:) = TBRefereceResult.Sensitivity;
RFSen(i,:)= RFRefereceResult.Sensitivity;

%Specificity
DTSpe(i,:) =  DTRefereceResult.Specificity; 
KNNSpe(i,:) = KNNRefereceResult.Specificity;
SVMSpe(i,:) = SVMRefereceResult.Specificity;
TBSpe(i,:) = TBRefereceResult.Specificity;
RFSpe(i,:)= RFRefereceResult.Specificity;

%Precision
DTPre(i,:)  =  DTRefereceResult.Precision; 
KNNPre(i,:) = KNNRefereceResult.Precision;
SVMPre(i,:) = SVMRefereceResult.Precision;
TBPre(i,:) = TBRefereceResult.Precision;
RFPre(i,:)= RFRefereceResult.Precision;

%F1_score
DTF1(i,:)  = DTRefereceResult.F1_score; 
KNNF1(i,:)= KNNRefereceResult.F1_score;
SVMF1(i,:) = SVMRefereceResult.F1_score;
TBF1(i,:) = TBRefereceResult.F1_score;
RFF1(i,:)= RFRefereceResult.F1_score;

%AUC
DTAUCTotal(i,:)  =DTAUC; 
KNNAUCTotal(i,:) = KNNAUC;
SVMAUCTotal(i,:) =SVMAUC;
TBAUCTotal(i,:) = TBAUC;
RFAUCTotal(i,:)= RFAUC;

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

