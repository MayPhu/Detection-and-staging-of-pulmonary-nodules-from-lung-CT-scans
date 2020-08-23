clear all;
close all;
clc;
%% Rewrite Dicom File sorting instance Number the lung ct scans from a single case 

srcFiles = dir('F:\Second Journal\NSCLC-Radiomics\NSCLC-Radiomics\LUNG1-003\01-01-2014-StudyID-34270\1-28595\*.dcm');  % the folder in which images exists
n = length(srcFiles);
% 
% X = zeros(512,512,n);
for i = 1:n
    fileName = ['F:\Second Journal\NSCLC-Radiomics\NSCLC-Radiomics\LUNG1-003\01-01-2014-StudyID-34270\1-28595\',srcFiles(i).name];
    info(i)=dicominfo(fileName);
    instanceNo(i)= info(i).InstanceNumber;
    I=dicomread(fileName);
    wfilename= ['F:\Second Journal\NSCLC-Radiomics\NSCLC-Radiomics\LUNG1-003\CT\',num2str(instanceNo(i)),'.dcm'];
    dicomwrite(I, wfilename);
%     seg= SegmentParenchyma(I); 
%      X(:,:,i) = seg; 
end
