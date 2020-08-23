clc;
clear all;
close all;
tic

%% Getting DICOM information of the current input CT exam

% eg-LIDC-IDRI-0001 is a folder of CT slcies for one patient, all silces in the folder has the same dicom info
% so read one CT slice (000001.dcm) to get dicom information  
info= dicominfo('F:\Experiment\LIDC-IDRI-0001\000001.dcm'); 

% m is the rescaleslope and b is the rescaleintercept
% they are used to assign threshold value for segmentation

m= info.RescaleSlope;
b= info.RescaleIntercept;


%%  Performing segmentation slice by slice and generate 3D data

% All slices in the input CT exam are rewirted (sequentially by instance number) into a folder named "Write" 
% using WritebyInstanceNumber.m

srcFiles = dir('F:\Experiment\LIDC-IDRI-0001\Write\*.dcm');  % the folder in which images exists
n = length(srcFiles);

for q =1:n
     fileName = ['F:\Experiment\LIDC-IDRI-0001\Write\',num2str(q),'.dcm'];
     I=dicomread(fileName);    
     VolumeCT(:,:,q) = I;   % to show stacked I CT scans)
          
     %% segmentation 
     [body,lung,trachea,vertebra,lungmask,lungmaskedge,finalcandidate] = Segmentation(I);
     VolumeBody (:,:,q)= body;
     VolumeLung(:,:,q) = lung;
     RoughTrachea(:,:,q) = trachea;
     VolumeVertebra(:,:,q) = vertebra;
     VolumeLungMask(:,:,q) = lungmask;
     VolumeLungMaskEdge(:,:,q) = lungmaskedge;
     Volumefinalcandidate(:,:,q) = finalcandidate;
%    Volumefissures(:,:,q) = fissures;
%    Volumevessels(:,:,q) = vessels;   
end
VolumeLungMask= flipdim(VolumeLungMask,3); 
figure,isosurface(VolumeLungMask,0.3),title('Volume Lung Mask');
camlight; lighting gouraud  
% set(gca,'XColor', 'none','YColor','none','ZColor','none')

%    figure,isosurface(Volumefinalcandidate,0.3),title('Volume Lung Mask Edge');

Volumefinalcandidate= flipdim(Volumefinalcandidate,3); 
figure,isosurface(Volumefinalcandidate,1),title('Volume Candidates');
camlight; lighting gouraud  
% set(gca,'XColor', 'none','YColor','none','ZColor','none')

%% Nodule detection 

VolumeNodule= squeeze(num2cell(permute(Nodule,[1,2,4,3]),1:3));
VolumeNodule= VolumeNodule{1};
labelNodule = bwconncomp(VolumeNodule, 26);
lNodule= labelmatrix(labelNodule);
statusNodule = regionprops(Nodule,'Area');      
NoduleArea= extractfield(statusNodule,'Area');
[biggestNoduleArea,indexNodule]= max(NoduleArea); 
FinalNodule  = (lNodule==indexNodule);

%% trachea segmentation   
 statusTrachea = regionprops(RoughTrachea,'Area');      
 AreaTrachea= extractfield(statusTrachea,'Area');
 [biggestAreaTrachea,indexTrachea]= max(AreaTrachea);      
 labelTrachea = bwconncomp(RoughTrachea,26);
 ltrachea = labelmatrix(labelTrachea);      
 VolumeTrachea =(ltrachea==indexTrachea);
 
 %% bronchi (downcarina)  and trachea (uptrachea)  segmentation
 VolumeTrachea = padarray(VolumeTrachea,[3 3],0,'both');
 skel = Skeleton3D(VolumeTrachea);
 BP = branchpoints3(skel);
 labelBP = bwconncomp(BP, 26);
 lbp = labelmatrix(labelBP); 
 point= (lbp==1);
 statuspoint = regionprops(point,'Centroid');  
 CentroidBP= extractfield(statuspoint,'Centroid');
 x= CentroidBP(1,1);
 y= CentroidBP(1,2);
 z= CentroidBP(1,3);
 upCarina = VolumeTrachea(:,:,1:z);
 downCarina  = VolumeTrachea(:,:,z:end); 

%% right and left segmentation ( from the centriod point of the Carina) 
 labelLung = bwconncomp(VolumeLung, 26);
 lLung= labelmatrix(labelLung);
 statusLung = regionprops(lLung,'Area');  
 AreaLung= extractfield(statusLung,'Area');
 [biggestAreaLung,indexLung]= sort(AreaLung,'descend');
 indexes  = indexLung(1:2);
 VolumeLung1 =  (lLung==indexes(1));
 VolumeLung2 =  (lLung==indexes(2));
 VolumeLungClear = imadd ( VolumeLung1,VolumeLung2);
 % segment left and right by centriod of carnia 
 % [A(1:B(1),1) ; zeros(5-B(1),1)];
  
 
 y=round(y);
 leftLung = VolumeLungClear; 
 leftLung(:,y:end,:) = 0;

rightLung= VolumeLungClear;
rightLung(:,1:y,:) = 0;


 %% rib segmentaiton        
  labelVetrebra = bwconncomp(VolumeVertebra, 26);
  lVertebra = labelmatrix(labelVetrebra);  
  statusVertebra = regionprops(lVertebra,'Area');  
  AreaVertebra= extractfield(statusVertebra,'Area');
  [biggestArea,index]= max(AreaVertebra); 
  VolumeRib =  (lVertebra==index);
  
  %% check nodule is inside right lung 
  labelrightLung = bwconncomp(rightLung, 26);
  lrightLung= labelmatrix(labelrightLung); 
  statusrightLung = regionprops(lrightLung,'Area');  
  ArearightLung= extractfield(statusrightLung,'Area');
  [biggestArearigthLung,indexrightLung]= max(ArearightLung); 
  rightLungtoCheck =  (lrightLung==indexrightLung);
  ResultInRightLung  = CheckInsideOrOutside(rightLungtoCheck,FinalNodule);
  
  %% check nodule is inside left lung 
  labelleftLung = bwconncomp(leftLung, 26);
  lleftLung= labelmatrix(labelleftLung); 
  statusleftLung = regionprops(lleftLung,'Area');  
  ArealleftLung= extractfield(statusleftLung,'Area');
  [biggestArealleftLung,indexlleftLung]= max(ArealleftLung); 
  leftLungtoCheck =  (lleftLung==indexlleftLung);
  ResultInLeftLung  = CheckInsideOrOutside(leftLungtoCheck,FinalNodule);
  
%% check nodule is inside or outside of UpCarina( Trachea) 
labelupCarina = bwconncomp(upCarina, 26);
lupCarina= labelmatrix(labelupCarina); 
ResultInvadeTrachea  = CheckInsideOrOutside(lupCarina,FinalNodule);

%% check nodule is inside or outside of DownCarina(Bronchi,Hilar Regions) 
labeldownCarina = bwconncomp(downCarina, 26);
ldownCarina= labelmatrix(labeldownCarina); 
ResultInvadeBronchus  = CheckInsideOrOutside(ldownCarina,FinalNodule);

%%  check nodule is inside or outside of rib
ResultInvadeRibs  = CheckInsideOrOutside(VolumeRib,FinalNodule);
  
%% check nodule is pleural attached or not 

ResultPleuralAttach = CheckPleuralAttachment(rightLungtoCheck,leftLungtoCheck,FinalNodule);
  
%% final localization features 

% LungLocation ( determine RightLung(1) or LeftLung(0))
if ResultInRightLung == 1 
    LungLocation = 1;
elseif ResultInLeftLung==1
    LungLocation =2;
else 
    LungLocation =0;
end

% PleuralAttach ( determine PleuralAttach(1) or NonPleuralAttach(0))

if ResultPleuralAttach== 1
    Pleural=1;
else 
    Pleural=0;
end

% Zone ( Zone1(Bronchus,Hilar) Zone2( Rib) Zone3(Trachea,Carnia,Epsophugus)
  
if ResultInvadeBronchus==1
    Zone=1;
elseif  ResultInvadeTrachea==1
    Zone=2;
elseif  ResultInvadeRibs==1
    Zone=3;
else 
    Zone=0;
end
LocalizationFeatures = [ LungLocation Pleural Zone];
Features = featureExtraction(FinalNodule,m,b); 

FinalFeatures = [LocalizationFeatures Features];
save(saveResultfilename);
timeElapsed= toc;


%% 3D visualinzation 
% figure
% surf1=isosurface(rightLungtoCheck,0);
% p1 = patch(surf1);
% isonormals(rightLung,p1);
% set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
% daspect([1,1,1])
% view(3); axis tight
% camlight; lighting gouraud
% hold on
% 
% surf1=isosurface(leftLungtoCheck,0);
% p1 = patch(surf1);
% isonormals(leftLungtoCheck,p1);
% set(p1,'FaceColor','magenta','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
% daspect([1,1,1])
% view(3); axis tight
% camlight; lighting gouraud
% hold on
% 
% 
% surf1=isosurface(FinalNodule,0);
% p1 = patch(surf1);
% isonormals(FinalNodule,p1);
% set(p1,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
% daspect([1,1,1])
% view(3); axis tight
% camlight; lighting gouraud
% hold on
% 
% surf1=isosurface(upCarina,0);
% p1 = patch(surf1);
% isonormals(upCarina,p1);
% set(p1,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
% daspect([1,1,1])
% view(3); axis tight
% camlight; lighting gouraud
% hold on
% 
% surf1=isosurface(downCarina,0);
% p1 = patch(surf1);
% isonormals(downCarina,p1);
% set(p1,'FaceColor','green','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
% daspect([1,1,1])
% view(3); axis tight
% camlight; lighting gouraud
% hold on
% 
% 
% surf1=isosurface(VolumeRib,0);
% p1 = patch(surf1);
% isonormals(VolumeRib,p1);
% set(p1,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
% daspect([1,1,1])
% view(3); axis tight
% camlight; lighting gouraud
% hold on
%       

%% plot small and big

% surf1=isosurface(BigObj,0);
% p1 = patch(surf1);
% isonormals(BigObj,p1);
% set(p1,'FaceColor','green','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
% daspect([1,1,1])
% view(3); axis tight
% camlight; lighting gouraud
% hold on
% 
% 
% surf1=isosurface(SmallObj,0);
% p1 = patch(surf1);
% isonormals(SmallObj,p1);
% set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
% daspect([1,1,1])
% view(3); axis tight
% camlight; lighting gouraud
% hold on

% 
% s=load('C:\Users\MYCOM\Desktop\data.mat');SmallObj=s.SmallObj;BigObj=s.BigObj;
% temp_small=SmallObj;
% temp_small=temp_small(~BigObj);
% if any(temp_small(:))
%     disp('some part of SmallObj is outside BigObj')
% else
%     temp_big = bwperim(BigObj);%get the shell voxels
%     temp_small=SmallObj;
%     temp_small=temp_small(temp_big);
%     if any(temp_small(:))
%         disp('some part of SmallObj along the edge of BigObj')
%     else
%         disp('SmallObj is completely inside BigObj')
%     end
% end


%% 

% surf1=isosurface(Volumefissures,0);
% p1 = patch(surf1);
% isonormals(Volumefissures,p1);
% set(p1,'FaceColor','g','EdgeColor','none'); % set the color, mesh and transparency level of the surface
% daspect([1,1,1])
% view(3); axis tight
% camlight; lighting gouraud
% hold on