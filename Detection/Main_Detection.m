clc;
clear;
close all;

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
for q = 1:n
     fileName = ['F:\Experiment\LIDC-IDRI-0001\Write\',num2str(q),'.dcm'];
     I=dicomread(fileName);
     %% Segment the Parenchyma
     segParenchyma = SegmentParenchyma(I);
     %figure(),imshow(segParenchyma,[]),title('Segmented Parenchyma');

    %% Border Reconstruction by bidirectional chain code
    BorderReconstructed = BorderReconstruct( segParenchyma );
    %figure(),imshow(BorderReconstructed,[]),title('Border Refined Segmented Parenchyma');

    %% Nodule candidate detection 
    finalsegParenchyma = double(segParenchyma) .* double(I);
    % figure(),imshow(finalsegParenchyma,[]),title('Final egmented Parenchyma');
    
    thresholdHU = -375;  
    thresholdGray = thresholdHU + 1000;%     
    InternalObjects = finalsegParenchyma>thresholdGray;
    
   
    finalcandidate = double(InternalObjects).*double(I); % This is the final nodule candidates
    lungVolume(:,:,q) = finalcandidate;    % create 3D structure of the nodule candidates
  
    
end
figure(),isosurface(lungVolume);  % show 3D structure of the nodule candidates


%% Preliminary screening   (this is the task of removing objects having high possiblity to be non-tumors)

% Eliminate blobs by checking area  (very small and very large objects may
% not be tumors)
        LB =10;
        UB= 3700;
        RemoveSmallblob = xor(bwareaopen(lungCandidate,LB,26),  bwareaopen(lungCandidate,UB,26));
        %figure(),isosurface(RemoveSmallblob,0),title('Nodule candidates after eliminating by Area');
    
    % Eliminate blobs by checking eccentricity (very tubular objects may
    % not be tumors)
        CC = bwconncomp(RemoveSmallblob, 26);
        labelVolume2 = bwlabeln(RemoveSmallblob, 26);
        LL= labelmatrix(CC);
        NumObj2= CC.NumObjects

        for i=1 :NumObj2 
            a= (LL==i);
            status2= regionprops3(a,'AllAxes');
            Ax=status2.FirstAxisLength;
            Ay=status2.SecondAxisLength;
            Az=status2.ThirdAxisLength;
    
            tmpMajor = max(Ax,Ay);
            major=max(tmpMajor,Az);  
    
            tmpMinor = min(Ax,Ay);
            minor=min(tmpMajor,Az);
            elong = minor/major;
            Elongation (i,:) = [i elong];
        end
        %slect blobs that are greater than 0.8 in  eccentricity
            lines=  Elongation(Elongation(:,2)<0.3,:);
            lineslabel= lines(:,1);
            RemovedLine = ismember(labelVolume2, lineslabel)>0;
            %figure,isosurface(RemovedLine,0),title('Nodule candidates after eliminating by Eccentricity');

        % Eliminate blobs by checking VRR 
         CCC = bwconncomp(RemovedLine, 26);
        labelVolume3 = bwlabeln(RemovedLine, 26);
         LLL= labelmatrix(CCC);
         NumObj3= CCC.NumObjects%
         
        for i=1 :NumObj3 
            obj= (LLL==i);
            status3 = regionprops(obj,'Area');
            NV(i,:) = status3.Area;   
            skel{i} = Skeleton3D(obj);
            sobj= skel{i};   
            status4= regionprops(sobj,'Area'); 
            size=length(status4);
%           row=size(:,1);
            if size==0
            NVS(i,:)= 0; 
            else
            NVS(i,:)= status4.Area;
            end
            VRRTemp = 1- (NVS(i,:)./ NV(i,:));
            VRR(i,:)= [i VRRTemp];
    
        end

        %slect blobs that are greater than 0.9 in  VRR rate
        vessles=  VRR(VRR(:,2)<0.95,:);
        vessellabel= vessles(:,1);
        RemovedVRR = ismember(labelVolume3, vessellabel)>0;
        figure,isosurface(RemovedVRR,0),title('Nodule candidates after eliminating by VRR Rates');
    




%% Juxta-vascular nodule detection (To cut the tumors from attached vessels )

I = permute(RemovedVRR,[3,1,2]);
I = I(end:-1:1,:,:);

%normalize input a little bit
I = I - min(I(:));
I = I / prctile(I(I(:) > 0.5 * max(I(:))),90);
I(I>1) = 1;

% compute enhancement for two different tau values
[V,I1,I2,I3] = blobness3D(I, 1:5, [1;1;1], 0.75, false);
figure(),isosurface(V);



%% Feature extraction  (nodule features)
    
    title = {'Lable' 'Volume' 'X' 'Y' 'Z' 'Perimeter' 'Radius' 'SurfaceArea' 'Density' 'Sphericity' 'MajorAxis' 'MinorAxis' 'Elongation' 'Eccentricity' 'Mean' 'Std' 'Skewness' 'Kurtosis'};
    extractedfeatures = FeaturesExtraction(V,xSpacing,ySpacing,zSpacing );
    FinalFeatures= [cellstr(title); num2cell(extractedfeatures)];
    xlswrite(xlsfilename,FinalFeatures);
    save(workspacefilename)




