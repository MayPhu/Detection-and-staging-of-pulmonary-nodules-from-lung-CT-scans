function [ mask ] = SegmentParenchyma( I )
%This function trys to segment the Parenchyma from the lung CT scan image.
%   I= lung ct scan (Dicomformat);
%% Setting Threshold value using CT value 
thresholdHU = -375;  
thresholdGray = thresholdHU + 1000;

%% Generates the first mask using Gray Threshold.
mask = false(size(I));
mask(find(I < thresholdGray)) = 1;

mask = imclearborder(mask); 
border = edge(mask); % guess the border of the parenchyma
% figure,imshow(border);

%% find the coordinate of first white pixel
[row, column] = find(border, 1, 'first');
%  Assign the initial seeplot(binLocations, pdf, 'b-', 'LineWidth', 2);
grid on;d point          
    seedpointR = row;    
    seedpointC = column; 
 


%% Applying Region growing 
W = graydiffweight(I, seedpointC, seedpointR, 'GrayDifferenceCutoff', 8);
thresh = 0.01;
[mask, D] = imsegfmm(W, mask, thresh);


%% Change from double
mask = int16(mask);
mask = imclearborder(mask);


%% Remove trachea 
% mask = bwareafilt(logical(mask),2,'largest');

% figure,imshow(removeTrachea,[]),title('ClearBorder and Remove Trachea Image');

% figure (2),imshow(mask,[]),title('Binary mask');


%% Segmented parenchyma
% Parenchyma =double(mask).* double(I);
% figure (3),imshow(extracted1,[]),title('Segmented Parenchyma');


end