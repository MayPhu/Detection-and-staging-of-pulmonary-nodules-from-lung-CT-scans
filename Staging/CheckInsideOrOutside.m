function [result] = CheckInsideOrOutside(BigObj,SmallObj)
%% This function check a small 3D image is inside or outside of a big 3D image using convex hull.

% Check the convex hull of big object ( be area only one label martix) 

convexhulllBig = regionprops3(BigObj,'ConvexHull');
convexhulllBig = table2array(convexhulllBig.ConvexHull);

%%% sort the convexhull by X coordinate ( largest and smallest X coordinate pairs are extracted);
[~,idx] = sort(convexhulllBig(:,1));
sortedbyX = convexhulllBig(idx,:);
BigXmin = sortedbyX (1,1);
BigXmax = sortedbyX (end,1);

%%% sort the convexhull by Y coordinate ( largest and smallest Y coordinate pairs are extracted);
[~,idx] = sort(convexhulllBig(:,2));
sortedbyY = convexhulllBig(idx,:);
BigYmin = sortedbyY (1,2);
BigYmax = sortedbyY (end,2);

%%%  sort the convexhull by Z coordinate ( largest and smallest Z coordinate pairs are extracted);
[~,idx] = sort(convexhulllBig(:,3));
sortedbyZ = convexhulllBig(idx,:);
BigZmin = sortedbyZ (1,3);
BigZmax = sortedbyZ (end,3);

%% 
% Check the convex hull of big object ( be area only one label martix) 
 
convexhulllSmall = regionprops3(SmallObj,'ConvexHull');
convexhulllSmall = table2array(convexhulllSmall.ConvexHull);

%%% sort the convexhull by X coordinate ( largest and smallest X coordinate pairs are extracted);
[~,idx] = sort(convexhulllSmall(:,1));
sortedbyX = convexhulllSmall(idx,:);
SmallXmin = sortedbyX (1,1);
SmallXmax = sortedbyX (end,1);

%%% sort the convexhull by Y coordinate ( largest and smallest Y coordinate pairs are extracted);
[~,idx] = sort(convexhulllSmall(:,2));
sortedbyY = convexhulllSmall(idx,:);
SmallYmin = sortedbyY (1,2);
SmallYmax = sortedbyY (end,2);

%%%  sort the convexhull by Z coordinate ( largest and smallest Z coordinate pairs are extracted);
[~,idx] = sort(convexhulllSmall(:,3));
sortedbyZ = convexhulllSmall(idx,:);
SmallZmin = sortedbyZ (1,3);
SmallZmax = sortedbyZ (end,3);


%% compare coordinate values 

if SmallXmin> BigXmin && SmallXmax< BigXmax && SmallYmin>BigYmin && SmallYmax<BigYmax && SmallZmin>BigZmin && SmallZmax<BigZmax
    result= 1;
else 
    result = 0;
end






end

