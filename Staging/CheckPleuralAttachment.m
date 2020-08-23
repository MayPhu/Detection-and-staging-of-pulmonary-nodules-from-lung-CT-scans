function [result] = CheckPleuralAttachment(RightLung,LeftLung,Nodule)
%% be sure each of the three inputs is a label matrix (only one object) 

convexhulllRightLung = regionprops3(RightLung,'ConvexHull');
convexhulllRightLung = table2array(convexhulllRightLung.ConvexHull);

convexhulllLeftLung = regionprops3(LeftLung,'ConvexHull');
convexhulllLeftLung = table2array(convexhulllLeftLung.ConvexHull);

convexhulllNodule = regionprops3(Nodule,'ConvexHull');
convexhulllNodule = table2array(convexhulllNodule.ConvexHull);


RightPleuralAttach  = intersect(convexhulllRightLung,convexhulllNodule,'rows');
LeftPleuralAttach  = intersect(convexhulllLeftLung,convexhulllNodule,'rows');

if ~isempty(RightPleuralAttach) || ~isempty(LeftPleuralAttach)
result = 1;
else 
result = 0;
end



end

