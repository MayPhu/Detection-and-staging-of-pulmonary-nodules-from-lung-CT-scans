unction [ features ] = FeaturesExtraction( V,xSpacing,ySpacing,zSpacing )
%V is the volume image

CCCC = bwconncomp(V, 26);
LLL = labelmatrix(CCCC);
 for i = 1:length(CCCC.PixelIdxList)
  
     a = (LLL==i);    
    
    %% Shape Features     
    status = regionprops(a,'Area','Centroid');
    V= status.Area*xSpacing*ySpacing*zSpacing*1e-6;
    Lable(i,:) = i;
    Volume(i,:) = V;
    Centroid (i,:) = status.Centroid;
    
    
    %perimeter
    boundary= bwperim(a,26);
    p=nnz(boundary);
    Perimeter(i,:) = p;
    
    % Radius 
    tmp= 3*V/4*pi;
    R(i,:) = p* nthroot(tmp,3);
    
    
    SA= imSurface(a);
    SurfaceArea(i,:) = SA;
    Density(i,:) = imVolumeDensity(a);
    

    
    % Sphericity V
    Sphericity(i,:) = V/SA;
    
    status2= regionprops3(a,'AllAxes','Eccentricity');
    Ax=status2.FirstAxisLength;
    Ay=status2.SecondAxisLength;
    Az=status2.ThirdAxisLength;
    
    tmpMajor = max(Ax,Ay);
    major=max(tmpMajor,Az);
    MajorAxis(i,:)= major;
    
    tmpMinor = min(Ax,Ay);
    minor=min(tmpMajor,Az);
    MinorAxis(i,:)= minor;
    
   Elongation (i,:) = minor/major;
   Eccentricity(i,:) =status2.EquatorialEccentricity;
     
     
       
   %% Histogram Features
   Mean(i,:) = mean2(a);
   Std(i,:) = std2(a);
   % Skewness and Kurtiosis
   [pixelCounts GLs] = imhist(a(:)); 
   [skew kurtosis] = GetSkewAndKurtosis(GLs, pixelCounts);
   Skewness (i,:) = skew;
   Kurtosis(i,:)= kurtosis;
   
   end
%     
 features = [Lable Volume Centroid Perimeter R SurfaceArea Density Sphericity MajorAxis MinorAxis Elongation Eccentricity Mean Std Skewness Kurtosis];


end