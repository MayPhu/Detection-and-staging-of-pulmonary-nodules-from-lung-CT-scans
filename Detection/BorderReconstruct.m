function [ RefinedParenchyma ] = BorderReconstruct( SegmentedParenchyma )
rng(1);
border = imfill(SegmentedParenchyma,'holes');
border =edge(border,'Canny');
% border=imdilate(border,strel('disk',1));
% figure,imshow(border,[]),title('Border');


cc = bwconncomp(border,8);
l = labelmatrix(cc);
numlobe= cc.NumObjects;  


if numlobe == 0
      RefinedParenchyma = SegmentedParenchyma;
else       
       
     for li=1:numlobe
        currentlobe = (l==li);
%         figure(),imshow(currentlobe),title('Lobe');
        boundaries = bwboundaries(currentlobe,8);
        %% trace the pixels values in current lobe boundary      
        
        for kb = 1 : length(boundaries)
        % Get the k'th boundary.
        thisBoundary = boundaries{kb};
        % Get the x and y coordinates.
        x = thisBoundary(:, 2);
        y = thisBoundary(:, 1);
        % Plot this boundary over the original image.
        hold on;
        plot(x, y, 'w-', 'LineWidth', 0.5);
        
        for p = 2 : length(y)
            %  Horizontal direction-to-code convention is:  1  0  -1
        %                                                \ | /
        %                                             1 -- P -- -1
        %                                                / | \
        %                                               1  0  -1

        %    Vertical direction-to-code convention is:  1  1  1
        %                                                \ | /
        %                                             0 -- P -- 0
        %                                                / | \
        %                                              -1 -1  -1
        coordinate{p-1} = [x(p);y(p)];
        
        
            if x(p) < x(p-1) && y(p) < y(p-1)
			% Moved to the upper left
			HChainCode(p) = -1;
            VChainCode(p) = 1;
            
            elseif x(p) == x(p-1) && y(p) < y(p-1)
			% Moved straight up
			HChainCode(p) = 0;
            VChainCode(p) = 1;
            
            elseif x(p) > x(p-1) && y(p) < y(p-1)
			% Moved to the upper right
			HChainCode(p) = 1;
            VChainCode(p) = 1;
            
            elseif x(p) > x(p-1) && y(p) == y(p-1)
			% Moved right
			HChainCode(p) = 1;
            VChainCode(p) = 0;
		
            elseif x(p) > x(p-1) && y(p) > y(p-1)
			% Moved right and down.
			HChainCode(p) = 1;
            VChainCode(p) = -1;
		
            elseif x(p) == x(p-1) && y(p) > y(p-1)
			% Moved down
			HChainCode(p) = 0;
            VChainCode(p) =-1;
            
            elseif x(p) <= x(p-1) && y(p) > y(p-1)
			% Moved down and left
			HChainCode(p) = -1;
            VChainCode(p) = -1;
            
            elseif x(p) <= x(p-1) && y(p) == y(p-1)
			% Moved left
			HChainCode(p) = -1;
            VChainCode(p) = 0;
            end
        
        end
        
        % Save the chain code that we built up for this boundary
        % into the cell array that will contain all chain codes.
        Coordinates{kb}= cell2mat(coordinate);
        HchainCodes{kb} = HChainCode;    
        VchainCodes{kb} = VChainCode;
        
        %% Gaussian low-pass filter to reduce noise
        %Dilate the edges
        HGaussianChainCodes{kb}= imdilate(HchainCodes{kb},strel('disk',10));
        HGaussianChainCodes{kb} = round(imfilter(HGaussianChainCodes{kb},fspecial('gaussian')));
        VGaussianChainCodes{kb}= imdilate(VchainCodes{kb},strel('disk',10));
        VGaussianChainCodes{kb} = round(imfilter(VGaussianChainCodes{kb},fspecial('gaussian')));
        
        %% Differential operation to calculate infection points
        HDifferentiate{kb}= diff( HGaussianChainCodes{kb});
        VDifferentiate{kb}= diff( VGaussianChainCodes{kb});
        
        end
       
        %% Find infection points ( non-zero vlaues)
        HInfection = [cell2mat(HDifferentiate);cell2mat(Coordinates)];
        VInfection = [cell2mat(VDifferentiate);cell2mat(Coordinates)];
        
        
        %% Horizontal Infection Points
        Xindices = find(HInfection(1,:)== 0);
        HInfection(:,Xindices) = [];
        % to know the coordinate points
        Hx= HInfection(2,:);
        Hy=HInfection(3,:);
        Hinfections= [Hx' Hy'];
        
       %% Vertical Infection Points
        Yindices = find(VInfection(1,:)==0);
        VInfection(:,Yindices) = [];
        Vx= VInfection(2,:);
        Vy=VInfection(3,:);
        Vinfections= [Vx' Vy'];
        
        %% Show Infections Points on the image
%         hold on;
%         plot(Hinfections(:,1), Hinfections(:,2), 'b*'); 
%         plot(Vinfections(:,1),Vinfections(:,2), 'g*');

        %% Assume all infectin points as data and perform Kmeans Clustering      
        x = Hinfections ; 
        y = Vinfections ;
        Input = [x;y];
        H=Input(:,1);
        V=Input(:,2);
        % apply kmeans 
        
        %% Finding optimal K value
        myfunc = @(X,K)(kmeans(X, K));
        [r c]=size(Input);
        % krange= round(sqrt(r)/2);
        % eva = evalclusters(x,clust,criterion)
        eva = evalclusters(Input,myfunc,'silhouette',...
        'klist',[1:3]);
%         k=eva.OptimalK;
k=3;

        %% Perform K-mean Clustering
        [idx,C,sumd,D] = kmeans(Input,k) ;
        % get each cluster
        data = cell(k,1) ;
        distancebetweenpoints=cell(k,1);
        tightness= cell(k,1);
        ratio = cell(k,1);
        variation= cell(k,1);
        convex= cell(k,1);
        sumcon= cell(k,1);
        area= cell(k,1);
        tmp = zeros;
        mrk={'o','s','*','v','+','^','o','s','*','v','+','^'};
        
        %% Plotting according to the cluster
        for i = 1:k
        data{i} = [H(idx==i),V(idx==i)] ; 
        %[convex{i},area{i}]= convhull(data{i});
        sumcon{i}= sum(cell2mat(convex(i,:)));
        distancebetweenpoints{i}=squareform(pdist(data{i}));
        dtmp = cell2mat(distancebetweenpoints(i,:));
        maxdistance = round(max(dtmp(:)));
        [row col]=size(dtmp);
        avgdistance= round(sum(dtmp(:))/ row);
        ratio{i} = round([sumd(i,:)/row]);
        tightness{i}= [maxdistance avgdistance];
        variation{i}= sum(var(data{i}));    
        tmp = cell2mat(data(i,:)); 
        plot(tmp(:,1),tmp(:,2),'-');
        hold on
        
        %% Reduce Infection points 
        % [min_num,min_idx] = min(cell2mat(sumcon))
        % data(min_idx, :) = [];
        for i= 1:k
        tmp = cell2mat(data(i,:)); 
        % plot(tmp(:,1),tmp(:,2), mrk{i});
        end
        cluster1= cell2mat(data(1));
        cluster2=cell2mat(data(2));
        X= ([cluster1(:,1)]);
        Y=([cluster1(:,2)]);
        % X= ([cluster1(:,1);cluster2(:,1)]);
        % Y=([cluster1(:,2);cluster2(:,2)]);
        [r c]= size (X);
        % for i=1:r
        %     for j=1:r-1
        %         hLine = imline(gca,[X(i,:) X(j+1,:)],[Y(i,:) Y(j+1,:)]);
        %        
        %     line ([X(i,:),X(j+1,:)],[Y(i,:),Y(j+1,:)],...
        %      'Color','b')
        %     end
        % end
        end
hold off
F = getframe ;
RL = F.cdata;
RL= rgb2gray(RL);
RL = imresize(RL, [512 512]);
RL(:, 1) = zeros(1, 512);
RL(512,:)= zeros(1,512);
mask= imfill(logical(RL),'holes');
refinemask{li}= mask;
     end
clearvars -except li refinemask numlobe
leftlobemask =cell2mat(refinemask(:,1));
if numlobe==1
      rightlobemask =zeros(size(leftlobemask));
 else
rightlobemask =cell2mat(refinemask(:,2));
 end
RefinedParenchyma= imadd(double(leftlobemask),double(rightlobemask));
close all;
end
end

