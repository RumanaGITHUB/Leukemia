close all;
clear all;

% rgb = imread('bld.jpg');
%   for z = 1:200
[filename pathname] = uigetfile({'*.jpg';'*.jpeg''*.bmp';'*.png'},'Select Microscopic image');
inputimage=strcat(pathname, filename);
rgb= imread(inputimage);
% %  imshow(rgb);

myimage=imread(inputimage);

%% Blood cells count

%%Extracting the blue plane 
bPlane = myimage(:,:,3)  - 0.5*(myimage(:,:,1)) - 0.5*(myimage(:,:,2));
%figure
%imshow(bPlane)
%title('Extracted White Blood Cells','fontsize',14);
%%Extract out purple cells
%figure
BW = bPlane > 29;
%imshow(BW)
%%Remove noise 100 pixels or less
BW2 = bwareaopen(BW, 100);
%imshow(BW2)
%%Calculate area of regions
L = bwlabel(BW2);

%% create new figure to output superimposed images
% first display the original image
% %  figure; 
% %  imshow(myimage); hold on

%% Label connected components
[L, Ne]=bwlabel(BW2);
propied=regionprops(L,'BoundingBox'); 
 himage = imshow(BW2);

%% Get the total number of cells that have been added with bounding box
whitecount = size(propied,1);

%% Added bounding box to the white blood cells
hold on
for n=1:whitecount
  rectangle('Position',propied(n).BoundingBox,'EdgeColor','g','LineWidth',2)
end
hold off

%% Superimpose the two image
set(himage, 'AlphaData', 0.5);

%% Output total white blood cells onto the figure
title(sprintf('%i White Blood Cells Detected',whitecount),'fontsize',14);


%% RED BLOOD CELLS
%% Extracting the red plane 
rPlane = myimage(:,:,1)- 0.4*(myimage(:,:,3)) - 0.6*(myimage(:,:,2));
%figure
%imshow(rPlane)
%title('Extracted Red Blood Cells','fontsize',14);
%% Extract out red cells
BWr = rPlane > 19;
%figure
%imshow(BW)
%%Remove noise 100 pixels or less
BWr2 = bwareaopen(BWr, 100);
%imshow(BW2)
%%Calculate area of regions
cellStatsr = regionprops(BWr2, 'all');
cellAreasr = [cellStatsr(:).Area];

%% create new figure to output superimposed images
% first display the original image
% % figure;
% % imshow(myimage); 
hold on

%% Label connected components
[Lr, Ner]=bwlabel(BWr2);
propiedr=regionprops(Lr,'BoundingBox'); 
himager = imshow(BWr2);

%% Get the total number of cells that have been added with bounding box
redcount = size(propiedr,1);

%% Added bounding box to the red blood cells
hold on
for n=1:redcount
  rectangle('Position',propiedr(n).BoundingBox,'EdgeColor','r','LineWidth',2)
end
hold off

%% Superimpose the two image
set(himager, 'AlphaData', 0.5);

%% Output total red blood cells onto the figure
title(sprintf('%i Red Blood Cells Detected',redcount),'fontsize',14);

%% Calculate percentages
totalCells = whitecount + redcount;
wbcPercent = (whitecount ./ totalCells) .* 100;
rbcPercent = (redcount ./ totalCells) .* 100;

disp('total blood cells percentage:')
disp(totalCells)
disp('white blood cells percentage:')
disp(wbcPercent)
disp('red blood cells percentage:')
disp(rbcPercent )

if vpa(wbcPercent) >= 20
    disp( 'POTENTIAL LEUKEMIA DETECTED');
else
   disp( 'Normal');
 end



%% extracting leukemia cell

% pause

[centers, radii] = imfindcircles(rgb,[35 60],'ObjectPolarity','dark','Sensitivity',0.9)
% %  imshow(rgb);
% pause
h = viscircles(centers,radii);
cell=length(centers);
% pause
red=rgb(:,:,1);green= rgb(:,:,2); blue= rgb(:,:,3);
%p=impixel(rgb);
out=red>25 & red<123 &green<135 & blue>167 & blue<201;
% %  imshow(out);
% pause

% % %  hole filling and morphological operation
subplot(2,2,1);
% %  imshow(out);
title('extracted cell');
subplot(2,2,2);
out1=imfill(out,'holes');
% %  imshow(out1);
title('hole filling');
% pause
subplot(2,2,3);
out2=bwmorph(out1,'erode');
% %  imshow(out2);
title('eroding');
subplot(2,2,4);
% pause
out3=bwmorph(out2,'dilate',4);
out3=imfill(out3,'holes');
% %  imshow(out3);
title('dilation');
out4=out3;
% %  figure
% %  imshow(out4)
title('dilation and hole fill');

gray_image = rgb2gray(rgb);
% imshow(gray_image);
threshold=gray_image;
bw=BW2;
%Feature


glcm = graycomatrix(threshold);
stats = graycoprops(glcm,'Contrast Correlation Energy Homogeneity');
Contrast = stats.Contrast;
Correlation = stats.Correlation;
Energy = stats.Energy;
Homogeneity = stats.Homogeneity;
Entropy = entropy(bw);
Mean = mean2(bw);
Standard_Deviation = std2(bw);
Variance = mean2(var(double(bw)));
s = sum(double(bw(:)));
Smoothness = 1-(1/(1+s));
Kurtosis = kurtosis(double(bw(:)));
Skewness = skewness(double(bw(:)));
coefficient_of_variance=Standard_Deviation/Mean;

RMS=rms(bw*1);

%R=randn(bw*1);
Moment=moment(bw,5);
% Eigen=eig(bw*1);
% x1=mean(Eigen);


%%%% Create Training Set

group= [Contrast,Correlation,Energy,Homogeneity, Mean, Standard_Deviation, Entropy,Variance, Smoothness, Kurtosis, Skewness,coefficient_of_variance ,RMS, Moment];


%  group(z,:) = [Contrast,Correlation,Energy,Homogeneity, Mean, Standard_Deviation, Entropy,Variance, Smoothness, Kurtosis, Skewness,coefficient_of_variance ,RMS, Moment ];
%     
%      str =  0;
%      str1 = 1;
%      
% 
%      if z<120
%          class(z,:) = str;
%      
%      else 
%          class(z,:)= str1;
%      
%           end
%     end
%      
%      save('feature.mat', 'group','class');
    Train1 = array2table(group)

%   Train.class = categorical(class)

   load acute.mat
 yfit = acute.predictFcn(Train1);
  t1=double(yfit);
  t2=string(yfit);
  
   if t2 == '1'
     disp('AML')   
 end
 if t2 == '0'
     disp('ALL')   
  end
%  hom=Homogeneity;
% % c=Contrast;
% en=Energy;
% co=Correlation;
% e=Entropy;
% fismat=readfis('cancer.fis')
% input=[hom en co e];
%   output=evalfis(input,fismat)
% if(output>=0)&&(output<=0.4)
%     disp('ALL')
% else(output<=0.4)&&(output>=0.7)
%      disp('AML')
% end
%   end
% if (output1>=output2)
%     disp('AML')
% 
% else 
%     disp('ALL')
% end