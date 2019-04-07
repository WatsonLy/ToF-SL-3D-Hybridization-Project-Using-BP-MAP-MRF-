clear all
close all


% vidv2 = videoinput('kinect',2,'Depth_640x480');
% vidv1 = videoinput('kinect',1,'RGB_1280x960');
% vidv2 = videoinput('kinect',1,'Infrared_640x480');
% vidv2 = videoinput('kinect',3,'BGR_1920x1080');
%vidv2 = videoinput('kinect',4,'Depth_512x424');
% src1 = getselectedsource(vidv1);
% src2 = getselectedsource(vidv2);

% vidv1.ROIPosition = [597 304 897 436];
% vidv2.ROIPosition = [126 132 304 152];

% vidv1.FramesPerTrigger = 1;
% vidv2.FramesPerTrigger = 1;


v1 = zeros(960, 1280);
v2 = zeros(480,640);
% v2 = zeros(424,512);
% v2original = zeros(1080,1920);
% v2distorted = zeros(424,512);


% start(vidv1);
% v1 = getdata(vidv1);
% stop(vidv1);
% delete(vidv1);

% start(vidv2);
% v2=getdata(vidv2);
% %Attempt to eliminate noise in IR stream by averaging
% sumFrame = zeros(424,512);
% 
% for i = 1:vidv2.TriggerRepeat
%     v2=getdata(vidv2);
%     if i == 1
%         sumFrame = double(v2);
%     else
%         sumFrame = sumFrame+double(v2);
%     end
%     v2 = zeros(480,640);
% end
% v2 = uint16(sumFrame / vidv2.TriggerRepeat);
% stop(vidv2);
% delete(vidv2);

% v1 = imresize(v1,[480 640]);
v1 = imread('InfraredV1.tif');
v2 = imread('InfraredV2.tif');
% v1 = imread('V1infra1.tiff');
% v2 = imread('v2infra1.tif');

%v1 = rgb2gray(v1);
%v2 = rgb2gray(v2);

%Histogram equalization
% v1 = histeq(v1);
% v2 = histeq(v2);

%Gamma transformation (v2 = 0.2, v1 = 0.4)
v1 = imadjust(v1 ,[], [], 0.36);
v2 = imadjust(v2 ,[], [], 0.21);

%Log transformation
% v1 = c*(log(1+v1));
% v2 = c*(log(1+v2));

%Experimental RGB-IR transformation matrix for V2
%https://link.springer.com/content/pdf/10.1007%2F978-3-662-47487-7_17.pdf
T = [0.3584, 0.0089, 0.000;0.0031,0.3531,0.0001;-101.5934,13.6311,0.9914];
% T = [0.3584, 0.0089, 0.000;0,0.3531,0.0001;-101.5934,13.6311,1];
% tform = projective2d(T);
% outputView = imref2d(size(v2));
% Ir = imwarp(v1,tform,'OutputView',outputView);
% figure(2)
% imshow(Ir); 
% title('Recovered image');

figure(1)
subplot(1,2,1)
imshow(v1);
subplot(1,2,2)
imshow(v2);


% Detects matching features from both images (V1 and V2) and stores them in
% ptsOriginal, ptsDistorted.
ptsOriginal  = detectSURFFeatures(v1);
ptsDistorted = detectSURFFeatures(v2);
[featuresOriginal,validPtsOriginal] = ...
    extractFeatures(v1,ptsOriginal);
[featuresDistorted,validPtsDistorted] = ...
    extractFeatures(v2,ptsDistorted);




%Matches detected surface points
index_pairs = matchFeatures(featuresOriginal,featuresDistorted, "MatchThreshold", 1, "MaxRatio", 0.9);
matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));

% f = cpselect(v1,v2);

%Displays matched points including outliers
figure; 
showMatchedFeatures(v1,v2,...
    matchedPtsOriginal,matchedPtsDistorted);
title('Matched SURF points,including outliers');



%Checks matched surface points and chooses valid pairs
% Outputs homography matrix for registration
[tform,inlierPtsDistorted,inlierPtsOriginal] = ...
    estimateGeometricTransform(matchedPtsDistorted,matchedPtsOriginal,...
    'projective', "MaxNumTrials", 1000, "Confidence", 99, "MaxDistance",6.5);
figure; 

%Displays inlier matched points
showMatchedFeatures(v1,v2,...
    inlierPtsOriginal,inlierPtsDistorted);
title('Matched inlier points');

T = [1.5277089	-0.026776986 -0.0000358;	-0.010839382 1.5996544 -0.0001383637;-76.933472	-61.545399	1];
tform = projective2d(T);
 

outputView = imref2d(size(v1));
Ir = imwarp(v2,tform,'OutputView',outputView);
figure; imshow(Ir); 
title('Recovered image');


