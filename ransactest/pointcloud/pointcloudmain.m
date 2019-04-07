clear;
close all;

v1 = imread('Cropped KinectV1.tif');
v2 = imread('Cropped KinectV2.tif');
restored = imread('Restored.tif');

imwrite(histeq(v1), "HistEqV1.tif");
imwrite(histeq(v2), "HistEqV2.tif");
imwrite(histeq(restored), "HistEqRestored.tif");

d1 = abs(v1-v2);


[n,m] = size(v1);
v1 = padarray(v1,[480-n,640-m],0,'post');
v2 = padarray(v2,[480-n,640-m],0,'post');
restored = padarray(restored,[480-n,640-m],0,'post');

d1 = padarray(d1,[480-n,640-m],0,'post');
% image = uint16(zeros(640,480));

% image(1:n,1:m) = uint16(v1);



% most of the following was taken from https://www.mathworks.com/help/vision/ref/pcfromkinect.html
depthDevice = imaq.VideoDevice('kinect',2);
depthImage = step(depthDevice);

ptCloud = pcfromkinect(depthDevice,v1);
ptCloud2 = pcfromkinect(depthDevice,v2);
ptCloud3 = pcfromkinect(depthDevice,restored);
ptCloud4 = pcfromkinect(depthDevice,d1);

% player = pcplayer([0,1.273803949356079],[-0.955353021621704,0],[0,2.289000034332275]);
 player = pcplayer(ptCloud.XLimits,ptCloud.YLimits,ptCloud.ZLimits,...
        'VerticalAxis','y','VerticalAxisDir','down');
xlabel(player.Axes,'X (m)');
ylabel(player.Axes,'Y (m)');
zlabel(player.Axes,'Z (m)');

player2 = pcplayer(ptCloud2.XLimits,ptCloud2.YLimits,ptCloud2.ZLimits,...
        'VerticalAxis','y','VerticalAxisDir','down');
xlabel(player2.Axes,'X (m)');
ylabel(player2.Axes,'Y (m)');
zlabel(player2.Axes,'Z (m)');

player3 = pcplayer(ptCloud3.XLimits,ptCloud3.YLimits,ptCloud3.ZLimits,...
        'VerticalAxis','y','VerticalAxisDir','down');
xlabel(player3.Axes,'X (m)');
ylabel(player3.Axes,'Y (m)');
zlabel(player3.Axes,'Z (m)');

player4 = pcplayer(ptCloud4.XLimits,ptCloud4.YLimits,ptCloud4.ZLimits,...
        'VerticalAxis','y','VerticalAxisDir','down');
xlabel(player4.Axes,'X (m)');
ylabel(player4.Axes,'Y (m)');
zlabel(player4.Axes,'Z (m)');
 


view(player,ptCloud);
view(player2,ptCloud2);
view(player3,ptCloud3);
view(player4,ptCloud4);




release(depthDevice);

