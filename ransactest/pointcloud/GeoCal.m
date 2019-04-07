vid = videoinput('kinect', 2, 'Depth_640x480');
src = getselectedsource(vid);

vid.FramesPerTrigger = 50;
vid.ROIPosition = [318 238 42 35];

start(vid);
v1 = getdata(vid);
stop(vid);

vid2 = videoinput('kinect', 4, 'Depth_512x424');
src = getselectedsource(vid2);

vid2.FramesPerTrigger = 50;
vid2.ROIPosition = [252 192 23 21];

start(vid2);
v2 = getdata(vid2,50);
stop(vid2);

% for a = 1:50
% 
%     if a == 1
%         temp1 = v1(:,:,1,a);
%     else
%         temp1 = [temp1; v1(:,:,1,a)];
%     end
% end
% 
% 
% for a = 1:50
%     
%     if a == 1
%         temp2 = v2(:,:,1,a);
%     else
%         temp2 = [temp2; v2(:,:,1,a)];
%     end
% end
% 
% resultv1 = typecast(temp1(:), "double");
% resultv2 = typecast(temp2(:), "double");


C = 36.015;
if(v2 > C)
    v2 = v2 - C;
end

% v1mean = [v1a, v1b, v1c];
v1mean = mean(mean(mean(v1)))
v1std = std2(v1)
% v1std = std(v1mean)

% v2mean = [v2a, v2b, v2c];
v2mean = mean(mean(mean(v2)))
v2std = std2(v2)
% v2std = std(v2mean)

diff = abs(v1mean - v2mean)