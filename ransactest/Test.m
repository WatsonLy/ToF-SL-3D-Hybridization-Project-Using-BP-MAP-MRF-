vid1 = videoinput('kinect', 1, 'RGB_1280x960');
vid2 = videoinput('kinect', 3, 'BGR_1920x1080');

vid1.ReturnedColorspace = 'grayscale';
vid2.ReturnedColorspace = 'grayscale';

vid1.FramesPerTrigger = 3;
vid2.FramesPerTrigger = 3;

start(vid1);
start(vid2);
v1 = getdata(vid1);
v2 = getdata(vid2);
stop(vid1);
stop(vid2);
