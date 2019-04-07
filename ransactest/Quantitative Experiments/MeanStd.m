r1 = imread('Restored1.tif');
r2 = imread('Restored2.tif');
r3 = imread('Restored2.tif');
r4 = imread('Restored3.tif');
r5 = imread('Restored4.tif');

res = [r1;r2;r3;r4;r5];
res = res(res>0);

resMean = mean(mean(res))
resStd = std2(res)