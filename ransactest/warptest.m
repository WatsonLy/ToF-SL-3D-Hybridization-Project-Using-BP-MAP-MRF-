clear all
close all
J=zeros(2048);
I = imread('lena.tif');
figure(1)
subplot(1,2,1)
imshow(I)
% H=[-0.99,-0.24,276;-0.91,-0.21,250;-0.003,-.000836,1];
% H=[-0.737116997695806,-0.586938027527255,2.795264252743421e+02;-0.706864083490700,-0.538386725962351,2.625341743720673e+02;-0.002691747751094,-0.002062544269873,1];
H=[-0.498287107387665,-0.399879290357225,2.075598456517095e+02;-0.795818279837687,-0.610727015825541,3.232719773994849e+02;-0.002441006229952,-0.001899934850893,1]
A=transpose(inv(H));
tform = projective2d(A);

J = imwarp(I,tform);
subplot(1,2,2)
imshow(J)