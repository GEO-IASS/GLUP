
clear all ; close all ; clc; 
 
load('200pxl8end40dB1.mat');

rho = 1;
mu = 1;
epsabs = 1E-2;
epsrel = 1E-2;
 
t = clock;
[X, Niter] = GLUP(Data,Data,rho,mu,epsabs,epsrel);
t = etime(clock,t);  

figure;
imshow(mat2gray(X));
set(gca,'XAxisLocation','top','ydir','reverse');
xlabel('Num. of column','Interpreter','Latex','FontSize',14)
ylabel('Num. of line','Interpreter','Latex','FontSize',14)
title(['Grayscale image of estimated abundance matrix $X$, SNR=' num2str(40) 'dB.'],'Interpreter','Latex','FontSize',16);
colormap(bone)
colormap(flipud(colormap))
colorbar

