
clear all ; close all ; clc;

%% Load HSI
demoFile =   'Pavia.mat'; 
vars =  {'pavia'};
load(demoFile, vars{:});
M3d = pavia; 
M3d = M3d(:,224:end,:); % keep right part of image
[nl,nc,p] = size(M3d);
nt = nc*nl;
M = reshape(M3d, nt, p).';


%% Normalise HSI data cube
temp = max(M(:));
M = M./temp;


%% Initialise H by taking 200 samples from original image
n1 = 200;
tmp = randperm(nt);
H = M(:,tmp(1:n1));
H_idx = tmp(1:n1)';


%% GLUP
mu = 10; % regularization parameter 
rho =  1; % penalty parameter for ADMM 
epsabs = 1E-2;
epsrel = 1E-2;
t = clock;
[X, Niter] = GLUP(H,H,rho,mu,epsabs,epsrel);
t_GLUP = etime(clock,t);

meanX = mean(X,2);
[NX_sort, idx] = sort(meanX,'descend');
tol = 0.02;
R_tol = sum(meanX>tol);
idx = idx(1:R_tol);

H2 = H(:,idx);
H2_idx = H_idx(idx,1);


%% Perform FCLS // GLUP
R_det = length(H2_idx); 
X_est = zeros(R_det, nt);
C = H2'*H2;
d = - H2'*M;
l = zeros(R_det,1);
u = ones(R_det,1);
A = ones(1,R_det);
b =1;

for n1 = 1:nt
    X_est(:, n1) = qpas(C,d(:,n1),[],[],A,b,l,u);
end


%% Evaluation metrics // GLUP 
M_est = H2*X_est;
RMSE_GLUP = Compute_RMSE(M,M_est); 
[max_ang, min_angle, avg_angle_GLUP] = Angle_Difference(M,M_est); 


%% NFINDR
t = clock ; 
[E,E_indicie] = ZMD_NFINDER(H,R_det);
t_NFNDR = etime(clock,t); 


%% Perform FCLS // NFINDR
X_est2 = zeros(R_det, nt);
C = E'*E;
d = - E'*M;
l = zeros(R_det,1);
u = ones(R_det,1);
A = ones(1,R_det);
b =1;

for n1 = 1:nt
    X_est2(:, n1) = qpas(C,d(:,n1),[],[],A,b,l,u);
end


%% Evaluation metrics // NFINDR
M_est2 = E*X_est2;
RMSE_NFINDR = Compute_RMSE(M,M_est2); 
[max_ang2, min_angle2, avg_angle_NFINDR] = Angle_Difference(M,M_est2); 


%% SDSOMP
t = clock ; 
[E3,E_indicie] = SDSOMP(H,length(H2_idx));
t_SDSOMP = etime(clock,t); 


%% Perform FCLS // SDSOMP
X_est3 = zeros(R_det, nt);
C = E3'*E3;
d = - E3'*M;
l = zeros(R_det,1);
u = ones(R_det,1);
A = ones(1,R_det);
b =1;

for n1 = 1:nt
    X_est3(:, n1) = qpas(C,d(:,n1),[],[],A,b,l,u);
end


%% Evaluation metrics // SDSOMP
M_est3 = E*X_est3;
RMSE_SDSOMP = Compute_RMSE(M,M_est3); 
[max_ang3, min_angle3, avg_angle_SDSOMP] = Angle_Difference(M,M_est3); 


%% Draw abundance map for each identified endmember

%% Using Subplot
% figure; 
% for i=1:R_det
%     % // GLUP
%     Map = X_est(i,:)';
%     Map = reshape(Map,nl,nc);
%     subplot(3,R_det,i); imshow(mat2gray(Map)); box off ; axis off; 
%     
%     % // NFINDR
%     Map = X_est2(i,:)';
%     Map = reshape(Map,nl,nc);
%     subplot(3,R_det,i+R_det); imshow(mat2gray(Map));  box off ; axis off; 
%     
%     % // SDSOMP
%     Map = X_est3(i,:)';
%     Map = reshape(Map,nl,nc);
%     subplot(3,R_det,i+2*R_det); imshow(mat2gray(Map));  box off ; axis off; 
% end

%% In seperate figures 
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

for i=1:R_det
    
    fig1= figure; 
    
    position = get(fig1,'Position');
    outerpos = get(fig1,'OuterPosition');
    borders = outerpos - position;
    
    edge = -borders(1)/2;
    
    pos1 = [edge,scnsize(4) * (2/3),outerpos(3),outerpos(4)];
    pos2 = [edge + outerpos(3) + 100 ,pos1(2),outerpos(3),outerpos(4)];
    pos3 = [edge + 2*outerpos(3) + 200,pos1(2),outerpos(3),outerpos(4)];
    
    % // GLUP
    Map = X_est(i,:)';
    Map1 = reshape(Map,nl,nc);
    set(fig1,'OuterPosition',pos1) ; imshow(mat2gray(Map1)); title('GLUP');  box off ; axis off;

    % // NFINDR
    Map = X_est2(i,:)';
    Map2 = reshape(Map,nl,nc);
    fig2 = figure;     
    set(fig2,'OuterPosition',pos2) ; imshow(mat2gray(Map2)); title('NFINDR'); box off ; axis off; 
    
    % // SDSOMP
    Map = X_est3(i,:)';
    Map3 = reshape(Map,nl,nc);
    fig3 = figure;      
    set(fig3,'OuterPosition',pos3) ; imshow(mat2gray(Map3)); title('SDSOMP'); box off ; axis off;
end

% load test5_PaviaPartGLUP.mat;
