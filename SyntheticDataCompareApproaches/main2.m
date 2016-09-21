%% Monte Carlo for Comparing the three approaches
%% GLUP - NFINDR - SDSOMP

clear all ; close all ; clc

N_total = 100; % Nbr of repetitions for MonteCarlo

%% Create synthetic data

R = 8; % nbr of endmembers
N = 192; % nbr of mixed pixels
Nt = N+R; % total nbr of pixels

SNR = 20;

load endmembers.mat;
E = M(:,1:R);
A = randg(1,R, N);
A = A./(ones(R,1)*sum(A));
S1 = [E E*A];

%% Parameters for GLUP 
mu = 100; % regularisation parameter for GLUP 
rho = 100; % penalty parameter ADMM GLUP 

epsabs = 1E-2; % primal tolerance GLUP 
epsrel = 1E-2; % dual tolerance GLUP 

%% Memory allocation
acc1_ = zeros(1,N_total);
acc2_ = zeros(1,N_total);
acc3_ = zeros(1,N_total);

idx1_ = zeros(N_total,R); 
idx2_ = zeros(N_total,R);
idx3_ = zeros(N_total,R);

t1_total = 0;
t2_total = 0;
t3_total = 0;

for i =1:N_total
    disp(i);
    
    %% Add Gaussian noise
    sigma2 = 10^(-SNR/10)*mean(mean(S1.^2));
    S = S1 + sqrt(sigma2)*randn(size(S1));
         
    %% GLUP
    t = clock;
    [X, Niter] = GLUP(S,S,rho,mu,epsabs,epsrel);
    t = etime(clock,t);
    t1_total = t1_total + t ;
    
    tmp = mean(X,2); % figure ; bar(tmp); 
    [tmp_sort, E1_indx] = sort(-tmp);
    
    acc1_(i) = length(intersect(E1_indx(1:R),1:R))/R;
    idx1_(i,:) = E1_indx(1:R);
    
    %% NFINDR
    t = clock ; 
    [~,E2_indx] = ZMD_NFINDER(S,R);
    t = etime(clock,t);
    t2_total = t2_total + t ; 
    acc2_(i) = length(intersect(E2_indx,1:R))/R;
    idx2_(i,:) = E2_indx; 
    
    %% SDSOMP
    t = clock;
    [~, E3_indx] = SDSOMP(S,R);
    t = etime(clock,t);
    t3_total = t3_total + t ;
    acc3_(i) = length(intersect(E3_indx,1:R))/R;
    idx3_(i,:) = E3_indx;
    
end

t1_total = t1_total/N_total;
t2_total = t2_total/N_total;
t3_total = t3_total/N_total;

acc1 = sum(acc1_)/N_total;
acc2 = sum(acc2_)/N_total;
acc3 = sum(acc3_)/N_total;

%% write in file --- G = GLUP, N = NFINDR, S = SDSOMP

filename = strcat('GNS_test1.txt');
fid = fopen(filename,'a');

fprintf(fid,'(N_total, R, N, Nt, SNR, rho, mu, epsabs, epsrel) = ( %f, %f, %f, %f, %f, %f, %f, %f, %f )',N_total, R, N, Nt, SNR, rho, mu, epsabs, epsrel);fprintf(fid,'\n');

fprintf(fid,'acc(G,N,S) = (%f,%f,%f)',acc1,acc2,acc3);fprintf(fid,'\n');

fprintf(fid,'time(G,N,S) = (%f,%f,%f)',t1_total,t2_total,t3_total);fprintf(fid,'\n');

fprintf(fid,'\n'); fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');fprintf(fid,'\n');

fclose(fid);
