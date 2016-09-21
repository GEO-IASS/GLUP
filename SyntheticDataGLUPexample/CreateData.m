function [] = CreateData(SNR, ne, Nt, opt)

% Create data with the given:
% SNR, number of endmembers ne, and a total number of pixels Nt
% opt places the in random columns if 1, or in begining of data if 0 or other
% Saves :
% Data, abundances A, endmembers E, and mutual coherence of endmembers eta_E, and mutual coherence of the Data eta_Data
% Filename : (Nt)pxl(ne)end(SNR)dB(opt) (ex: 100pxl3end30dB1.mat)

load endmembers.mat;
M(:,1) = []; % eliminate the first row which contains the frequencies

% randomly select ne endmembers

[~, My] = size(M); 
x = rand(1,My);
[~, indx] = sort(x);
indx = indx(1:ne);
E = M(:,indx);

% Create abundances with Dirichlet and unit parameter

A = randg(1,ne, Nt - ne);
A = A./(ones(ne,1)*sum(A));

% Inject endmembers in data 

Data = [E E*A];
if opt ==1
    tmp = round(Nt*rand(Nt,1)); 
    tmp = tmp(1:ne,1);
    cpy = Data(:,tmp); 
    Data(:,tmp)=Data(:,1:ne); 
    Data(:,1:ne)=cpy;
end

% Add Gaussian noise

sigma2 = 10^(-SNR/10)*mean(mean(Data.^2));
Data = Data + sqrt(sigma2)*randn(size(Data));

% Mutual coherence of Endmembers

eta_E = zeros(1,(ne*ne-ne)/2);
k = 1;
for i=1:ne
    for j=i+1:ne
        ai = E(:,i); aj = E(:,j);
        eta_E(1,k) =abs(ai'*aj)/(norm(ai)*norm(aj)); 
        k = k + 1; 
    end
end
eta_E = max(eta_E); 

% Mutual coherence of Data
eta_Data = zeros(1,(Nt*Nt-Nt)/2);
k = 1;
for i=1:Nt
    for j=i+1:Nt
        ai = Data(:,i); aj = Data(:,j);
        eta_Data(1,k) =abs(ai'*aj)/(norm(ai)*norm(aj)); 
        k = k + 1; 
    end
end
eta_Data = max(eta_Data); 

% save variable
if opt == 1
    file_name = strcat(sprintf('%ld',Nt),'pxl',sprintf('%ld',ne),'end',sprintf('%ld',SNR),'dB',sprintf('%ld',opt),'.mat');
else
    file_name = strcat(sprintf('%ld',Nt),'pxl',sprintf('%ld',ne),'end',sprintf('%ld',SNR),'dB',sprintf('%ld',0),'.mat');
end
    
save(file_name,'Data', 'A', 'E', 'eta_E', 'eta_Data', 'indx');
