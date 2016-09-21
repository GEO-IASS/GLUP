function [Xs, indx_] = SDSOMP(X,N)
% Reference : Greedy algorithms for pure pixels identification in hyperspectral
% unmixing : a multiple measurement view point
% Input : X MMV matrix (column wise data )- N Number of endmembers
% Output : Xs matrix of endmembers - indx_ index of selected endmembers (columns) in X 

[M, L] = size(X);
Xs = zeros(M,N);
indx_ = zeros(1, N);
R = X;

for k = 1:N
    corr = sum(R'*X);
    [~, indx] = max(corr);
    
    Xs(:,k) = X(:,indx);
    indx_ (1,k) = indx;
    
    R = X - Xs*pinv(Xs)*X;
end
