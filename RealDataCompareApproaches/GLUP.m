function [X,Niter] = GLUP(H,Y,rho,mu,epsprim,epsdual)

% ADMM algorithm for Group Lasso with Unit sum and Positivity constraints
%
%     variable X
%     minimize( norm(Y-AX,'fro')^2 + mu * sum(norm(X,2,2)^1) );
%     subject to
%         X >= 0;
%         1 * X = 1;

[lH,cH] = size(H);
[lY,cY] = size(Y);
if lY ~= lH;error('taille de H er Y'); end

% init and precompute matrices

A = [eye(cH,cH);ones(1,cH)];
B = [-eye(cH,cH);zeros(1,cH)];
C = [zeros(cY,cY);ones(1,cY)];

H2A2 = H'*H + rho*(A'*A);
HY = H'*Y;
AB = A'*B;
AC = A'*C;

X = zeros(cH,cY);
Z = zeros(cH,cY);

% update stopping criteria parameters

Lambdam = zeros(cH+1,cY);
Lambda = zeros(cH+1,cY);
Zm = zeros(cH,cY);

Niter = 0;
err1 = inf;
err2 = inf;

while (err1 > epsprim || err2 > epsdual) && Niter < 2000;
    Niter = Niter + 1;
    X = H2A2\(HY - rho*AB*Zm + AC - A'*Lambdam);
    
    for k = 1:cH;
        
        v = max(X(k,:) + (1/rho) * Lambda(k,:),0);
        
        if norm(v) == 0;
            Z(k,:) = zeros(1,cY);
        else
            Z(k,:) = max(1 - (mu/rho)/norm(v),0) * v;
        end
        
    end
 
    R = A*X + B*Z - C;
    P = rho * AB *(Z-Zm);
    Lambda = Lambdam +rho * R;
    
    Lambdam = Lambda;
    Zm = Z;
    
    err1 = norm(R); 
    err2 = norm(P); 
end

