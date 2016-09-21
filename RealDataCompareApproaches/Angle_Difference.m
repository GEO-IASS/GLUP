function [max_ang, min_angle, avg_angle] = Angle_Difference(M,M_est)

% M = rand(2,3); M_est = rand(2,3); 
[p,n] = size(M); 

temp = sqrt(sum(M.^2)); temp = repmat(temp,p,1); 
M = M./temp ; 

temp = sqrt(sum(M_est.^2)); temp = repmat(temp,p,1); 
M_est = M_est./temp ;

temp = M.*M_est;
temp = sum(temp); 
ang = acos(temp)*180/pi; 

max_ang = max(ang);
min_angle = min(ang);
avg_angle = sum(ang)/n; 
