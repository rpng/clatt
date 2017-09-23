function [zhat] = robot_measurement_model(x)

fpos = 3 + 1;

for i = 1:size(x,2)
    
    C = [cos(x(3,i)) -sin(x(3,i)); sin(x(3,i)) cos(x(3,i)) ];
    
    k_xL = C'*(x(fpos:fpos+1,i)-x(1:2,i));
    
    rho = norm(k_xL);
    th = atan2(k_xL(2),k_xL(1));
    zhat(:,i)  = [rho; th];
    
end