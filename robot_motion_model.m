function [y] = robot_motion_model(x,dt)

for i = 1:size(x,2)
    
    v = x(end-1,i);
    omega = x(end,i);
    
    y(:,i) = x(1:end-2,i);

    y(1:3,i) = [ x(1,i) + v*dt*cos(x(3,i));
        x(2,i) + v*dt*sin(x(3,i));
        pi_to_pi(x(3,i) + omega*dt) ];

end
