%% Real-world simulation
function [v_m,omega_m,v,omega,xR_true,zr,Rr, zl,Rl,  xT_true,PHI,Qd,zt,Rt ] = rws(nR,nSteps, nL,xL_true, dt,  ...
    v_true,omega_true,sigma_v,sigma_w, sigma_r,sigma_th,sigma_p, ...
    nT, vt, sigma_a, at, sigma_j,dim_target, ...
    max_range,min_range, r_max,omega_max,DORANDOM,SIGPERCENT)


%% robots' starting poses
x0 = zeros(3,nR);
xyinit = 0;
if nR>0, x0(:,1) = [xyinit, xyinit, 0]; end
if nR>1, x0(:,2) = [xyinit, -xyinit, 0]; end
if nR>2, x0(:,3) = [-xyinit, -xyinit, 0]; end
if nR>3, x0(:,4) = [-xyinit, xyinit, 0]; end


%% generate robot odometry
for ell = 1:nR
    
    %true velocities w/ noise at each nSteps
    v(ell,:) = v_true(ell,:);
    omega(ell,:) = omega_true(ell,:);
    
    for k = 1:nSteps
        if k==1
            xR_true(:,ell,k) = x0(:,ell);
        else
            xR_true(1,ell,k) = xR_true(1,ell,k-1)+v(ell,k-1)*dt*cos(xR_true(3,ell,k-1));
            xR_true(2,ell,k) = xR_true(2,ell,k-1)+v(ell,k-1)*dt*sin(xR_true(3,ell,k-1));
            xR_true(3,ell,k) = pi_to_pi(xR_true(3,ell,k-1)+omega(ell,k-1)*dt);
            
            if DORANDOM
                max_change = omega_max*dt;
                
                if norm(xR_true(1:2,ell,k),2)>0.98*r_max
                    %then turn back to the center:
                    %the direction towards the center:
                    phi_c=atan2(-xR_true(2,ell,k),-xR_true(1,ell,k));
                    
                    %we turn as much as the robot can towards the center:
                    phi1=atan2(sin(xR_true(3,ell,k-1)),cos(xR_true(3,ell,k-1)));
                    
                    if abs(phi1-phi_c)<max_change
                        xR_true(3,ell,k)=phi_c;
                    else
                        xR_true(3,ell,k)=xR_true(3,ell,k)-max_change*sign(phi1-phi_c);
                        
                        %crude handling of a special case...
                        if ((xR_true(1,ell,k)>0.98*r_max) && abs(phi_c)>pi/10)
                            xR_true(3,ell,k)=phi_c;
                        end
                    end
                    
                    omega(ell,k-1) = pi_to_pi(xR_true(3,ell,k)-xR_true(3,ell,k-1))/dt;
                end
            end
            
        end
    end %k
    
    %noisy odom measurements
    nv = sigma_v*randn(size(v(ell,:)));
    v_m(ell,:) = v(ell,:) + nv;
    nw = sigma_w*randn(size(omega(ell,:)));
    omega_m(ell,:) = omega(ell,:) + nw;
    
end%ell



%% target kinematic model: driven by continuous-time white noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % **constant velocity model / zero acceleration model**
% %  see Bar-Shalom pp.269
if dim_target ==4
    Q = zeros(2,2,nT,nSteps); %uncertainty of acceleration
    for i = 1:nT
        for k = 1:nSteps
            Q(:,:,i,k)  = sigma_a^2*eye(2); % Q can also used for time-varing case
        end
    end
    
    Qd = zeros(4,4, nT,nSteps); % The Discrete State Propogate Uncertainty
    for i=1:nT
        for k = 1:nSteps
            Qd(:,:,i,k) = kron( [  1/3*(dt^3) ,   1/2*(dt^2);    1/2*(dt^2) ,   (dt) ]  ,  Q(:,:,i,k) ); %\int [t;1]'*Q*[t;1] dt
        end
    end
    
    % Generate the Real State of Target
    xT_true = zeros(4,nT,nSteps);
    for i = 1:nT
        ax = zeros(1,nSteps-1);
        ay = zeros(1,nSteps-1);
        for k = 1:nSteps-1
            a = sqrtm(Q(:,:,k))*randn(2,1) ;
            ax(k) = a(1);
            ay(k) = a(2);
        end
        
        PT_init = 1e0*eye(4); %uncertainty of target's initial state
        xT_init_true = [10;-10;-vt;vt]; % initial state of targets
        xT_true(:,i,1) = xT_init_true * (-1)^0; % mvnrnd(xT_init_true, PT_init)';
        for k=2:nSteps
            xT_true(1,i,k) = xT_true(1,i,k-1)+xT_true(3,i,k-1)*(dt)  +1/2*ax(k-1)*(dt^2);
            xT_true(2,i,k) = xT_true(2,i,k-1)+xT_true(4,i,k-1)*(dt)  +1/2*ay(k-1)*(dt^2);
            xT_true(3,i,k) = xT_true(3,i,k-1)+ax(k-1)*(dt);
            xT_true(4,i,k) = xT_true(4,i,k-1)+ay(k-1)*(dt);
        end
    end
    
    % The State Transition Matrix
    PHI = zeros(4,4,nT,nSteps);
    FF = [1 dt; 0 1];
    perPhi = kron(FF, eye(2));
    for i=1:nT
        for k=1:nSteps
            PHI(:,:,i,k) = perPhi;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % **constant acceleration model**
    % % see Bar-Shalom pp.271
elseif dim_target==6
    Qk = sigma_j^2*eye(2); % Q can also used for time-varing case
    QQ = [ 1/20*dt^5, 1/8*dt^4, 1/6*dt^3 ;
        1/8*dt^4, 1/3*dt^3, 1/2*dt^2 ;
        1/6*dt^3, 1/2*dt^2, dt ];
    Qd = zeros(6,6, nT,nSteps); % The Discrete State Propogate Uncertainty
    for i=1:nT
        for k = 1:nSteps
            Qd(:,:,i,k) = kron(QQ, Qk);
        end
    end
    
    % Generate the Real State of Target
    xT_true = zeros(6,nT,nSteps);
    for i = 1:nT
        jx = zeros(1,nSteps-1);
        jy = zeros(1,nSteps-1);
        for k = 1:nSteps-1
            jerk = sqrtm(Qk)*randn(2,1) ;
            jx(k) = jerk(1);
            jy(k) = jerk(2);
        end
        
        PT_init = 1e0*eye(6); %uncertainty of target's initial state
        xT_init_true = [10;-10;-vt;vt; -at; at]; % initial state of targets
        xT_true(:,i,1) = xT_init_true * (-1)^0; % mvnrnd(xT_init_true, PT_init)';
        for k=2:nSteps
            xT_true(1,i,k) = xT_true(1,i,k-1)+xT_true(3,i,k-1)*(dt) +1/2*xT_true(5,i,k-1)*(dt^2)   +1/6*jx(k-1)*dt^3;
            xT_true(2,i,k) = xT_true(2,i,k-1)+xT_true(4,i,k-1)*(dt) +1/2*xT_true(6,i,k-1)*(dt^2) +1/6*jy(k-1)*dt^3;
            xT_true(3,i,k) = xT_true(3,i,k-1)+xT_true(5,i,k-1)*(dt)  +1/2*jx(k-1)*dt^2;
            xT_true(4,i,k) = xT_true(4,i,k-1)+xT_true(5,i,k-1)*(dt) +1/2*jy(k-1)*dt^2;
            xT_true(5,i,k) = xT_true(5,i,k-1)+jx(k-1)*dt;
            xT_true(6,i,k) = xT_true(6,i,k-1)+jy(k-1)*dt;
        end
    end
    
    % The State Transition Matrix
    FF = [1 dt 1/2*dt^2 ; 0 1 dt; 0 0 1];
    perPhi = kron(FF, eye(2));
    PHI = zeros(6,6,nT,nSteps);
    for i=1:nT
        for k=1:nSteps
            PHI(:,:,i,k) = perPhi;
        end
    end
end



%% generate measurements
% %Note that only first curr_meas_num nonzeros are actual measuremetns
zr = zeros(nR,3, nR, nSteps); %robot-to-robot
zl = zeros(nR,3, nL, nSteps); %robot-to-landmark
zt = zeros(nR,3, nT, nSteps); %robot-to-target

for ell = 1:nR
    for k = 1:nSteps
        Rr{ell,k} = [];
        Rl{ell,k} = [];
        Rt{ell,k} = [];
    end
end


for ell = 1:nR
    for k = 1:nSteps
        
        
        %% robot-to-robot measurements
        
        curr_meas_num = 0;
        
        for j = 1:nR
            
            if j==ell, continue, end
            
            % measurement wrt *global* frame
            [th,r] = cart2pol(xR_true(1,j,k)-xR_true(1,ell,k), xR_true(2,j,k)-xR_true(2,ell,k));
            %measurement wrt robot
            th = pi_to_pi(th-xR_true(3,ell,k));
            
            if SIGPERCENT, sigma_r = sigma_p*r; end %sigma_p percentage of range
            
            %use measurement only if landmark is closer than max_range
            if r<max_range && r>min_range
                
                curr_meas_num = curr_meas_num+1;
                
                %distance-bearing meausement
                Rii = diag([sigma_r^2, sigma_th^2]);
                Rr{ell,k} = blkdiag(Rr{ell,k},Rii);
                %noise = mvnrnd([0;0],Rii);
                r = r + sigma_r*randn; %noise(1);%
                th = th + sigma_th*randn; %noise(2);%
                
                %store measurement, and landmark id
                zr(ell,1:2,curr_meas_num,k) = [r;th]; %[dx;dy];
                zr(ell,3,curr_meas_num,k) = j;
                
            end
            
        end%j=nR
        
        
        
        
        %% robot-to-landmark measurements
        
        curr_meas_num = 0;
        
        for j = 1:nL
            
            % measurement wrt *global* frame
            [th,r] = cart2pol(xL_true(1,j)-xR_true(1,ell,k), xL_true(2,j)-xR_true(2,ell,k));
            %measurement wrt robot
            th = pi_to_pi(th-xR_true(3,ell,k));
            
            if SIGPERCENT, sigma_r = sigma_p*r; end  %sigma_p percentage of range
            
            %use measurement only if landmark is closer than max_range
            if r<max_range && r>min_range
                
                curr_meas_num = curr_meas_num+1;
                
                %distance-bearing meausement
                Rii = diag([sigma_r^2,sigma_th^2]);
                Rl{ell,k} = blkdiag(Rl{ell,k},Rii);
                %noise = mvnrnd([0;0],Rii);
                r = r + sigma_r*randn; %noise(1);%
                th = th + sigma_th*randn; %noise(2);%
                                
                %store measurement, and landmark id
                zl(ell,1:2, curr_meas_num,k) = [r;th];%[dx;dy];%
                zl(ell,3, curr_meas_num,k) = j;
                
            end
            
        end%nL
        
        
        %% robot-to-target measurements
        
        curr_meas_num = 0;
        
        for j = 1:nT
            
            % measurement wrt *global* frame
            [th,r] = cart2pol(xT_true(1,j,k)-xR_true(1,ell,k), xT_true(2,j,k)-xR_true(2,ell,k));
            %measurement wrt robot
            th = pi_to_pi(th-xR_true(3,ell,k));
            
            if SIGPERCENT, sigma_r = sigma_p*r; end  %sigma_p percentage of range
            
            if 1 %r<max_range && r>min_range  %always measures targets!!!
                
                curr_meas_num = curr_meas_num+1;
                
                %distance-bearing meausement
                Rii = diag([sigma_r^2,sigma_th^2]);
                Rt{ell,k} = blkdiag(Rt{ell,k},Rii);
                %noise = mvnrnd([0;0],Rii);
                r = r + sigma_r*randn; %noise(1);%
                th = th + sigma_th*randn; %noise(2);%
                
                %store measurement, and landmark id
                zt(ell,1:2, curr_meas_num,k) = [r;th];%[dx;dy];%
                zt(ell,3, curr_meas_num,k) = j;
                
            end
            
        end%nT
        
        
    end%nSteps
    
end%ell=nR

