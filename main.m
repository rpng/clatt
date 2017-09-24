%% -----------------------------------------------------------------------------------------------
% %  Simulation script: 2D Cooperative Localization and Target Tracking
% %  Reference:
% %  - Guoquan Huang, Michael Kaess and John Leonard,
% %    Consistent Unscented Incremental Smoothing for Multi-robot Cooperative Target Tracking.
% %    Robotics and Autonomous Systems, vol. 69, pp. 52-67, July 2015.
% %  Notes:
% %  - SuiteSparse is required: http://faculty.cse.tamu.edu/davis/suitesparse.html
% %
% %  Copyright (C) Robot Perception and Navigation Group (RPNG) - All Rights Reserved
% %  Unauthorized copying of this file, via any medium is strictly prohibited
% %  Proprietary and confidential material are included.
% %
% %  Written by:
% %  Guoquan (Paul) Huang <ghuang@cail.mit.edu>
%% -----------------------------------------------------------------------------------------------


clear all
close all
clc


%%
randomselect='repeatlast'; % 'random' or 'repeatlast';
switch randomselect
    case 'repeatlast'
        load randstates;
        rand('state',randstate);
        randn('state',randnstate);
    case 'random'
        % do nothing
    otherwise
        rand('state',randomselect);
        randn('state',randomselect);
end
% %
randstate=rand('state');
randnstate=randn('state');
save randstates.mat randstate randnstate;


%% simulation parameters %%
dt = 1; %sample time rate
v_true = .5; % robot linear velocity
omega_true = 0.05; %robot angular velocity, 0 - line; otherwise - circle
DORANDOM = 1; %random traj or not
sigma = .02*v_true;
sigma_v = sigma/sqrt(2);
sigma_w = 2*sqrt(2)*sigma;

% % measurement parameters
global gRB
gRB = 2; %1 - range-only; 2 - bearing-only

SIGPERCENT = 1; %sigma_r is  percentagae of distance or not
sigma_p = .03; %noise is the percentagae of distance- measurement
sigma_r = .1; %range measmnt noise
sigma_th = 3*pi/180;%bearing measuremnt noise

% % target stochastic motion model: constant velocity or constant acceleration
dim_target = 6; %4 - constant velocity; 6 - constant accel.
vt = .1; %target velocity
sigma_a = .01; %acceleration noise in target motion
at = 1e-6; %target accel.
sigma_j = 2e-4; %jerk noise

nR = 6; %number of robots
nT = 3; %number of targets
nL = 0; %number of landmarks, if nonzero, this is CSLAM-TT - Not used for this work.
if nL>0, error('Currently we do not consider landmarks yet!'), end

nSteps = 120; %nubmer of time steps
nRuns = 1; %number of monte carlo runs

max_range = inf; %if inf (no sensing range), then measure everything all the time
min_range = 0;
env_size = 100;
if DORANDOM
    omega_max = .25; %.5;%used in random motion
    radius_max = 25;%the max radius of confined arena
else
    omega_max = .05; %.5;%used in random motion
    radius_max = 50;%the max radius of confined arena
end
n_init_steps = 0;


%% allocate memory for saving estimation results %%

% % naive UIS: naive unscented incremental smoothing
xRest_isam2 = zeros(3,nR,nSteps,nRuns); %estimated traj
xRerr_isam2 = zeros(3,nR,nSteps,nRuns); %all err state
Prr_isam2 = zeros(3,nR,nSteps,nRuns); %actually diag of Prr
rmsRp_isam2 =  zeros(nR,nSteps,nRuns); %rms of robot position
rmsRth_isam2 = zeros(nR,nSteps,nRuns); %rms of robot orientation
xTest_isam2 = zeros(dim_target,nT,nSteps,nRuns); %estimated target trajectory
xTerr_isam2 = zeros(dim_target,nT,nSteps,nRuns); %all target error state
PTT_isam2 = zeros(dim_target,dim_target,nT,nSteps,nRuns);
rmsT_isam2 =  zeros(nT,nSteps,nRuns); %rms of target position
rmsTv_isam2 = zeros(nT,nSteps,nRuns); %rms of target velocity

% % UIS: unscented incremental smoothing
xRest_isam3 = zeros(3,nR,nSteps,nRuns); %estimated traj
xRerr_isam3 = zeros(3,nR,nSteps,nRuns); %all err state
Prr_isam3 = zeros(3,nR,nSteps,nRuns); %actually diag of Prr
rmsRp_isam3 =  zeros(nR,nSteps,nRuns); %rms of robot position
rmsRth_isam3 = zeros(nR,nSteps,nRuns); %rms of robot orientation
xTest_isam3 = zeros(dim_target,nT,nSteps,nRuns); %estimated target trajectory
xTerr_isam3 = zeros(dim_target,nT,nSteps,nRuns); %all target error state
PTT_isam3 = zeros(dim_target,dim_target,nT,nSteps,nRuns);
rmsT_isam3 =  zeros(nT,nSteps,nRuns); %rms of target position
rmsTv_isam3 = zeros(nT,nSteps,nRuns); %rms of target velocity

% % batch-MAP
xRest_bmap = zeros(3,nR,nSteps,nRuns); %estimated traj
xRerr_bmap = zeros(3,nR,nSteps,nRuns); %all err state
Prr_bmap = zeros(3,nR,nSteps,nRuns); %actually diag of Prr
rmsRp_bmap =  zeros(nR,nSteps,nRuns); %rms of robot position
rmsRth_bmap = zeros(nR,nSteps,nRuns); %rms of robot orientation
xTest_bmap = zeros(dim_target,nT,nRuns);
xTerr_bmap = zeros(dim_target,nT,nRuns);
PTT_bmap = zeros(dim_target,nT,nRuns);
rmsT_bmap = zeros(nT,nRuns);
rmsTv_bmap = zeros(nT,nRuns);



%% generate landmarks (same in each run) %%
if nL==1
    xL_true_fixed = [0; v_true/omega_true]; %put the landmark at the center
else
    if ~DORANDOM
        xL_true_fixed = gen_map(nL,v_true,omega_true,min_range, max_range, nSteps,dt);
    else
        xL_true_fixed = [rand(2,nL)  -  0.5*ones(2,nL)] * env_size ;
    end
end


%% generate true odometry (noise-free) for all time steps %%
for i = 1:nR
    vtemp(i,:) = v_true*ones(1,nSteps);
    if ~DORANDOM
        omegatemp(i,:) = omega_true*ones(1,nSteps) ;
    else
        omegatemp(i,:) = ((rand(1,nSteps)  - .5))  *1;  %random
    end
end
v_true = vtemp;
omega_true = omegatemp;


%% Monte Carlo simulations %%
for kk = 1:nRuns
    
    kk
    
    %% generate real world simulation data % %
    xL_true(:,:,kk) = xL_true_fixed;
    [v_m,omega_m, v_true_all,omega_true_all, xR_true(:,:,:,kk), zr,Rr, zl,Rl,  xT_true(:,:,:,kk),PHIT, QTd,zt,Rt] = rws(nR,nSteps,nL, xL_true(:,:,kk), dt,v_true,omega_true,sigma_v,sigma_w,sigma_r,sigma_th,sigma_p,  nT,vt,sigma_a, at,sigma_j,dim_target,  max_range,min_range, radius_max,omega_max, DORANDOM,SIGPERCENT);
    
    
    %%  initialization % %
    % robot
    xRtrue0 = reshape( xR_true(:,:,1,kk), 3*nR,1);
    PR0 = eye(length(xRtrue0)) *1e-10;
    xR0 = xRtrue0 ;
    k = n_init_steps;
    
    % target
    xTtrue0 = reshape(xT_true(:,:,1,kk),dim_target*nT,1) ;
    if dim_target ==4
        PT0 = kron(eye(nT), blkdiag(eye(2)*1e0, eye(2)*1e-2)) ;
    elseif dim_target==6
        PT0 = kron(eye(nT),  blkdiag( blkdiag(eye(2)*1e0, eye(2)*1e-2) , eye(2)*1e-6 ) );
    else
        error('has to be constant velocity or constant acceleration!')
    end
    xT0 = mvnrnd(xTtrue0,PT0)';
    
    %initial state: augmenting robot and target
    x0 = [xR0; xT0];
    P0 = blkdiag(PR0, PT0);
    
    % naive UIS
    xe_isam2 = x0;
    Pe_isam2 = P0;
    xlin_isam2 = xe_isam2;
    Plin_isam2 = P0;
    xprior_isam2 = x0;
    Pprior_isam2 = P0;
    Rp_isam2 = chol(inv(P0),'lower');
    bp_isam2 = zeros(length(xe_isam2),1);
    qp_isam2 = 1:length(xe_isam2);
    %
    for i = 1:nR
        xRest_isam2(:,i,k+1,kk) = xe_isam2(3*i-2:3*i,1);
        err = xR_true(:,i,k+1,kk) - xe_isam2(3*i-2:3*i,1);
        err(3) = pi_to_pi(err(3));
        xRerr_isam2(:,i,k+1,kk) = err;
        Prr_isam2(:,i,k+1,kk) = diag(Pe_isam2(3*i-2:3*i,3*i-2:3*i));
        rmsRp_isam2(i,k+1,kk) = err(1:2,1)'*err(1:2,1);
        rmsRth_isam2(i,k+1,kk) = err(3,1)'*err(3,1);
    end
    for i = 1:nT
        xTest_isam2(:,i,k+1,kk) = xe_isam2(3*nR+[dim_target*(i-1)+1:dim_target*i],1);
        PTT_isam2(:,:,i,k+1,kk) = PT0(dim_target*(i-1)+1:dim_target*i,dim_target*(i-1)+1:dim_target*i);
        Pk_isam2 = PT0(dim_target*(i-1)+1:dim_target*i,dim_target*(i-1)+1:dim_target*i);
        xk_isam2 = xe_isam2(3*nR+[dim_target*(i-1)+1:dim_target*i],1);
        err = xT_true(:,i,k+1,kk)-xk_isam2;
        xTerr_isam2(:,i,k+1,kk) = err;
        rmsT_isam2(i,k+1,kk) = err(1:2,1)'*err(1:2,1);
        rmsTv_isam2(i,k+1,kk) = err(3:4,1)'*err(3:4,1);
    end
    
    % UIS
    xe_isam3 = x0;
    Pe_isam3 = P0;
    xlin_isam3 = xe_isam3;
    Plin_isam3 = P0;
    xprior_isam3 = x0;
    Pprior_isam3 = P0;
    Rp_isam3 = chol(inv(P0),'lower');
    bp_isam3 = zeros(length(xe_isam3),1);
    qp_isam3 = 1:length(xe_isam3);
    %
    for i = 1:nR
        xRest_isam3(:,i,k+1,kk) = xe_isam3(3*i-2:3*i,1);
        err = xR_true(:,i,k+1,kk) - xe_isam3(3*i-2:3*i,1);
        err(3) = pi_to_pi(err(3));
        xRerr_isam3(:,i,k+1,kk) = err;
        Prr_isam3(:,i,k+1,kk) = diag(Pe_isam3(3*i-2:3*i,3*i-2:3*i));
        rmsRp_isam3(i,k+1,kk) = err(1:2,1)'*err(1:2,1);
        rmsRth_isam3(i,k+1,kk) = err(3,1)'*err(3,1);
    end
    for i = 1:nT
        xTest_isam3(:,i,k+1,kk) = xe_isam3(3*nR+[dim_target*(i-1)+1:dim_target*i],1);
        PTT_isam3(:,:,i,k+1,kk) = PT0(dim_target*(i-1)+1:dim_target*i,dim_target*(i-1)+1:dim_target*i);
        Pk_isam3 = PT0(dim_target*(i-1)+1:dim_target*i,dim_target*(i-1)+1:dim_target*i);
        xk_isam3 = xe_isam3(3*nR+[dim_target*(i-1)+1:dim_target*i],1);
        err = xT_true(:,i,k+1,kk)-xk_isam3;
        xTerr_isam3(:,i,k+1,kk) = err;
        rmsT_isam3(i,k+1,kk) = err(1:2,1)'*err(1:2,1);
        rmsTv_isam3(i,k+1,kk) = err(3:4,1)'*err(3:4,1);
    end
    
    
    
    %% sequential estimation over time %%
    for k = n_init_steps+1:nSteps-1 %first n_init_steps for ekf propagation to produce nonzero init cov
        
        timestep = k+1
        
        dofk = 3*nR+dim_target*nT;
        
        %% Naive UIS %%
        if length(xe_isam2)~=3*nR*(k)+dim_target*nT*(k)
            error('Naive UIS: dimension of state is not correct!!!')
        end
        
        DO_EXACT_COV_REC = 0; %exact or approximate covariance recovery
        % *exact* recovering covariance used in unscented transformation, though this is expensive!
        if k+1>2 && DO_EXACT_COV_REC
            Pe_isam2 = (Rp_isam2\eye(size(Rp_isam2)))*(Rp_isam2'\eye(size(Rp_isam2))); %info matrix for the states without first robots' poses which are assumed to be known!!!
            Pe_isam2 = blkdiag(Pprior_isam2(1:3*nR,1:3*nR),  Pe_isam2(qp_isam2,qp_isam2) );
        else
            % already intialized
        end
        
        A_isam2 = Rp_isam2'*Rp_isam2;
        for i = 1:nR
            indxr = dofk*(k-1) + [ 3*i-2:3*i ];
            xRi = xe_isam2(indxr,1);
            xRk1(3*i-2:3*i,1) = [ xRi(1,1) + v_m(i,k)*dt*cos(xRi(3,1));
                xRi(2,1) + v_m(i,k)*dt*sin(xRi(3,1));
                pi_to_pi(xRi(3,1) + omega_m(i,k)*dt) ];
            
            PHI = [1 0 -v_m(i,k)*dt*sin(xRi(3,1));      0 1  v_m(i,k)*dt*cos(xRi(3,1));      0 0   1];
            G = [dt*cos(xRi(3,1))   0;        dt*sin(xRi(3,1))   0;        0     dt];
            Q = [sigma_v^2  0;   0    sigma_w^2];
            % *approx* covariance recovery for submatrix of Rp_isam2 w/o reordering
            if k+1>2 && ~DO_EXACT_COV_REC
                PRRk = inv(A_isam2(indxr-3*nR,indxr-3*nR) );
                %Rtem = Rp_isam2(indxr-3*nR,indxr-3*nR);
                %PRRk = (Rtem\eye(size(Rtem))) * (Rtem'\eye(size(Rtem)));
            else
                PRRk = Pe_isam2(indxr,indxr);
            end
            PRRk1(3*i-2:3*i,3*i-2:3*i) = PHI * PRRk *PHI' + G*Q*G';
        end
        
        for i = 1:nT
            indxt = dofk*(k-1)+3*nR+ [ dim_target*(i-1)+1:dim_target*i ];
            xTi = xe_isam2(indxt,1);
            xTk1(dim_target*(i-1)+1:dim_target*i,1) =  PHIT(:,:,i,k) * xTi;
            % *approx* covariance recovery for submatrix of Rp_isam2 w/o reordering
            if k+1>2 && ~DO_EXACT_COV_REC
                PTTk = inv(A_isam2(indxt-3*nR,indxt-3*nR) );
            else
                PTTk = Pe_isam2(indxt,indxt);
            end
            PTTk1(dim_target*(i-1)+1:dim_target*i,dim_target*(i-1)+1:dim_target*i) =  PHIT(:,:,i,k) *PTTk*PHIT(:,:,i,k)' + QTd(:,:,i,k);
        end
        xe_isam2 = [xe_isam2; xRk1; xTk1 ] ; %initial estimates, augmented by propagated new initails
        xlin_isam2 =  [xlin_isam2; xRk1; xTk1 ] ; %linearization points, augmented by propagated new initails
        Ptem = blkdiag(PRRk1,PTTk1);
        Plin_isam2 = blkdiag(Plin_isam2,Ptem) ;
        
        [xe_isam2,xlin_isam2, Rp_isam2,bp_isam2, qp_isam2] = uis(xe_isam2,xlin_isam2, Plin_isam2, xprior_isam2,Pprior_isam2, Rp_isam2,bp_isam2,qp_isam2, k+1, dt, nR, v_m, omega_m, sigma_v, sigma_w, nT, dim_target,PHIT, QTd, zr, Rr, zt, Rt, zl, Rl, nL, false);
        
        % % save results
        for i = 1:nR
            xRest_isam2(:,i,k+1,kk) = xe_isam2(dofk*k+[3*i-2:3*i],1);
            err = xR_true(:,i,k+1,kk) - xe_isam2(dofk*k+[3*i-2:3*i],1);
            err(3) = pi_to_pi(err(3));
            xRerr_isam2(:,i,k+1,kk) = err;
            rmsRp_isam2(i,k+1,kk) = err(1:2,1)'*err(1:2,1);
            rmsRth_isam2(i,k+1,kk) = err(3,1)'*err(3,1);
        end
        for i = 1:nT
            xTest_isam2(:,i,k+1,kk) = xe_isam2(dofk*k+3*nR+[dim_target*(i-1)+1:dim_target*i],1);
            xk_isam2 = xe_isam2(dofk*k+3*nR+[dim_target*(i-1)+1:dim_target*i],1);
            err = xT_true(:,i,k+1,kk)-xk_isam2;
            xTerr_isam2(:,i,k+1,kk) = err;
            rmsT_isam2(i,k+1,kk) = err(1:2,1)'*err(1:2,1);
            rmsTv_isam2(i,k+1,kk) = err(3:4,1)'*err(3:4,1);
        end
        
        
        %% UIS %%
        if length(xe_isam3)~=3*nR*(k)+dim_target*nT*(k)
            error('iSAM3: dimension of state is not correct!!!')
        end
        
        DO_EXACT_COV_REC = 0; %exact or approximate covariance recovery
        % *exact* recovering covariance used in unscented transformation though this is expensive
        if k+1>2 && DO_EXACT_COV_REC
            Pe_isam3 = (Rp_isam3\eye(size(Rp_isam3)))*(Rp_isam3'\eye(size(Rp_isam3))); %info matrix for the states without first robots' poses which are assumed to be known!!!
            Pe_isam3 = blkdiag(Pprior_isam3(1:3*nR,1:3*nR),  Pe_isam3(qp_isam3,qp_isam3) );
        else
            % already intialized
        end
        
        A_isam3 = Rp_isam3'*Rp_isam3;
        for i = 1:nR
            indxr = dofk*(k-1) + [ 3*i-2:3*i ];
            xRi = xe_isam3(indxr,1);
            xRk1(3*i-2:3*i,1) = [ xRi(1,1) + v_m(i,k)*dt*cos(xRi(3,1));
                xRi(2,1) + v_m(i,k)*dt*sin(xRi(3,1));
                pi_to_pi(xRi(3,1) + omega_m(i,k)*dt) ];
            
            PHI = [1 0 -v_m(i,k)*dt*sin(xRi(3,1));      0 1  v_m(i,k)*dt*cos(xRi(3,1));      0 0   1];
            G = [dt*cos(xRi(3,1))   0;        dt*sin(xRi(3,1))   0;        0     dt];
            Q = [sigma_v^2  0;   0    sigma_w^2];
            % *approx* covariance recovery for submatrix of Rp_isam3 w/o reordering
            if k+1>2 && ~DO_EXACT_COV_REC
                PRRk = inv(A_isam3(indxr-3*nR,indxr-3*nR) );
            else
                PRRk = Pe_isam3(indxr,indxr);
            end
            PRRk1(3*i-2:3*i,3*i-2:3*i) = PHI * PRRk *PHI' + G*Q*G';
        end
        
        for i = 1:nT
            indxt = dofk*(k-1)+3*nR+ [ dim_target*(i-1)+1:dim_target*i ];
            xTi = xe_isam3(indxt,1);
            xTk1(dim_target*(i-1)+1:dim_target*i,1) =  PHIT(:,:,i,k) * xTi;
            % *approx* covariance recovery for submatrix of Rp_isam3 w/o reordering
            if k+1>2 && ~DO_EXACT_COV_REC
                PTTk = inv(A_isam3(indxt-3*nR,indxt-3*nR) );
            else
                PTTk = Pe_isam3(indxt,indxt);
            end
            PTTk1(dim_target*(i-1)+1:dim_target*i,dim_target*(i-1)+1:dim_target*i) =  PHIT(:,:,i,k) *PTTk*PHIT(:,:,i,k)' + QTd(:,:,i,k);
        end
        xe_isam3 = [xe_isam3; xRk1; xTk1 ] ; %initial estimates, augmented by propagated new initails
        xlin_isam3 =  [xlin_isam3; xRk1; xTk1 ] ; %linearization points, augmented by propagated new initails
        Ptem = blkdiag(PRRk1,PTTk1);
        Plin_isam3 = blkdiag(Plin_isam3,Ptem) ;
        
        [xe_isam3,xlin_isam3, Rp_isam3,bp_isam3, qp_isam3] = uis(xe_isam3,xlin_isam3, Plin_isam3, xprior_isam3,Pprior_isam3, Rp_isam3,bp_isam3,qp_isam3, k+1, dt, nR, v_m, omega_m, sigma_v, sigma_w, nT, dim_target,PHIT, QTd, zr, Rr, zt, Rt, zl, Rl, nL, true);
        
        % % save results
        for i = 1:nR
            xRest_isam3(:,i,k+1,kk) = xe_isam3(dofk*k+[3*i-2:3*i],1);
            err = xR_true(:,i,k+1,kk) - xe_isam3(dofk*k+[3*i-2:3*i],1);
            err(3) = pi_to_pi(err(3));
            xRerr_isam3(:,i,k+1,kk) = err;
            rmsRp_isam3(i,k+1,kk) = err(1:2,1)'*err(1:2,1);
            rmsRth_isam3(i,k+1,kk) = err(3,1)'*err(3,1);
        end
        for i = 1:nT
            xTest_isam3(:,i,k+1,kk) = xe_isam3(dofk*k+3*nR+[dim_target*(i-1)+1:dim_target*i],1);
            xk_isam3 = xe_isam3(dofk*k+3*nR+[dim_target*(i-1)+1:dim_target*i],1);
            err = xT_true(:,i,k+1,kk)-xk_isam3;
            xTerr_isam3(:,i,k+1,kk) = err;
            rmsT_isam3(i,k+1,kk) = err(1:2,1)'*err(1:2,1);
            rmsTv_isam3(i,k+1,kk) = err(3:4,1)'*err(3:4,1);
        end
        
        
    end %k
    
    
    
    %% batch MAP %%
    disp('------batch MAP------')
    % generate initial guess: ideally use ground truth
    xinit_bmap(1:dofk,1) = x0;
    for k = 1:nSteps-1
        for i = 1:nR
            xinit_bmap(dofk*k+[3*i-2:3*i],1) = xR_true(:,i,k+1,kk) ;
        end
        for i = 1:nT
            xinit_bmap(dofk*k+3*nR+[dim_target*(i-1)+1:dim_target*i],1) =  xT_true(:,i,k+1,kk) ;
        end
    end
   
    xe_bmap = bmap(xinit_bmap, x0,P0, k+1, dt, nR, v_m, omega_m, sigma_v, sigma_w, nT, dim_target,PHIT, QTd, zr, Rr, zt, Rt, zl, Rl, nL);
    
    % save batch results
    for k = 1:nSteps-1
        for i = 1:nR
            xRest_bmap(:,i,k+1,kk) = xe_bmap(dofk*k+[3*i-2:3*i],1);
            err = xR_true(:,i,k+1,kk) - xe_bmap(dofk*k+[3*i-2:3*i],1);
            err(3) = pi_to_pi(err(3));
            xRerr_bmap(:,i,k+1,kk) = err;
            rmsRp_bmap(i,k+1,kk) = err(1:2,1)'*err(1:2,1);
            rmsRth_bmap(i,k+1,kk) = err(3,1)'*err(3,1);
        end
        for i = 1:nT
            xTest_bmap(:,i,k+1,kk) = xe_bmap(dofk*k+3*nR+[dim_target*(i-1)+1:dim_target*i],1);
            xk_bmap = xe_bmap(dofk*k+3*nR+[dim_target*(i-1)+1:dim_target*i],1);
            err = xT_true(:,i,k+1,kk)-xk_bmap;
            xTerr_bmap(:,i,k+1,kk) = err;
            rmsT_bmap(i,k+1,kk) = err(1:2,1)'*err(1:2,1);
            rmsTv_bmap(i,k+1,kk) = err(3:4,1)'*err(3:4,1);
        end
    end
    
end%kk, end of monte carlo runs



%% % average results over all runs
% % naive UIS
rmsRp_avg_isam2 = sqrt(sum(rmsRp_isam2,3)/nRuns);
rmsRth_avg_isam2 = sqrt(sum(rmsRth_isam2,3)/nRuns);
rmsT_avg_isam2 = sqrt(sum(rmsT_isam2,3)/nRuns);
rmsTv_avg_isam2 = sqrt(sum(rmsTv_isam2,3)/nRuns);
xRerr_avg_isam2 = sum(xRerr_isam2,4)/nRuns;

% % UIS
rmsRp_avg_isam3 = sqrt(sum(rmsRp_isam3,3)/nRuns);
rmsRth_avg_isam3 = sqrt(sum(rmsRth_isam3,3)/nRuns);
rmsT_avg_isam3 = sqrt(sum(rmsT_isam3,3)/nRuns);
rmsTv_avg_isam3 = sqrt(sum(rmsTv_isam3,3)/nRuns);
xRerr_avg_isam3 = sum(xRerr_isam3,4)/nRuns;

% % batch MAP
rmsRp_avg_bmap = sqrt(sum(rmsRp_bmap,3)/nRuns);
rmsRth_avg_bmap = sqrt(sum(rmsRth_bmap,3)/nRuns);
rmsT_avg_bmap = sqrt(sum(rmsT_bmap,3)/nRuns);
rmsTv_avg_bmap = sqrt(sum(rmsTv_bmap,3)/nRuns);
xRerr_avg_bmap = sum(xRerr_bmap,4)/nRuns;


%% % plotting figures
plot_all

