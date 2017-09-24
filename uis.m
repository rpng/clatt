function [x, xlin, Rp,bp,qp] = uis(x,xlin, Plin, xprior,Pprior, Rp,bp,qp, curr_step, dt, nR, v_m, omega_m, sigma_v, sigma_w, nT, dim_target, PHIT, QTd, zr, Rr,  zt, Rt, zl, Rl, nL, doOC)
% % unscented incremental smoothing (UIS) for CLATT (cooperative localization and target tracking)

%% parameters/constants
global gRB
J = [0 -1; 1 0];

% during batch update: do unscented transform or not?
% it should be NOT, since 1) it is expensive to recover covariance 2) batch MAP is consistent (relinearize anyway)
noUT = 1;

epsilon = 1e-4;
itermax = 10; %max number of iterations
nrelin = 10; %frequency of relinearization

nddt = 3;
ddt = dt/nddt;
Q = [ sigma_v^2 0; 0 sigma_w^2 ]; %odometry noise

nsteps = curr_step; %number of total time steps

dofR = 3*nR;
dofT = dim_target*nT;
dofk = dofR+dofT;

% % first robot state is known, so excluded from the state, while first target state has prior!
% state = [robots-targets; robots-targets;...]
x0 = x(1:dofR,1);
xlin0 = xlin(1:dofR,1);
x = x(dofR+1:end,1);
xlin = xlin(dofR+1:end,1);
Plin = Plin(dofR+1:end, dofR+1:end);


%% Gauss-Newton
if mod(curr_step,nrelin) && curr_step>2 %incremental QR update
    
    % NOTE: run batch update for the first step in order to easily get Rp, bp, qp by QR (though not necessary)
    
    ndim = dofT + dofk*(nsteps-1); %first robot pose is known, hence not included in the state
    if gRB==3
        nmeas = dofT + 3*nR*(nsteps-1)+dim_target*nT*(nsteps-1)  + 2*nR*nR* (nsteps-1) + 2*nR*nT*(nsteps-1);  %may be converative!!!
    else
        nmeas = dofT + 3*nR*(nsteps-1)+dim_target*nT*(nsteps-1)  + nR*nR* (nsteps-1) + nR*nT*(nsteps-1);  %may be converative!!!
    end
    
    A = sparse(nmeas,ndim); %jacobian
    b = zeros(nmeas,1); %residual
    
    k = curr_step; %only compute jacobians for current time step
    
    % % proprioception measurements % %
    % robot propagation
    for ell = 1:nR
        if k==2 %this actuall is the first robot in the state
            x1k = xlin0(3*ell-2:3*ell,1); %known xlin0 or x0
            P1k = zeros(3);
            for i=1:nddt
                xk = [ x1k(1) + v_m(ell,k-1)*ddt*cos(x1k(3));
                    x1k(2) + v_m(ell,k-1)*ddt*sin(x1k(3));
                    pi_to_pi(x1k(3) + omega_m(ell,k-1)*ddt) ];
                F1k = [ 1 0 -v_m(ell,k-1)*ddt*sin(x1k(3));
                    0 1  v_m(ell,k-1)*ddt*cos(x1k(3));
                    0 0   1 ];
                Gk = [ ddt*cos(x1k(3))   0;
                    ddt*sin(x1k(3))   0;
                    0     ddt ];
                P1k = F1k*P1k*F1k'+Gk*Q*Gk';
                x1k = xk;
            end
            Qk = P1k;
            L = chol(Qk,'lower');
            invQk = (L'\eye(size(L)))*(L\eye(size(L)));
            sqrtQk = sqrtm(invQk);
            
            indr2 = dofT + [3*ell-2:3*ell];
            indxm =dofT + [3*ell-2:3*ell ];
            
            Fk =  - eye(3);
            
            ak = xlin(indr2,1) - xk;
            ak(3) = pi_to_pi(ak(3));
            
            b(indxm,1) = sqrtQk'*ak ;
            A(indxm,indr2) = sqrtQk'*Fk ;
            
        elseif k>2
            indr1 =  dofT + dofk*(k-3)+ [3*(ell)-2 : 3*(ell)] ;  %note the first pose (k=1) is excluded in x, i.e., x_k actually corresponds to xinit_k+1
            indr2 =  dofT + dofk*(k-2)+[3*(ell)-2 : 3*(ell)] ;
            indxm =  dofT + dofk*(k-2)+[ 3*ell-2:3*ell ];
            
            [F1k,sqrtQk,xk] = prop_jacob_ut(@robot_motion_model, xlin(indr1,1),  xlin([indr1, indr2 ],1), Plin(indr1,indr1),v_m(ell,k-1),omega_m(ell,k-1),Q,dt, doOC);
            
            Fk = - eye(3);
            
            ak = xlin(indr2,1) - xk;
            ak(3) = pi_to_pi(ak(3));
            
            b(indxm,1) = sqrtQk'*ak ;
            A(indxm,indr1) = sqrtQk'*F1k ;
            A(indxm,indr2) = sqrtQk'*Fk ;
        end
    end
    
    
    % target propagation: since it is linear, unscented transform is same as the analytical jacobian
    for ell = 1:nT
        L = chol(QTd(:,:,ell,k-1),'lower');
        invQk = (L'\eye(size(L)))*(L\eye(size(L)));
        sqrtQk = sqrtm(invQk);
        
        indt1 =  dofT + dofk*(k-3)+ 3*nR +[dim_target*(ell-1)+1:dim_target*ell];
        indt2 =  dofT + dofk*(k-2)+ 3*nR +[dim_target*(ell-1)+1:dim_target*ell];
        indxm =  dofT + dofk*(k-2)+ 3*nR +[dim_target*(ell-1)+1:dim_target*ell];
        
        F1k = PHIT(:,:,ell,k-1);
        Gk = - eye(dim_target);
        ak = xlin(indt2,1) - F1k*x(indt1,1);
        
        b(indxm,1) = sqrtQk'*ak ;
        A(indxm, indt1) = sqrtQk'*F1k ;
        A(indxm, indt2) = sqrtQk'*Gk ;
    end
    
    
    % % exterocpetive measurements % %
    
    % robot-to-robot
    imeas = 0;
    for ell = 1:nR
        for i = 1:size(zr(ell,:,:,k),3)
            if zr(ell,3,i,k)<=0, continue, end
            
            Ri = Rr{ell,k}(2*i-1:2*i, 2*i-1:2*i);
            sqrtRk = sqrtm(inv(Ri));
            
            imeas = imeas+1;
            indr1 =  dofT + dofk*(k-2)+[3*ell-2:3*ell];
            indr2 =  dofT + dofk*(k-2)+[3*zr(ell,3,i,k)-2 : 3*zr(ell,3,i,k)];
            
            if gRB==3
                indxm =  dofT + dofk*(nsteps-1) + 2*nR*nR*(k-2) + [ 2*imeas-1:2*imeas ] ;
            elseif gRB==2
                indxm =  dofT + dofk*(nsteps-1) + nR*nR*(k-2) + imeas;
                sqrtRk = sqrtRk(2,2);
            elseif gRB==1
                indxm =  dofT + dofk*(nsteps-1) + nR*nR*(k-2) + imeas;
                sqrtRk = sqrtRk(1,1);
            end
            
            [r,H1,H2] = meas_jacob_ut(@robot_measurement_model, xlin([indr1,indr2],1),  xlin([indr1,indr2],1),Plin([indr1,indr2],[indr1,indr2]),zr(ell,1:2,i,k),Ri,  doOC);
            
            A(indxm, [indr1, indr2])  = [sqrtRk'*H1 , sqrtRk'*H2] ;
            b(indxm,1) = sqrtRk'*r;
        end
    end
    
    % robot-to-target
    imeas = 0;
    for ell = 1:nR
        for i = 1:size(zt(ell,:,:,k),3)
            if zt(ell,3,i,k)<=0, continue, end
            
            Ri = Rt{ell,k}(2*i-1:2*i, 2*i-1:2*i);
            sqrtRk = sqrtm(inv(Ri));
            
            imeas = imeas+1;
            indr1 =  dofT + dofk*(k-2)+ [3*ell-2:3*ell];
            indr2 =  dofT + dofk*(k-2)+ 3*nR+ [dim_target*(zt(ell,3,i,k)-1)+1 : dim_target*zt(ell,3,i,k)];
            
            if gRB==3
                indxm =  dofT + dofk*(nsteps-1) + 2*nR*nR*(nsteps-1) + 2*nR*nT*(k-2) + [2*imeas-1:2*imeas] ;
            elseif gRB==2
                indxm =  dofT + dofk*(nsteps-1) + nR*nR*(nsteps-1) + nR*nT*(k-2) + imeas ;
                sqrtRk = sqrtRk(2,2);
            elseif gRB==1
                indxm =  dofT + dofk*(nsteps-1) + nR*nR*(nsteps-1) + nR*nT*(k-2) + imeas ;
                sqrtRk = sqrtRk(1,1);
            end
            
            [r,H1,H2] = meas_jacob_ut(@robot_measurement_model,xlin([indr1,indr2],1), xlin([indr1,indr2],1),Plin([indr1,indr2],[indr1,indr2]),zt(ell,1:2,i,k),Ri, doOC);
            
            A(indxm, [indr1, indr2])  = [sqrtRk'*H1 , sqrtRk'*H2] ;
            b(indxm,1) = sqrtRk'*r;
        end
    end
    
    
    % robot-to-landmark
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reuse previously-computed jacobian and residual, and augment by the new ones
    qp = [qp, length(qp)+1:length(qp)+dofk];
    Ak = A(:,qp); % reordering, qp, which was obtained from batch step
    A = sparse([ Rp, zeros(length(Rp),dofk) ;  Ak ]); %match the dimension, as new robot and target states added
    b = [bp; b];
    
    % solve linear system using QR
    % using sparse QR with reordering
    [~, n] = size(A);
    [V,beta,p,R] = cs_qr(A) ;
    d = cs_qleft(V, beta, p, b) ;
    Rp = R(1:n,1:n);
    bp = d(1:n,1);
    
    delta_x =  Rp \ bp;
    delta_x(qp) = delta_x; %back to original ordering
    
    if ~doOC
        ndx2 = norm(delta_x);
    else
        ndx3 = norm(delta_x);
    end
    
    % % *linearization points, xlin* updated by deltax to be the new estimate
    x = xlin + delta_x;
    
    
else %%%batch every *nrelin* steps (including first step)%%%
    
    
    for iter = 1:itermax
        
        ndim = dofT + dofk*(nsteps-1); %first robot pose is known, hence not included in the state
        if gRB==3
            nmeas = dofT + 3*nR*(nsteps-1)+dim_target*nT*(nsteps-1)  + 2*nR*nR* (nsteps-1) + 2*nR*nT*(nsteps-1);  %may be converative!!!
        else
            nmeas = dofT + 3*nR*(nsteps-1)+dim_target*nT*(nsteps-1)  + nR*nR* (nsteps-1) + nR*nT*(nsteps-1);  %may be converative!!!
        end
        
        A = sparse(nmeas, ndim); %jacobian
        b = zeros(nmeas,1); %residual
        
        % % prior for the first target's state
        xTprior = xprior(3*nR+[1:dofT],1);
        PTprior = Pprior(3*nR+[1:dofT] , 3*nR+[1:dofT]);
        sqrtQk =  real(sqrtm(inv(PTprior)));
        ak = xTprior - x(1:dofT,1);
        Gk = eye(dofT);
        
        indt2 = 1:dofT;
        indxm = 1:dofT;
        
        b(indxm,1) = sqrtQk'*ak ;
        A(indxm, indt2) = sqrtQk'*Gk ;
        
        % % all measurements
        for k = 2:nsteps %because first robot pose is known
            
            % % proprioception measurements % %
            % robot propagation
            for ell = 1:nR
                if k==2 %this actuall is the first robot in the state
                    x1k = x0(3*ell-2:3*ell,1); %xprior(3*ell-2:3*ell,1); %x0 should be same as xprior
                    P1k = zeros(3);
                    for i=1:nddt
                        xk = [ x1k(1) + v_m(ell,k-1)*ddt*cos(x1k(3));
                            x1k(2) + v_m(ell,k-1)*ddt*sin(x1k(3));
                            pi_to_pi(x1k(3) + omega_m(ell,k-1)*ddt) ];
                        F1k = [ 1 0 -v_m(ell,k-1)*ddt*sin(x1k(3));
                            0 1  v_m(ell,k-1)*ddt*cos(x1k(3));
                            0 0   1 ];
                        Gk = [ ddt*cos(x1k(3))   0;
                            ddt*sin(x1k(3))   0;
                            0     ddt ];
                        P1k = F1k*P1k*F1k'+Gk*Q*Gk';
                        x1k = xk;
                    end
                    Qk = P1k;
                    L = chol(Qk,'lower');
                    invQk = (L'\eye(size(L)))*(L\eye(size(L)));
                    sqrtQk = sqrtm(invQk);
                    
                    indr2 = dofT + [3*ell-2:3*ell];
                    indxm =dofT + [3*ell-2:3*ell ];
                    
                    Fk =  - eye(3);
                    ak = x(indr2,1) - xk;
                    ak(3) = pi_to_pi(ak(3));
                    
                    b(indxm,1) = sqrtQk'*ak ;
                    A(indxm,indr2) = sqrtQk'*Fk ;
                    
                elseif k>2
                    indr1 =  dofT + dofk*(k-3)+ [3*(ell)-2 : 3*(ell)] ;  %note the first pose (k=1) is excluded in x, i.e., x_k actually corresponds to xinit_k+1
                    indr2 =  dofT + dofk*(k-2)+[3*(ell)-2 : 3*(ell)] ;
                    indxm =  dofT + dofk*(k-2)+[ 3*ell-2:3*ell ];
                    
                    if noUT
                        x1k = x(indr1,1);
                        P1k = zeros(3);
                        for i=1:nddt
                            xk = [ x1k(1) + v_m(ell,k-1)*ddt*cos(x1k(3));
                                x1k(2) + v_m(ell,k-1)*ddt*sin(x1k(3));
                                pi_to_pi(x1k(3) + omega_m(ell,k-1)*ddt) ]; %note the odometry actually is from time k
                            F1k = [1 0 -v_m(ell,k-1)*ddt*sin(x1k(3));
                                0 1  v_m(ell,k-1)*ddt*cos(x1k(3));
                                0 0   1];
                            Gk = [ddt*cos(x1k(3))   0;
                                ddt*sin(x1k(3))   0;
                                0     ddt];
                            P1k = F1k*P1k*F1k'+Gk*Q*Gk';
                            x1k = xk;
                        end
                        Qk = P1k;
                        L = chol(Qk,'lower');
                        invQk = (L'\eye(size(L)))*(L\eye(size(L)));
                        sqrtQk = sqrtm(invQk);
                        
                        F1k =  [ eye(2)   J*(x(indr2(1:2),1)-x(indr1(1:2),1));   zeros(1,2)  1 ];%beware of sign
                        Fk = - eye(3);
                        
                    else %unscented
                        [F1k,sqrtQk,xk] = prop_jacob_ut(@robot_motion_model, x(indr1,1),   x([indr1, indr2 ],1), Plin(indr1,indr1),v_m(ell,k-1),omega_m(ell,k-1),Q,dt, doOC);
                        Fk = - eye(3);
                    end
                    
                    ak = xlin(indr2,1) - xk;
                    ak(3) = pi_to_pi(ak(3));
                    
                    b(indxm,1) = sqrtQk'*ak ;
                    A(indxm,indr1) = sqrtQk'*F1k ;
                    A(indxm,indr2) = sqrtQk'*Fk ;
                end
            end
            
            
            % target propagation
            for ell = 1:nT
                indt1 =  dofT + dofk*(k-3)+ 3*nR +[dim_target*(ell-1)+1:dim_target*ell];
                indt2 =  dofT + dofk*(k-2)+ 3*nR +[dim_target*(ell-1)+1:dim_target*ell];
                indxm =  dofT + dofk*(k-2)+ 3*nR +[dim_target*(ell-1)+1:dim_target*ell];
                
                L = chol(QTd(:,:,ell,k-1),'lower');
                invQk = (L'\eye(size(L)))*(L\eye(size(L)));
                sqrtQk = sqrtm(invQk);
                
                F1k = PHIT(:,:,ell,k-1);
                Gk = - eye(dim_target);
                ak = x(indt2,1) - F1k*x(indt1,1);
                
                b(indxm,1) = sqrtQk'*ak ;
                A(indxm, indt1) = sqrtQk'*F1k ;
                A(indxm, indt2) = sqrtQk'*Gk ;
            end
            
            
            % % exterocpetive measurements % %
            
            % robot-to-robot
            imeas = 0;
            for ell = 1:nR
                for i = 1:size(zr(ell,:,:,k),3)
                    if zr(ell,3,i,k)<=0, continue, end
                    
                    Ri = Rr{ell,k}(2*i-1:2*i, 2*i-1:2*i);
                    sqrtRk = sqrtm(inv(Ri));
                    
                    imeas = imeas+1;
                    indr1 =  dofT + dofk*(k-2)+[3*ell-2:3*ell];
                    indr2 =  dofT + dofk*(k-2)+[3*zr(ell,3,i,k)-2 : 3*zr(ell,3,i,k)];
                    
                    if noUT
                        C = [cos(x(indr1(3))) -sin(x(indr1(3)));  sin(x(indr1(3))) cos(x(indr1(3))) ];
                        i_p_r = C'*(x(indr2(1:2))-x(indr1(1:2)));
                        rho = norm(i_p_r);
                        th = atan2(i_p_r(2),i_p_r(1));
                        zhat = [rho; th];
                        
                        if gRB==3
                            indxm =  dofT + dofk*(nsteps-1) + 2*nR*nR*(k-2) + [ 2*imeas-1:2*imeas ] ;
                            r =  [ zr(ell,1,i,k)-zhat(1,1);  pi_to_pi(zr(ell,2,i,k)-zhat(2,1))  ] ; %beware of sign
                            Htem = [ 1/norm(i_p_r)*i_p_r'; 1/norm(i_p_r)^2*i_p_r'*J' ];
                            H1 = - Htem* C'*[ eye(2)  J*(x(indr2(1:2)) - x(indr1(1:2))) ];
                            H2 = Htem* C'*[eye(2) zeros(2,1)];
                        elseif gRB==2
                            indxm =  dofT + dofk*(nsteps-1) + nR*nR*(k-2) + imeas ;
                            r =  pi_to_pi(zr(ell,2,i,k)-zhat(2,1)); %beware of sign
                            Htem = [ 1/norm(i_p_r)*i_p_r'; 1/norm(i_p_r)^2*i_p_r'*J' ];
                            H1 = - Htem* C'*[ eye(2)  J*(x(indr2(1:2)) - x(indr1(1:2))) ];
                            H2 = Htem* C'*[eye(2) zeros(2,1)];
                            H1 = H1(2,:);
                            H2 = H2(2,:);
                            sqrtRk = sqrtRk(2,2);
                        elseif gRB==1
                            indxm =  dofT + dofk*(nsteps-1) + nR*nR*(k-2) + imeas ;
                            r =  zr(ell,1,i,k)-zhat(1,1); %beware of sign
                            Htem = [ 1/norm(i_p_r)*i_p_r'; 1/norm(i_p_r)^2*i_p_r'*J' ];
                            H1 = - Htem* C'*[ eye(2)  J*(x(indr2(1:2)) - x(indr1(1:2))) ];
                            H2 = Htem* C'*[eye(2) zeros(2,1)];
                            H1 = H1(1,:);
                            H2 = H2(1,:);
                            sqrtRk = sqrtRk(1,1);
                        end
                        
                    else %unscented
                        if gRB==3
                            indxm =  dofT + dofk*(nsteps-1) + 2*nR*nR*(k-2) + [ 2*imeas-1:2*imeas ] ;
                        elseif gRB==2
                            indxm =  dofT + dofk*(nsteps-1) + nR*nR*(k-2) + imeas;
                            sqrtRk = sqrtRk(2,2);
                        elseif gRB==1
                            indxm =  dofT + dofk*(nsteps-1) + nR*nR*(k-2) + imeas;
                            sqrtRk = sqrtRk(1,1);
                        end
                        
                        [r,H1,H2] = meas_jacob_ut(@robot_measurement_model, x([indr1,indr2],1),  x([indr1,indr2],1),Plin([indr1,indr2],[indr1,indr2]),zr(ell,1:2,i,k),Ri,  doOC);
                    end
                    
                    A(indxm, [indr1, indr2])  = [sqrtRk'*H1 , sqrtRk'*H2] ;
                    b(indxm,1) = sqrtRk'*r;
                    
                end
            end
            
            
            % robot-to-target
            imeas = 0;
            for ell = 1:nR
                for i = 1:size(zt(ell,:,:,k),3)
                    if zt(ell,3,i,k)<=0, continue, end
                    
                    Ri = Rt{ell,k}(2*i-1:2*i, 2*i-1:2*i);
                    sqrtRk = sqrtm(inv(Ri));
                    
                    imeas = imeas+1;
                    indr1 =  dofT + dofk*(k-2)+ [3*ell-2:3*ell];
                    indr2 =  dofT + dofk*(k-2)+ 3*nR+ [dim_target*(zt(ell,3,i,k)-1)+1 : dim_target*zt(ell,3,i,k)];
                    
                    if noUT
                        C = [cos(x(indr1(3))) -sin(x(indr1(3)));  sin(x(indr1(3))) cos(x(indr1(3))) ];
                        i_p_t = C'*(x(indr2(1:2))-x(indr1(1:2)));
                        rho = norm(i_p_t);
                        th = atan2(i_p_t(2),i_p_t(1));
                        zhat = [rho; th];
                        
                        if gRB==3
                            indxm =  dofT + dofk*(nsteps-1) + 2*nR*nR*(nsteps-1) + 2*nR*nT*(k-2) + [2*imeas-1:2*imeas] ;
                            r =  [ zt(ell,1,i,k)-zhat(1,1);  pi_to_pi(zt(ell,2,i,k)-zhat(2,1))  ] ; %beware of sign
                            Htem = [ 1/norm(i_p_t)*i_p_t'; 1/norm(i_p_t)^2*i_p_t'*J' ];
                            H1 = - Htem* C'*[ eye(2)  J*(x(indr2(1:2)) - x(indr1(1:2))) ];
                            H2 = Htem* C'*[eye(2) zeros(2,dim_target-2)];
                        elseif gRB==2
                            indxm =  dofT + dofk*(nsteps-1) + nR*nR*(nsteps-1) + nR*nT*(k-2) + imeas ;
                            r =  pi_to_pi(zt(ell,2,i,k)-zhat(2,1)) ; %beware of sign
                            Htem = [ 1/norm(i_p_t)*i_p_t'; 1/norm(i_p_t)^2*i_p_t'*J' ];
                            H1 = - Htem* C'*[ eye(2)  J*(x(indr2(1:2)) - x(indr1(1:2))) ];
                            H2 = Htem* C'*[eye(2) zeros(2,dim_target-2)];
                            H1 = H1(2,:);
                            H2 = H2(2,:);
                            sqrtRk = sqrtRk(2,2);
                        elseif gRB==1
                            indxm =  dofT + dofk*(nsteps-1) + nR*nR*(nsteps-1) + nR*nT*(k-2) +imeas ;
                            r =  zt(ell,1,i,k)-zhat(1,1);  %beware of sign
                            Htem = [ 1/norm(i_p_t)*i_p_t'; 1/norm(i_p_t)^2*i_p_t'*J' ];
                            H1 = - Htem* C'*[ eye(2)  J*(x(indr2(1:2)) - x(indr1(1:2))) ];
                            H2 = Htem* C'*[eye(2) zeros(2,dim_target-2)];
                            H1 = H1(1,:);
                            H2 = H2(1,:);
                            sqrtRk = sqrtRk(1,1);
                        end
                        
                    else %unscented
                        if gRB==3
                            indxm =  dofT + dofk*(nsteps-1) + 2*nR*nR*(nsteps-1) + 2*nR*nT*(k-2) + [2*imeas-1:2*imeas] ;
                        elseif gRB==2
                            indxm =  dofT + dofk*(nsteps-1) + nR*nR*(nsteps-1) + nR*nT*(k-2) + imeas ;
                            sqrtRk = sqrtRk(2,2);
                        elseif gRB==1
                            indxm =  dofT + dofk*(nsteps-1) + nR*nR*(nsteps-1) + nR*nT*(k-2) + imeas ;
                            sqrtRk = sqrtRk(1,1);
                        end
                        
                        [r,H1,H2] = meas_jacob_ut(@robot_measurement_model,x([indr1,indr2],1), x([indr1,indr2],1),Plin([indr1,indr2],[indr1,indr2]),zt(ell,1:2,i,k),Ri, doOC);
                    end
                    
                    A(indxm, [indr1, indr2])  = [sqrtRk'*H1 , sqrtRk'*H2] ;
                    b(indxm,1) = sqrtRk'*r;
                    
                end
            end
            
            % robot-to-landmark
            
        end%k
        
        
        %%%%%%%%%%%%%%%%%%%%%
        % % solve linear system : Ax = b
        % using sparse QR with reordering
        [~,n] = size(A);
        qp = 1:n;
        [V,beta,p,R] = cs_qr(A) ;
        d = cs_qleft(V, beta, p, b) ;
        Rp = R(1:n,1:n);
        bp = d(1:n,1);
        delta_x =  Rp \ bp;
        delta_x(qp) = delta_x ;
        
        x = x + delta_x;
        xlin = x;
        
        if ~doOC
            ndx2 = norm(delta_x);
            if ndx2<epsilon, break, end
        else
            ndx3 = norm(delta_x);
            if ndx3<epsilon, break, end
        end
        
        
    end%end iteration
    
end

% % augment the first robot/target state
x = [x0; x];
xlin = [xlin0; xlin];

