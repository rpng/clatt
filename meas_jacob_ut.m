function [r, H1,H2] = meas_jacob_ut(hmeas, x, xlin,P, zm, R, doOC)

J = [0 -1; 1 0];
dim = length(x);
npnts = dim*2 + 1;
scale = 3;      % scale = dim+lambda
kappa = scale-dim;

% generate samples
Pi = chol(P)' * sqrt(scale); 
xi = [x, repvec(x,dim)+Pi, repvec(x,dim)-Pi];

% transform the samples through measuremnt function
zi = feval(hmeas,xi);
z = repvec(zi(:,1),npnts);  % set first transformed sample as the base
ri = z-zi;
ri(2,:) = pi_to_pi(ri(2,:)); % compute correct residual
zi = z - ri; % offset zi from base according to correct residual

global gRB
if gRB==3
    zhat = (kappa*zi(:,1) + 0.5*sum(zi(:,2:end), 2)) / scale;
    r =  [ zm(1)-zhat(1,1);  pi_to_pi(zm(2)-zhat(2,1))  ] ; %beware of sign
elseif gRB==2
    zi = zi(2,:); %bearing only
    zhat = (kappa*zi(:,1) + 0.5*sum(zi(:,2:end), 2)) / scale;
    r =  pi_to_pi(zm(2)-zhat)  ; %beware of sign
elseif gRB==1
    zi = zi(1,:); %range only
    zhat = (kappa*zi(:,1) + 0.5*sum(zi(:,2:end), 2)) / scale;
    r =  zm(1)-zhat;   %beware of sign
end

dx = xi - repvec(x,npnts);
dz = zi - repvec(zhat,npnts);

Pxz = (2*kappa*dx(:,1)*dz(:,1)' + dx(:,2:end)*dz(:,2:end)') / (2*scale);
Pzz = (2*kappa*dz(:,1)*dz(:,1)' + dz(:,2:end)*dz(:,2:end)') / (2*scale);

Hproj = Pxz' / P; 

% %
if doOC    
    % %construct nullspace  % %
    if length(x)==6 %robot-to-robot
        N = [ eye(2), J*xlin(1:2,1);  0  0 1   ;    eye(2), J*xlin(4:5,1); 0 0 1 ];
    else %robot-to-target
        dim_target = length(x)-3;
        N = [ eye(2), J*xlin(1:2,1); 0 0 1 ;    eye(2), J*xlin(4:5,1);  zeros(dim_target-2,2),  kron(eye((dim_target-2)/2) , J)*xlin(6:end,1) ];
    end
    % % projection
    Hproj = Hproj * (eye(size(x,1))  -  N/(N'*N)*N') ;
end

H1 = Hproj(:,1:3);
H2 = Hproj(:,4:end);

