function [PHI, sqrtQ,yhat] = prop_jacob_ut(fmotion,x, xlin, P,v,omega,Q,dt, doOC)

% sample entire state
x = [x; v;omega]; %
P = blkdiag(P,Q); %

dim = length(x);
npnts = dim*2 + 1;
scale = 3;      %scale = dim+lambda
kappa = scale-dim;

% generate samples
Pi = chol(P)' * sqrt(scale); 
xi = [x, repvec(x,dim)+Pi, repvec(x,dim)-Pi];

% transform the samples through measuremnt function
yi = feval(fmotion,xi,dt);
y = repvec(yi(:,1),npnts);  % set first transformed sample as the base
ri = y - yi;
ri(3,:) = pi_to_pi(ri(3,:));  % compute correct residual
yi = y - ri; % offset yi from base according to correct residual

yhat = (kappa*yi(:,1) + 0.5*sum(yi(:,2:end),2)) / scale;

dx = xi-repvec(x,npnts);
dy = yi-repvec(yhat,npnts);

Pxy = (2*kappa*dx(:,1)*dy(:,1)' + dx(:,2:end)*dy(:,2:end)') / (2*scale);
Pyy = (2*kappa*dy(:,1)*dy(:,1)' + dy(:,2:end)*dy(:,2:end)') / (2*scale);

F = (P\Pxy)';

Pee = Pyy - F*P*F'; %Qstar
if min(eig(Pee))<1e-6
    Pee = Pee + eye(size(Pee))*1e-6;
end
sqrtQ =  sqrtm(eye(size(Pee))/Pee);


% correspondence to the normal ekf propagation
PHI = F(:,1:end-2);
G = F(:,end-1:end);

% % 
if doOC
    J = [0 -1; 1 0];
    PHIproj = [ PHI,  - eye(3)] ;
    % %construct nullspace
    N = [ eye(2), J*xlin(1:2,1);  0  0 1   ;    eye(2), J*xlin(4:5,1); 0 0 1 ];
    % % projection
    PHIproj = PHIproj * (eye(size(xlin,1))  -  N/(N'*N)*N') ;
    PHI = PHIproj(:,1:3);
end


