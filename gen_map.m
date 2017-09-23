function xL = gen_map(nL,v,omega,rmin,rmax, nSteps,dt)

if omega~=0

    % radius
    r = v/omega;

    for i=1:nL

        if rand  > .5 %outside of circle
            disp('landmarks: ourside of circle')
            rho = r + rmin + (rmax-rmin)*rand;
            th = 2*pi*rand;
        else %inside of circle
            disp('landmarks: inside of circle')
            rho = r - rmax + (rmax-rmin)*rand;
            th = 2*pi*rand;
        end

        [x,y] = pol2cart(th,rho);
        xL(:,i) = [x;y+r]; %shift y w/ r since robot starts at (0,0)

    end
    
end
%
xL(:,1) = [0;r];%one landmark at center


% % %exploration case
if omega==0 
       
    for i=1:nL
        if round(rand)
            rho = rmin + (rmax-rmin)*rand;
            th = 2*pi*rand;
        else
            rho = - rmax + (rmax-rmin)*rand;
            th = 2*pi*rand;
        end
        [x,y] = pol2cart(th,rho);
        yy(1,i) = rho;%y;
    end
    
    dtraj = v*dt*nSteps;
    
    xL = [dtraj*rand(1,nL); yy ];
    

end
