function [x,theta] = torus_generate(DATA_INDEX,am,Ntheta,Nphi)

% number of pts
N = Ntheta*Nphi;

%%% intrinsic data
%%% 1 == well-sampled data
%%% 2 == random data
if DATA_INDEX == 1 % well-sampled

    theta2 = [0:2*pi/Ntheta:2*pi-2*pi/Ntheta]'; % always [0,2*pi]
    phi2 = [0:2*pi/Nphi:2*pi-2*pi/Nphi]'; % maybe [0 <2*pi]
    
    [THETA,PHI] = meshgrid(theta2,phi2);
    THETA = reshape(THETA,N,1);
    PHI = reshape(PHI,N,1);
    
elseif DATA_INDEX == 2 % random

    THETA = rand(N,1)*2*pi;
    PHI = rand(N,1)*2*pi;
    
end

%%% intrinsic data
theta = [THETA, PHI];
%%% extrinsic data
x = [(am+cos(THETA)).*cos(PHI), (am+cos(THETA)).*sin(PHI), sin(THETA)];

 

end


