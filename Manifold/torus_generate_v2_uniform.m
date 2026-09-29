function [x,theta] = torus_generate_v2_uniform(DATA_INDEX,am,Ntheta,Nphi)

% number of pts
N = Ntheta*Nphi;

%%% intrinsic data
%%% 1 == well-sampled data
%%% 2 == random data
if DATA_INDEX == 1 % well-sampled

    %%% uniform sampling 1
%     theta2 = [0:1/Ntheta:1-1/Ntheta]'+1/Ntheta/2; % maybe [0,pi]
%     options = optimoptions('fsolve','Display','none');
%     for ii = 1:Ntheta
%         % for sphere: theta2 = acos(1-2*theta2);
%         x0 = fsolve(@(x) am*x+sin(x)-2*pi*am*theta2(ii),[pi],options);
%         theta2(ii) = x0;
%     end    
%     % theta2 = [0:2*pi/Ntheta:2*pi-2*pi/Ntheta]'; % always [0,2*pi]
%     phi2 = [0:2*pi/Nphi:2*pi-2*pi/Nphi]'; % maybe [0 <2*pi]
    
    %%% sampling 2 such that theta close to pi more data
    theta2 = [0:1/Ntheta:1-1/Ntheta]'+1/Ntheta/2; % maybe [0,pi]
    options = optimoptions('fsolve','Display','none');
    for ii = 1:Ntheta
        % for sphere: theta2 = acos(1-2*theta2);
        x0 = fsolve(@(x) am*x-sin(x)-2*pi*am*theta2(ii),[pi],options);
        theta2(ii) = x0;
    end    
    % theta2 = [0:2*pi/Ntheta:2*pi-2*pi/Ntheta]'; % always [0,2*pi]
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



%% ±¸·Ý´úÂë
% c = -1;
% x = fsolve(@(x) x+c,[-5]);

% function F = myfun(x,c)
% F = [2*x(1) - x(2) - exp(c*x(1))
%     -x(1) + 2*x(2) - exp(c*x(2))];
% end