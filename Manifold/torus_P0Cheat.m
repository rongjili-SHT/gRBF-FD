function [P0,tvec] = torus_P0Cheat(theta,x,am)

N = size(x,1);
n = size(x,2);
d = size(theta,2);

%%% intrinsic data
% theta = [THETA, PHI];
%%% extrinsic data
% x = [(am+cos(THETA)).*cos(PHI), (am+cos(THETA)).*sin(PHI), sin(THETA)];
THETA = theta(:,1);
PHI = theta(:,2);

%%% tangent direction
tvec = zeros(d,n,N);
% tangent 1
tvec(1,1,:) = -sin(THETA).*cos(PHI);
tvec(1,2,:) = -sin(THETA).*sin(PHI);
tvec(1,3,:) = cos(THETA);
% tangent 2
tvec(2,1,:) = (am+cos(THETA)).*(-sin(PHI));
tvec(2,2,:) = (am+cos(THETA)).*cos(PHI);
tvec(2,3,:) = cos(THETA)*0;

%%% analytic true projection matrix 
P0 = zeros(n,n,N);
for ii = 1:N
    P0(:,:,ii) = tvec(:,:,ii)'*pinv(tvec(:,:,ii)*tvec(:,:,ii)')*tvec(:,:,ii);
end
P0 = permute(P0,[3,2,1]); % P0 is N*n*n

end


% theta1 = (0:1/Ntime_r:1-1/Ntime_r)'*2*pi;
% theta2 = (0:1/Ntime_r:1-1/Ntime_r)'*2*pi;
% 
% Theta1 = repmat(theta1,Ntime_r,1);
% Theta2 = kron(theta2,ones(Ntime_r,1));
% 
% x1 = (2+cos(Theta1)).*cos(Theta2);
% y1 = (2+cos(Theta1)).*sin(Theta2);
% z1 = sin(Theta1);
% 
% % tangent 1
% x2 = -sin(Theta1).*cos(Theta2);
% y2 = -sin(Theta1).*sin(Theta2);
% z2 = cos(Theta1);
% t2 = [x2 y2 z2];
% 
% % tangent 2
% x3 = (2+cos(Theta1)).*(-sin(Theta2));
% y3 = (2+cos(Theta1)).*cos(Theta2);
% z3 = cos(Theta1)*0;
% t3 = [x3 y3 z3];
% 
% % orthogonal 
% x4 = y2.*z3 - y3.*z2;
% y4 = x3.*z2 - x2.*z3;
% z4 = x2.*y3 - x3.*y2;
% t4 = [x4 y4 z4];




