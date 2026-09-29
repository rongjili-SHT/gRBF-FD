function [f,c,u,bc1,bc2] = torus_fcompute(THET,PHI,a,sol,cvar,bc1,bc2)

%%% Input
%%% theta       - intrinsic theta
%%% a           - a=1 is circle and a\neq 1 is ellipse, which analytic solution needs to be recomputed for ellipse
%%% bc1         - one-side boundary, 3x1 vector, contains
%%%               1. index of boundary; 2. beta1; 3. alpha1.
%%%               where BC is given as beta1*du/dv + alpha1*u = h1
%%% bc2         - the other side boundary, contains similar data
%%% cvar        - 1 is variable, 0 is constant

%%% Output
%%% f           - force
%%% h1          - nonhomogeneous boundary 1
%%% h2          - nonhomogeneous boundary 2


%%% 
N = size(THET,1);

%%% parameters for arbitrary solution
k = sol.k;
beta1 = sol.beta1;
alpha1 = sol.alpha1;

if cvar == 1
    %%% analytic results
    c = ones(N,1);
    c = 1.1*c + sin(THET).^2.*cos(PHI).^2; % nonconstant diffusion coefficient
    %%%%% analytic results
    % derivative of c 
    c_t = 2*sin(THET).*cos(THET).*cos(PHI).^2;   % this will give you div(c \grad)
    c_p = -2*sin(THET).^2.*sin(PHI).*cos(PHI);
else
    c = ones(N,1);
    c_t = zeros(N,1);
    c_p = zeros(N,1);
end



u = (alpha1*sin(k*PHI)-beta1./(a+cos(THET))*k.*cos(k*PHI)).*cos(THET);     % true solution u(theta)
u_t = -beta1*k*cos(k*PHI).*sin(THET).*cos(THET)./(a+cos(THET)).^2-(alpha1*sin(k*PHI)-beta1*k*cos(k*PHI)./(a+cos(THET))).*sin(THET);
u_tt = -2*beta1*k*cos(k*PHI).*sin(THET).^2.*cos(THET)./(a+cos(THET)).^3 ...
    -beta1*k*cos(k*PHI).*cos(THET).^2./(a+cos(THET)).^2 ...
    +2*beta1*k*cos(k*PHI).*sin(THET).^2./(a+cos(THET)).^2 ...
    -(alpha1*sin(k*PHI)-beta1*k*cos(k*PHI)./(a+cos(THET))).*cos(THET);
u_p = (alpha1*k*cos(k*PHI)+beta1*k^2*sin(k*PHI)./(a+cos(THET))).*cos(THET);
u_pp = (-alpha1*k^2*sin(k*PHI)+beta1*k^3*cos(k*PHI)./(a+cos(THET))).*cos(THET);

% g = [1 0; 0 (2+cos(theta))^2];   % Riemannian metric
gsin = sin(THET);
gcos = a + cos(THET);

% Analytic Lu = f = div(c \grad)
f = ( -gsin.*c.*u_t + gcos.*c_t.*u_t + gcos.*c.*u_tt + c_p./gcos.*u_p + c./gcos.*u_pp ) ./ gcos;


%%% nonhomogeneous boundary 
%%% Left is beta*du/dn + alpha*u = g
if exist('bc1','var') % exit or not
    indx_1 = bc1.indx_1;
    beta1_L = bc1.beta1_L;
    alpha1_L = bc1.alpha1_L;
    
    bc1.h1 = - beta1_L./(a+cos(THET(indx_1,:))).*u_p(indx_1,:) + alpha1_L*u(indx_1,:);
end

%%% Right is beta*du/dn + alpha*u = g
if exist('bc2','var') % exit or not
    indx_2 = bc2.indx_2;
    beta1_R = bc2.beta1_R;
    alpha1_R = bc2.alpha1_R;
    
    bc2.h2 = + beta1_R./(a+cos(THET(indx_2,:))).*u_p(indx_2,:) + alpha1_R*u(indx_2,:);
end

end



