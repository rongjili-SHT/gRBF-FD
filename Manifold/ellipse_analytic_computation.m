function [Lgf_true,func] = ellipse_analytic_computation(theta,am)

%%% analytic computation of solution and LB acting on solution

%%% Inputs
    %%% theta   - dataset for intrinsic coordinate 
    %%% am      - radius of ellipse
        
%%% Outputs
    %%% P0_tilde   - Tangential projection matrix, N*n*n
    
    
%%% Riemaniann metric
g = sin(theta).^2+am^2*cos(theta).^2;   % Riemannian metric
gp = (2-2*am^2)*sin(theta).*cos(theta);  % g'
% gpp = (1-am^2)*2.*cos(2*theta);  % g''

%%% test for Laplace-Beltrami
%%% test function
a0 = 1;
func = sin(a0*theta);
funcp = a0*cos(a0*theta);   % u'(theta)
funcpp = -a0^2*sin(a0*theta);  % u"(theta)

% Laplace-Beltrami of func
Lgf_true = (funcpp -1/2*gp./g.*funcp)./g; % u_ij * g^ij 
% LB_f = -a0^2*sin(a0*theta); % for circle

end
