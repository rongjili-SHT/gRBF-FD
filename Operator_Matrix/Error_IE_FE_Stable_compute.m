function [method] = Error_IE_FE_Stable_compute(Matrix_Lap1_Local,cm,func,Rhs,Lgf_true,INDEX_ERROR)

%%% Inputs
    %%% Matrix_Lap1_Local  - Laplacian matrix N*N
    %%% cm                 - parameter of PDE
    %%% func               - true solution
    %%% Rhs                - right hand side for PDE
    %%% Lgf_true           - true Delta u
    %%% INDEX_ERROR        - FE,IE,Stable
    
%%% Outputs
    %%% method   - append to method structure
        %%% LBf_col_num1   - discrete Lu 
        %%% Lgf_true       - true Delta u 
        %%% FE             - max(abs(Lu - Delta*u))
        %%% sol1           - discrete (cm*I-L)\Rhs
        %%% func           - true sol.
        %%% err_sol_Linf   - max(abs(sol1-func)) 
        %%% Stable         - inf norm of Laplace matrix
        %%% timcount       - count the running time
    
%%% Created by Shixiao Willing Jiang Sep/04/2025    

% number of data
N = size(Matrix_Lap1_Local,1);

%%% numerical error of Laplace-Beltrami
if INDEX_ERROR.FE == 1    
    LBf_col_num1 = Matrix_Lap1_Local*func;    
    errP_local_Linf1 = max(abs(LBf_col_num1-Lgf_true));
    disp(['consistency = ' num2str(errP_local_Linf1)]);
    method.FE = errP_local_Linf1;
    method.LBf_col_num1 = LBf_col_num1;
end

if INDEX_ERROR.IE == 1
    sol1 = (cm*speye(N)-Matrix_Lap1_Local)\Rhs;    
    err_sol_Linf1 = max(abs(sol1-func));
    disp(['convergence = ' num2str(err_sol_Linf1)]);
    method.IE = err_sol_Linf1;
    method.sol1 = sol1;
end

if INDEX_ERROR.Stable == 1
    LBf_invs1 = (cm*speye(N)-Matrix_Lap1_Local)\speye(N);    
    inf_norm1 = norm(LBf_invs1,Inf);
    disp(['stability = ' num2str(inf_norm1)]);
    method.Stable = inf_norm1;
end
        
        
end