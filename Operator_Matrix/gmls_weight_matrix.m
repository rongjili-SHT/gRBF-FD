function [W] = gmls_weight_matrix(k,weightflag,Am)

%%% Inputs
    %%% k          - k nearest neighbors
    %%% weightflag - 1 is 1/K,  
    
%%% Outputs
    %%% W      - K*K symmetric positive definite matrix 
    
    
if weightflag == 1
    W = speye(k,k)*1/k;
    W(1,1) = Am; % Am = 1 for first standard Knn iteration, Am = 10 for next
end

end