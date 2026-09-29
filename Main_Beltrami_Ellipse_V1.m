clear all,

addpath ./Manifold
addpath ./Operator_Matrix
addpath ./RBF_Kernel

%% given manifold parameters and information

%%% manifold parameter
d = 1;
n = 2;
N = 400;
% extrinsic data
am = 2;
%%% 1 == well-sampled data
%%% 2 == random data
DATA_INDEX = 2;

%% polynomial regression
%%% k-NN
k0 = 41;
%%% degree 
degr = 4;
%%% 1=LB, 2=Boch, 3=Hodge, 4=Lich
operator = 1;

%% generate required dataset
%%% generate data
rng(5);
[x,theta] = ellipse_generate(DATA_INDEX,am,N);
%%% projection matrix P0
[P0,tvec2] = ellipse_P0Cheat_v2(theta,x,am);
%%% analytic computation of solution and LB
[Lgf_true,func] = ellipse_analytic_computation(theta,am);

%% PDE solver

if operator == 1
    %%% constant for PDE
    cm = 1;
    %%% solve (cm*I-Delta)*u = f
    Rhs = cm*func - Lgf_true;
    
    %%% gmls, intrinsic method, 1/K weight 
    recurs = 1; % 1 is adaptive Knn, 0 is not.
    weightflag = 1; % 1 is 1/K weight.  
    % [Matrix_Lap1_Local] = Beltrami_v1p2_gmls_int_weightK(x,k0,degr,tvec2,operator,recurs,weightflag);
    [Matrix_Lap1_Local] = Beltrami_v1p2_gmls_int_weightK_2(x,k0,degr,tvec2,operator,recurs,weightflag);
    Matrix_Lap1_Local = Matrix_Lap1_Local - spdiags(sum(Matrix_Lap1_Local,2),0,N,N);
    
    %%% grbf, intrinsic method, 1/K weight, phs kernel
    recurs = 1; % 1 is adaptive Knn, 0 is not.
    weightflag = 1; % 1 is 1/K weight. 
    rbfshape.kerflag = 1; % 1 is phs.
    rbfshape.kappa = 3; % kappa for phs degree
    rbfshape.delta = 2e-3; % pinv or generalized Moore-Penrose inverse
    % [Matrix_Lap2_Local,iteration2] = Beltrami_v2p2_gmls_int_grbf(x,k0,degr,tvec2,operator,recurs,weightflag,rbfshape);
    [Matrix_Lap2_Local,iteration2] = Beltrami_v2p2_gmls_int_grbf_2(x,k0,degr,tvec2,operator,recurs,weightflag,rbfshape);
    Matrix_Lap2_Local = Matrix_Lap2_Local - spdiags(sum(Matrix_Lap2_Local,2),0,N,N);
    
    %%% (FE+LBf_col_num1), (IE+sol1), (Stable) is the inf norm
    INDEX_ERROR.FE = 1;
    INDEX_ERROR.IE = 1;
    INDEX_ERROR.Stable = 0;
    disp('gmls weight-1/K');
    [gmls1_1oK] = Error_IE_FE_Stable_compute(Matrix_Lap1_Local,cm,func,Rhs,Lgf_true,INDEX_ERROR);
    disp('grbf-fd weight-1/K+phs');
    [grbffd1_phs] = Error_IE_FE_Stable_compute(Matrix_Lap2_Local,cm,func,Rhs,Lgf_true,INDEX_ERROR);
end





