clear all,

addpath ./Manifold
addpath ./Operator_Matrix
addpath ./RBF_Kernel

%% given manifold parameters and information

%%% manifold parameter
d = 3;
m = 2;
n = 2*m*d;
Ntheta = 50;
Nphi = 50*4;
N = Ntheta*Nphi;
%%% 1 == well-sampled data
%%% 2 == random data
DATA_INDEX = 2;

%% polynomial regression
%%% k-NN
k0 = 75;
%%% degree 
degr = 4;
%%% 1=LB, 2=Boch, 3=Hodge, 4=Lich
operator = 1;

%% generate required dataset
%%% generate data
rng(5);
[x,theta] = flattorus_generate(DATA_INDEX,d,m,N);
%%% projection matrix P0
[P0,tvec2] = flattorus_P0Cheat(theta,x,d,m);
%%% analytic computation of solution and LB
func = sin(theta(:,1)).*sin(theta(:,2)).*sin(theta(:,3));
%%% analytic forcing term f
Lgf_true = -3*sin(theta(:,1)).*sin(theta(:,2)).*sin(theta(:,3));

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
    rbfshape.delta = 1e-5; % pinv or generalized Moore-Penrose inverse
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





