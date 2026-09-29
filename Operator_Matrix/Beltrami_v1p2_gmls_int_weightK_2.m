function [Matrix_Lap,iteration] = Beltrami_v1p2_gmls_int_weightK_2(x,k0,degr,tvec,operator,recurs,weightflag)

%%% Inputs
    %%% x          - N*n data set with N data points on d-dim M in R^n
    %%% k          - k-nearest-neighbors
    %%% degr       - degree of polynomials used
    %%% tvec       - tvec is d*n*N, tangent vectors 
    %%% operator   - 1=Beltrami, 2=Bochner, 3=Hodge, 4=Lich.
    %%% recurs     - recurssive or not for adaptive KNN
    
%%% Outputs
    %%% Matrix_Lap - N*N Laplacian matrix
    %%% iteration  - k2NN_save, iter_numb
    
%%% Created by Shixiao Willing Jiang and Rongji Li Sep/04/2025    

%% k-nearest-neighbors

N = size(x,1);
n = size(x,2);
d = size(tvec,1); % tvec is d*n*N

if d > 1
    index =  generatemultiindex(degr,d); % index is d*term
else
    index = 0:1:degr;
end
term = size(index,2);


% all index used for function components of vector fields
denorm_all = sum(index,1); % remove norm effect
ind_a = zeros(d,1); % find index with power of (2,0,0,...) which for LB
for jj = 1:d
    ind_a(jj,1) = find((sum(index,1)==2)&(index(jj,:)==2));
end
    
if operator == 1
    a1 = zeros(N*k0,1);
    a2 = zeros(N*k0,1);
    j3 = zeros(N*k0,1);
end

%% adaptive sufficiently large enough Knn to make matrix close to SDD
k2NN_save = zeros(N,1); % local k2-NN neighbors 
iter_numb = zeros(N,1); % iteration number
k2NNtotal = 0; % k2 increment

bestofknn = ones(N,1)*k0; % optimal Knn before quite while
run_ornot = zeros(N,1); % if run over all Knn and delta 
ratrecord = zeros(N,1); % record the best ratio of w1 / max(w2:wk)
ratio_ref = 3.0; % threashold for the ratio of w1 / max(w2:wk)

%% start for each base point
pp = 1;
while pp <= N    
    
    while (1)
        
        %% iterative for many k and knn for w1<0
        if run_ornot(pp) < 0.5 % run_ornot(pp) == 0
            k = k0 + iter_numb(pp)*2;
        elseif run_ornot(pp) >= 0.5 % run_ornot(pp) == 1
            k = bestofknn(pp);
        end
        % d is N*k and contains the distances to the k nearest neighbors.
        % inds is N*k and contains the indices of the k nearest neighbors.
        % [~,inds] = knnCPU(x,x,k);
        [~,inds2] = knnCPU(x,x(pp,:),k);
        
        %% generation of PHI
        xx = x(inds2(1,:),:); % knn points
        x0 = x(pp,:); % center point
        pvec = tvec(:,:,pp); % pvec is d*n at center point, tvec is d*n*N
        xnorm = sqrt( sum(sum((xx-repmat(x0,k,1)).^2))/n/k );
        yy = (xx-repmat(x0,k,1))/xnorm; % k*n normalized data
        inyy = yy*pvec'; % k*d intrinsic normalized polynomial
        
        %%% Phi = [1 inyy inyy.^2 inyy.^3 ...]
        Phi = ones(k,term); % k*term
        %%% contant has 1, deg 1 has d, deg 2 has (d+1)*d/2
        %%% index =  generatemultiindex(degr,d); % index is d*term
        for ss = 1:term
            for rr = 1:d
                Phi(:,ss) = Phi(:,ss).*inyy(:,rr).^index(rr,ss);
            end
        end
        
        %% regression for function components of vector fields
        %%% OLS regression
        %%% compute weight function
        % [W] = gmls_weight_matrix(k,weightflag);
        % W = speye(k,k)*1/k;
        % W(1,1) = Am; % Am = 1 for first or Am = 10 for all next
        if iter_numb(pp) < 0.5
            [W] = gmls_weight_matrix(k,weightflag,1);
        else
            [W] = gmls_weight_matrix(k,weightflag,10);
        end
        PhiInv = (Phi'*W*Phi) \ (Phi'*W); % term*k
        PhiInv = PhiInv./xnorm.^repmat(denorm_all',1,k); % remove normalization effect
        
        %% construction of Laplacian
        if operator == 1
            weight0 = zeros(1,k);
            % 1,1 entry
            co11 = zeros(1,term);
            for jj = 1:d
                co11(1,ind_a(jj,1)) = 2;
            end
            weight0(1,1:k) = co11*PhiInv;
            if (weight0(1,1) <= 0) && ( abs(weight0(1,1)) / max(abs(weight0(1,2:k))) > 1.8*ratrecord(pp) ) % 1st record ratio, knn, and delta
                ratrecord(pp) = abs(weight0(1,1)) / max(abs(weight0(1,2:k)));
                bestofknn(pp) = k;
            end
        end
        % keyboard,
        
        %% to stop while
        if recurs == 1
            if run_ornot(pp) >= 0.5 % run_ornot(pp) == 1
                if weight0(1,1) >= 0
                    disp(['warning: for ' num2str(pp) ' th point, w1 >= 0, for all Knn and delta, study further why!']);
                end
                break;
            % two cases here: 1. value of weight0(1,1)<=-1 2. weight0 <= all others
            elseif (weight0(1,1) < -1) && ( abs(weight0(1,1)) > ratio_ref*max(abs(weight0(1,2:k))) ) % (weight0(1,1) <= 0) && (weight0(1,1) <= all others)
                break; % break while
            elseif iter_numb(pp) <= 10
                % disp([num2str(pp) ' th point is very poor, exitflag not 1, increase localknn!']);
                iter_numb(pp) = iter_numb(pp) + 1; % continue to iterate            
            else
                disp([num2str(pp) ' th point is boundary or very poor singular, study further why!']);
                run_ornot(pp) = 1;
                % disp([num2str(pp) ' th point is boundary or very poor singular, study further why!']);keyboard,
                % return;
            end
        elseif recurs == 0
            break; % break while
        end
    end %%% end while(1)
    
    %% Laplace-Beltrami operator optimization construction 
    if operator == 1
        %%% give values to Beltrami Matrix 
        a1(k2NNtotal+(1:k),1) = pp*ones(k,1);
        a2(k2NNtotal+(1:k),1) = inds2(1,:);
        j3(k2NNtotal+(1:k),1) = weight0(1,1:k); % matvec2((jj-1)+1,(qq-1)*k+(1:k));        
        %%% update k2NNtotal
        k2NNtotal = k2NNtotal + k; % update j3count
        % save local KNN
        k2NN_save(pp,1) = k;
    end
    % disp([num2str(pp) 'point;' 'C=' num2str(stdminerr(pp))]);
pp = pp + 1;    
end %%% end of all points

%%% function output result 1
if operator == 1
    %%% Matrix_Lich local LB sparse
    Matrix_Lap = sparse(a1,a2,j3,N,N);
end


%%% function output result 2
iteration.k2NN_save = k2NN_save;
iteration.iter_numb = iter_numb;
        
end









%% 备份代码，vector参考Gan Zecheng课题代码
