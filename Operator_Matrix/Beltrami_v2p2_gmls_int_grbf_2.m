function [Matrix_Lap,iteration] = Beltrami_v2p2_gmls_int_grbf_2(x,k0,degr,tvec,operator,recurs,weightflag,rbfshape)

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

%% index of multipolynomial and degree 2 for beltrami

N = size(x,1);
% n = size(x,2);
d = size(tvec,1); % tvec is d*n*N

if d > 1
    index =  generatemultiindex(degr,d); % index is d*term
else
    index = 0:1:degr;
end
term = size(index,2);


% all index used for function components of vector fields
% denorm_all = sum(index,1); % remove norm effect
ind_a = zeros(d,1); % find index with power of (2,0,0,...) which for LB
for jj = 1:d
    ind_a(jj,1) = find((sum(index,1)==2)&(index(jj,:)==2));
end
    
if operator == 1
    a1 = zeros(N*k0,1);
    a2 = zeros(N*k0,1);
    j3 = zeros(N*k0,1);
end

%% rbf kernel, power, and shape parameter
if rbfshape.kerflag == 1 %%% polyharmonic spline
    grbf = phs();
end
% keyboard,
%% adaptive sufficiently large enough Knn to make matrix close to SDD
k2NN_save = zeros(N,1); % local k2-NN neighbors 
iter_numb = zeros(N,1); % iteration number for Knn
k2NNtotal = 0; % k2 increment
deltaiter = zeros(N,1); % adaptive delta for pseudo-inverse of rbf

bestofknn = ones(N,1)*k0; % optimal Knn before quite while
bestdelta = ones(N,1)*rbfshape.delta; % optimal delta before quite while
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
        
        %% generation of Phi_p, multipolynomial
        xx = x(inds2(1,:),:); % knn points
        x0 = x(pp,:); % center point
        pvec = tvec(:,:,pp); % pvec is d*n at center point, tvec is d*n*N
        
        yy = (xx-repmat(x0,k,1)); % k*n original unnormalized data
        inyy = yy*pvec'; % k*d intrinsic normalized polynomial
        rdist = pdist2(inyy,inyy);
        xnorm = max(max(rdist));
        inyy = yy*pvec'/ xnorm; % k*d intrinsic normalized polynomial
        
        %%% Phi = [1 inyy inyy.^2 inyy.^3 ...]
        Phi_p = ones(k,term); % k*term
        %%% contant has 1, deg 1 has d, deg 2 has (d+1)*d/2
        %%% index =  generatemultiindex(degr,d); % index is d*term
        for ss = 1:term
            for rr = 1:d
                Phi_p(:,ss) = Phi_p(:,ss).*inyy(:,rr).^index(rr,ss);
            end
        end
        
        %% gmls regression for function 
        %%% OLS regression
        %%% compute weight function
        % W = speye(k,k)*1/k;
        % W(1,1) = Am; % Am = 1 for first or Am = 10 for all next
        if iter_numb(pp) < 0.5
            [W] = gmls_weight_matrix(k,weightflag,1);
        else
            [W] = gmls_weight_matrix(k,weightflag,10);
        end
        PhiInv = (Phi_p'*W*Phi_p) \ (Phi_p'*W); % term*k
        GrbfInv = speye(k,k) - Phi_p*PhiInv; % residual of gmls
        % PhiInv = PhiInv./xnorm.^repmat(denorm_all',1,k); % remove normalization effect
        
        %% radial basis functoin
        rdist = pdist2(inyy,inyy);
        Phi_r = grbf.phi(rdist,rbfshape.kappa);
        Phi_r_Beltrm = grbf.beltrami(rdist(1,:),rbfshape.kappa,d);        
        % keyboard,
        
        %% construction of Laplacian
        if operator == 1
            weight0 = zeros(1,k);
            % 1,1 entry
            co11 = zeros(1,term);
            for jj = 1:d
                co11(1,ind_a(jj,1)) = 2;
            end
            % I_f=Pnorm*bnorm=Pnorm*PhiInv*f
            % Delta I_f(x0) = (2/xnorm^2*PhiInv)*f
            weight0(1,1:k) = co11*PhiInv/xnorm.^2 ; % gmls part
            weighttemp = weight0(1,1:k);
            
            if run_ornot(pp) < 0.5 % run_ornot(pp) == 0
                %%% pinv or generalized Moore-Penrose inverse
                % PhiInv_r = pinv(Phi_r,rbfshape.delta);
                PhiInv_r = (Phi_r'*W*Phi_r + rbfshape.delta.^2.*speye(k,k)) \ (Phi_r'*W);
                % If = Phi_r*a+Phi_p*b, b=PhiInv*f
                % a = PhiInv_r*(I-Pnorm*PhiInv)*f
                weight0(1,1:k) = weighttemp + Phi_r_Beltrm*PhiInv_r*GrbfInv/xnorm.^2 ;% rbf part
                if (weight0(1,1) <= 0) && ( abs(weight0(1,1)) / max(abs(weight0(1,2:k))) > 1.8*ratrecord(pp) ) % 1st record ratio, knn, and delta
                    ratrecord(pp) = abs(weight0(1,1)) / max(abs(weight0(1,2:k)));
                    bestofknn(pp) = k;
                    bestdelta(pp) = rbfshape.delta;
                end
            elseif run_ornot(pp) >= 0.5 % run_ornot(pp) == 1
                PhiInv_r = (Phi_r'*W*Phi_r + bestdelta(pp).^2.*speye(k,k)) \ (Phi_r'*W);
                weight0(1,1:k) = weighttemp + Phi_r_Beltrm*PhiInv_r*GrbfInv/xnorm.^2 ;% rbf part
            end
            
            %%% auto tune of delta for pseudo inverse of rbf
            if (recurs == 1) && (run_ornot(pp) < 0.5)
                deltaiter(pp) = 0;
                % if (weighttemp(1,1) < -1) && ( abs(weighttemp(1,1)) > ratio_ref*max(abs(weighttemp(1,2:k))) ) % gmls is stable such that Knn is proper
                while (1)
                    if (weight0(1,1) < -1) && ( abs(weight0(1,1)) > ratio_ref*max(abs(weight0(1,2:k))) ) % grbf-fd is stable
                        break; % break while
                    % elseif 2.^deltaiter(pp)*rbfshape.delta > 7e-3 % grbf-fd is not stable but delta is large enough
                    elseif deltaiter(pp) > 3 % grbf-fd is not stable but delta is large enough    
                        break;
                    else
                        deltaiter(pp) = deltaiter(pp) + 1;
                        % PhiInv_r = pinv(Phi_r,rbfshape.delta);
                        PhiInv_r = (Phi_r'*W*Phi_r + (2.^deltaiter(pp)*rbfshape.delta).^2.*speye(k,k)) \ (Phi_r'*W);
                        weight0(1,1:k) = weighttemp + Phi_r_Beltrm*PhiInv_r*GrbfInv/xnorm.^2 ;% rbf part
                        if (weight0(1,1) <= 0) && ( abs(weight0(1,1)) / max(abs(weight0(1,2:k))) > 1.8*ratrecord(pp) ) % 1st record ratio, knn, and delta
                            ratrecord(pp) = abs(weight0(1,1)) / max(abs(weight0(1,2:k)));
                            bestofknn(pp) = k;
                            bestdelta(pp) = 2.^deltaiter(pp)*rbfshape.delta;
                        end
                    end
                end
                % end % end for gmls
            end % end for recurs
            
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
iteration.deltaiter = deltaiter;
        
iteration.bestofknn = bestofknn; % optimal Knn before quite while
iteration.bestdelta = bestdelta; % optimal delta before quite while
iteration.run_ornot = run_ornot; % if run over all Knn and delta 
iteration.ratrecord = ratrecord; % record the best ratio of w1 / max(w2:wk)

end







