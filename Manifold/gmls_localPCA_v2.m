function [P0_tilde,tvec_tilde2] = gmls_localPCA_v2(x,k,d,n,degr,indxB)

%%% Inputs
    %%% x       - N-by-[n_R] data set with N data points in R^n_R
    %%% k       - number of nearest neighbors to use
    %%% d       - intrinsic dimension, must be integer
    %%% n       - ambient dimension
    %%% degr    - degree of polynomial used
    %%% indxB   - the indices of points to be estimated
        
%%% Outputs
    %%% P0_tilde   - Tangential projection matrix, N*n*n
   
% %%% Weight function with parameter
% % weight = @(d,D) (1-d/D)^4*((4*d)/D + 1);
% % weight = @(d,D) exp(-D^2*d^2);
% weight = @(d,D) 1;
% % weight = @(d,D) max((1-d/D)^4,0);
% Interval_para = 1;

% N = size(x,1); 
NB = length(indxB);

if d > n % d <= n
    disp('d must be <= n');
    return;
end
if d == n
    P0_tilde = zeros(n,n,NB);  
    for i = 1:NB
        P0_tilde(:,:,i) = eye(n);
    end
    return;
end
minterm = nchoosek(degr+d,d);
if k <= minterm % k > D = (d+1)*d/2
    disp(['k must be > ' num2str(minterm)]);
    return;
end

%%% k nearest neighbors
[~,inds] = knnCPU(x,x(indxB,:),k); % NB*k

%%% index for gmls
if d > 1
    index =  generatemultiindex(degr,d); % index is d*term
else
    index = 0:1:degr;
end
term = size(index,2);

% all index used for function components of vector fields
denorm = sum(index,1); % remove norm effect
ind_alin = zeros(d,1);
for jj = 1:d
    ind_alin(jj) = find((sum(index,1)==1)&(index(jj,:)==1));
end

%% local SVD 

P0_tilde = zeros(n,n,NB);    
tvec_tilde = zeros(d,n,NB); % tvec is d*n*N
tvec_tilde2 = zeros(d,n,NB); % tvec is d*n*N

for i = 1:NB
    
    temp = x(inds(i,:),:)'; % n*k
    tempmean = mean(temp,2); % n*1
    temp_rough = temp - repmat(tempmean,1,k); % n*k
    
    [U,S,~] = svd(temp_rough);
    [~,IX] = sort(diag(S),'descend');
    U = U(:,IX); % rotation matrix
    
    %%% rough tangent space for representation of geodesic normal
    %%% coordiante
    t_rough = U(:,1:d); % n*d
    n_rough = U(:,d+1:n); % n*(n-d)
    tvec_tilde(:,:,i) = t_rough'; % tvec is d*n*N
    
    %%% correction
    xx = temp'; % k*n, knn points
    x0 = x(indxB(i),:); % center point
    xnorm = sqrt( sum(sum((xx-repmat(x0,k,1)).^2))/n/k ); % xx is k*n
    yy = (xx-repmat(x0,k,1))/xnorm; % k*n normalized data
    inyy = yy*U; % k*n with leading d in tvec, latter n-d in nvec, intrinsic normalized polynomial
    
    %%% Phi = [1 yy yy.^2 yy.^3 ...]
    Phi = ones(k,term); % k*term
    
    %%% contant has 1, deg 1 has d, deg 2 has (d+1)*d/2
    %%% index =  generatemultiindex(degr,d); % index is d*term
    for ss = 1:term
        for rr = 1:d
            Phi(:,ss) = Phi(:,ss).*inyy(:,rr).^index(rr,ss);
        end
    end
        
    %%% regress for manifolds using PHI inverse, for latter n-d components
    for jj = 1:n-d
        
        bb = (Phi'*Phi) \ (Phi'*inyy(:,d+jj)); % regress for all deg
        bb = bb./xnorm.^(denorm'-1); % term*1, remove normalization effect
        
        for mm = 1:d
            tvec_tilde(mm,:,i) = tvec_tilde(mm,:,i) + bb(ind_alin(mm))*n_rough(:,jj)'; % tvec is d*n*N
        end
    end
    
    %%% generate projection matrix by tangent vectors
    P0_tilde(:,:,i) = tvec_tilde(:,:,i)'*pinv(tvec_tilde(:,:,i)*tvec_tilde(:,:,i)')*tvec_tilde(:,:,i);
    %%% generate orthonormal tangent vectors by projection matrix
    [U2,S2,~] = svd(P0_tilde(:,:,i));
    [~,IX2] = sort(diag(S2),'descend');
    U2 = U2(:,IX2); % eigenvector
    tvec_tilde2(:,:,i) = U2(:,1:d)'; % tvec is d*n*N
% keyboard,
   
end

P0_tilde = permute(P0_tilde,[3,2,1]); % P0 is N*n*n

    
end




%% ±¸·Ý´úÂë


