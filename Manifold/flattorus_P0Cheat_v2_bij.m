function [P0,tvec2,tvec,hess] = flattorus_P0Cheat_v2_bij(theta,x,d,m)

%%% Inputs
    %%% N    - number of pts
    %%% dim  - intrinsic dimension
    %%% m    - oscillation modes
    %%% n    - n = 2*m*d, ambient space dimension

%%% Outputs
    %%% tvec2   - d*n*N orthonormal randomly-oriented tvec2
    %%% tvec    - d*n*N standard parametrized manifold tangent basis for dx/dphi
    %%% hess    - d*d*n*N hessian for x wrt phi at point cloud X 
    %%% P0      - N*n*n, projection matrix
    
N = size(x,1);
n = size(x,2);

%%% intrinsic data
% theta = [THETA, PHI];
%%% extrinsic data
% for ii = 1:dim
%     for jj = 1:m
%         x(:,(ii-1)*2*m+(jj-1)*2+1) = cos(jj*theta(:,ii));
%         x(:,(ii-1)*2*m+(jj-1)*2+2) = sin(jj*theta(:,ii));
%     end
% end

%%% tangent direction
tvec = zeros(d,n,N);
for kk = 1:d
    for ii = 1:d
        for jj = 1:m
            if kk == ii
                % tangent kk
                tvec(kk,(ii-1)*2*m+(jj-1)*2+1,:) = -jj*sin(jj*theta(:,ii));
                tvec(kk,(ii-1)*2*m+(jj-1)*2+2,:) = +jj*cos(jj*theta(:,ii));
            end
        end
    end
end
tvec = tvec/sqrt(sum((1:m).^2));

%%% normalized orthonormal unit tangent
tvec2 = tvec;

%%% analytic true projection matrix 
P0 = zeros(n,n,N);
for ii = 1:N
    P0(:,:,ii) = tvec(:,:,ii)'*pinv(tvec(:,:,ii)*tvec(:,:,ii)')*tvec(:,:,ii);
end
P0 = permute(P0,[3,2,1]); % P0 is N*n*n


%%% hessian matrix
hess = zeros(d,d,n,N);
for ii = 1:d
    for jj = 1:m
        % tensor i,i
        hess(ii,ii,2*m*(ii-1)+2*(jj-1)+1,:) = -jj^2*cos(jj*theta(:,ii));
        hess(ii,ii,2*m*(ii-1)+2*(jj-1)+2,:) = -jj^2*sin(jj*theta(:,ii));    
    end
end
hess = hess/sqrt(sum((1:m).^2));

end






