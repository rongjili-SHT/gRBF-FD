function [P0,tvec] = flattorus_P0Cheat(theta,x,d,m)

%%% N    - number of pts
%%% dim  - intrinsic dimension
%%% m    - oscillation modes
%%% n    - n = 2*m*d, ambient space dimension

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

%%% analytic true projection matrix 
P0 = zeros(n,n,N);
for ii = 1:N
    P0(:,:,ii) = tvec(:,:,ii)'*pinv(tvec(:,:,ii)*tvec(:,:,ii)')*tvec(:,:,ii);
end
P0 = permute(P0,[3,2,1]); % P0 is N*n*n

end






