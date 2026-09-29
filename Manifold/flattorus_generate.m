function [x,theta] = flattorus_generate(DATA_INDEX,dim,m,N)

%%% N    - number of pts
%%% dim  - intrinsic dimension
%%% m    - oscillation modes
%%% n    - n = 2*m*d, ambient space dimension

%%% intrinsic data
theta = zeros(N,dim);

%%% intrinsic data
%%% 1 == well-sampled data
%%% 2 == random data
if DATA_INDEX == 1 % well-sampled

    
elseif DATA_INDEX == 2 % random

    for ii = 1:dim
        theta(:,ii) = rand(N,1)*2*pi;
    end
    
end

%%% extrinsic data
% x = 1/sqrt(1^2 + ... + m^2)*[cost1 sint1 ... cosmt1 sinmt1 ... costd sintd ... cosmtd sinmtd];
n = 2*m*dim;
x = zeros(N,n);
for ii = 1:dim
    for jj = 1:m
        x(:,(ii-1)*2*m+(jj-1)*2+1) = cos(jj*theta(:,ii));
        x(:,(ii-1)*2*m+(jj-1)*2+2) = sin(jj*theta(:,ii));
    end
end
x = x/sqrt(sum((1:m).^2));


end


