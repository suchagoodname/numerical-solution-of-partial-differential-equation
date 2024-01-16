function [X, T, U,o] = wave_equation_fd1d(NS, NT, pde, theta)
if nargin < 4
    theta = 0; % 默认用显格式
end
[X, h] = pde.space_grid(NS);
[T, tau] = pde.time_grid(NT);
N = length(X);
M = length(T);
r = pde.a()*tau/h;
% if r >=1 && theta==0
%    error('时间空间离散不满足显格式的稳定条件！') 
% end
r2 = r*r;
U = zeros(N,M);
% 初值条件
U(:,1) = pde.init_solution(X); 
U(2:end-1,2) = r2/2*(U(1:end-2,1)+U(3:end,1)) + (1-r2)*U(2:end-1,1)...
    + tau*pde.init_dt_solution(X(2:end-1));
% 边值条件
U(1,:) = pde.left_solution(T);
U(end,:) = pde.right_solution(T);

%% 隐格式
d = 1 + 2*ones(N-2,1)*r2*theta;
c = -ones(N-3,1)*r2*theta;
A2 = diag(c,-1) + diag(c,1)+diag(d);

d = 2 - 2*ones(N-2,1)*r2*(1-2*theta);
c = ones(N-3,1)*r2*(1-2*theta);
A1 = diag(c,-1) + diag(c,1)+diag(d);

d = -1 - 2*ones(N-2,1)*r2*theta;
c = ones(N-3,1)*r2*theta;
A0 = diag(c,-1) + diag(c,1)+diag(d);
for i=3:M
    RHS = tau*tau*pde.source(X,T(i));
    RHS(2) = RHS(2) + theta*r2*U(1,i) + ...
        (1-2*theta)*r2*U(1,i-1)+ theta*r2*U(1,i-2);
    RHS(end-1) = RHS(end-1) + theta*r2*U(end,i) + ...
        (1-2*theta)*r2*U(end,i-1)+ theta*r2*U(end,i-2);
    U(2:end-1,i) = A2\(A1*U(2:end-1,i-1) + A0*U(2:end-1,i-2)+RHS(2:end-1));
end
o=tau^2+h^2;
end
