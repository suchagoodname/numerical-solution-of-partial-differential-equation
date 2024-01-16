t0 = 0;
t1 = 1;
x0 = 0;
x1 = 1;
pde = model_data(t0, t1, x0, x1); %模型数据结构体
u_exact=@(x,t) cos(pi*t).*sin(pi*x);
n=[4,5,10,20,100];
for i=1:5
% 显格式
[X,T,U,o] = wave_equation_fd1d(n(i),(i+1)*n(i),pde);
showsolution(X,T,U,i); % 以二元函数方式显示数值解

% % 隐格式
% [X,T,U,r] = wave_equation_fd1d(n(i),2*n(i),pde,0.5);
% showsolution(X,T,U,i); % 以二元函数方式显示数值解
O(i)=o;
e(i)=getmaxerror(X,T,U,u_exact);
end

P=polyfit(O,e,1);
figure(6)
plot(O,e,'*',O,polyval(P,O),'b-')
xlabel('τ^2+h^2')
ylabel('误差')
title('误差与(τ^2+h^2)的关系')
p=P(1);%收敛阶
