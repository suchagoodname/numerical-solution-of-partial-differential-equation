pde = model_data(); %模型数据结构体
n=[4,5,10,20,100];
H=zeros(1,5);
e=zeros(1,5);
u_exact=@(x,t) exp(-pi^2*t)*sin(x*pi)+x*(1-x);
for i=1:5

% 向前差分格式
[X,T,U1,r] = heat_equation_fd1d(n(i),n(i),pde,'forward');
showsolution(X,T,U1,i); % 以二元函数方式显示数值解

% % 向后差分格式
% [X,T,U2,r] = heat_equation_fd1d(n(i),n(i),pde,'backward');
% showsolution(X,T,U2,i); % 以二元函数方式显示数值解

% % 六点对称格式，即 Crank-Nicholson 格式
% [X,T,U3,r] = heat_equation_fd1d(n(i),n(i),pde,'crank-nicholson');
% showsolution(X,T,U3,i); % 以二元函数方式显示数值解
H(i)=log(r);
e(i)=getmaxerror(X,T,U1,u_exact);
end

P=polyfit(H,e,1);
figure(6)
plot(H,e,'*',H,polyval(P,H),'b-')
xlabel('步长')
ylabel('误差')
title('误差与步长的关系')
p=P(1);%收敛阶

