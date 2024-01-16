t0 = 0;
t1 = 0.5;
x0 = 0;
x1 = 0.5;
pde = model_data(t0, t1, x0, x1); %模型数据结构体
u_exact=@(x,t) sin(t).*sin(4*pi*x);
n=[5,10,20,50,100];
for i=1:5
% 显格式
[X,T,U,o] = wave_equation_fd1d(n(i),(i+1)*n(i),pde);
showsolution(X,T,U,i); % 以二元函数方式显示数值解

% % 隐格式
% [X,T,U,o] = wave_equation_fd1d(n(i),2*n(i),pde,0.1);
% showsolution(X,T,U,i); % 以二元函数方式显示数值解
O(i)=o;
e(i)=getmaxerror(X,T,U,u_exact);
end

P=polyfit(O,e,1);
figure(6)
plot(O,e,'*',O,polyval(P,O),'b-')
xlabel('步长')
ylabel('误差')
title('误差与步长的关系')
p=P(1);%收敛阶
