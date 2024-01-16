t0 = 0;
t1 = 1;
x0 = 0;
x1 = 1;
pde = model_data(t0, t1, x0, x1); %ģ�����ݽṹ��
u_exact=@(x,t) cos(pi*t).*sin(pi*x);
n=[4,5,10,20,100];
for i=1:5
% �Ը�ʽ
[X,T,U,o] = wave_equation_fd1d(n(i),(i+1)*n(i),pde);
showsolution(X,T,U,i); % �Զ�Ԫ������ʽ��ʾ��ֵ��

% % ����ʽ
% [X,T,U,r] = wave_equation_fd1d(n(i),2*n(i),pde,0.5);
% showsolution(X,T,U,i); % �Զ�Ԫ������ʽ��ʾ��ֵ��
O(i)=o;
e(i)=getmaxerror(X,T,U,u_exact);
end

P=polyfit(O,e,1);
figure(6)
plot(O,e,'*',O,polyval(P,O),'b-')
xlabel('��^2+h^2')
ylabel('���')
title('�����(��^2+h^2)�Ĺ�ϵ')
p=P(1);%������
