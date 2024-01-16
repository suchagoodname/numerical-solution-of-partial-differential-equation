pde = model_data(); %ģ�����ݽṹ��
n=[4,5,10,20,100];
H=zeros(1,5);
e=zeros(1,5);
u_exact=@(x,t) exp(-pi^2*t)*sin(x*pi)+x*(1-x);
for i=1:5

% ��ǰ��ָ�ʽ
[X,T,U1,r] = heat_equation_fd1d(n(i),n(i),pde,'forward');
showsolution(X,T,U1,i); % �Զ�Ԫ������ʽ��ʾ��ֵ��

% % ����ָ�ʽ
% [X,T,U2,r] = heat_equation_fd1d(n(i),n(i),pde,'backward');
% showsolution(X,T,U2,i); % �Զ�Ԫ������ʽ��ʾ��ֵ��

% % ����ԳƸ�ʽ���� Crank-Nicholson ��ʽ
% [X,T,U3,r] = heat_equation_fd1d(n(i),n(i),pde,'crank-nicholson');
% showsolution(X,T,U3,i); % �Զ�Ԫ������ʽ��ʾ��ֵ��
H(i)=log(r);
e(i)=getmaxerror(X,T,U1,u_exact);
end

P=polyfit(H,e,1);
figure(6)
plot(H,e,'*',H,polyval(P,H),'b-')
xlabel('����')
ylabel('���')
title('����벽���Ĺ�ϵ')
p=P(1);%������

