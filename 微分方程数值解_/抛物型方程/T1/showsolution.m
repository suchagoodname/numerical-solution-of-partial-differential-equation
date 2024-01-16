function showsolution(X,T,U,i) 
figure(i)
[x,t] = meshgrid(X,T);
subplot(1,3,1)
mesh(x,t,U');
xlabel('X');
ylabel('T');
zlabel('U(X,T)');
title('��ֵ��')
subplot(1,3,2)
u_exact=2*exp(-pi^2*t./4).*sin(x.*2*pi);
mesh(x,t,u_exact)%��ȷ��ͼ��
xlabel('X');
ylabel('T');
zlabel('u_exact');
title('��ȷ��')
subplot(1,3,3)
mesh(x,t,abs(U'-u_exact))%���ͼ
title('���ͼ');
end


