function showsolution(X, T, U,i)
figure(i)
[x, t] = meshgrid(X, T);
subplot(1,3,1)
mesh(x, t, U');
xlabel('X');
ylabel('T');
zlabel('数值解');
subplot(1,3,2)
U_r=cos(pi*t).*sin(pi*x);
mesh(x, t, U_r);
xlabel('X');
ylabel('T');
zlabel('精确解');
subplot(1,3,3)
mesh(x,t,abs(U'-U_r))%误差图
title('误差图');
end

