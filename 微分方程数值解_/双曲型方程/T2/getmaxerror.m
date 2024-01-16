function e = getmaxerror(X,T,U,u)
[x,t] = meshgrid(X,T);
u = u(x',t');
e = max(max(abs(u - U)));
