function e = getmaxerror(X,T,U,u_exact)
[x,t] = meshgrid(X,T);
ue = u_exact(x',t');
e = max(max(abs(ue - U)));