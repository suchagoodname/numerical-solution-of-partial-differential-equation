%question1
H=[0.25;0.2;0.1;0.05;0.01];
e=zeros(5,1);
for k=1:5
    h=H(k);
    n=1/h-1;
    x=1:h:2;
    y=0:h:1;
    B=diag(4*ones(n,1))+diag(-ones(n-1,1),-1)+diag(-ones(n-1,1),1);
    A=cell(n,n);
    for i=1:n
        for j=1:n
            if i==j
                A(i,j)={B};
            elseif abs(i-j)==1
                A(i,j)={-eye(n)};
            else
                A(i,j)={zeros(n)};
            end
        end
    end
    A=cell2mat(A);
    F=zeros(n^2,1);
    for i=1:n%底边
        F(i)=F(i)+2*log(x(i+1));
    end
    for j=1:n%右边
        F(n*j)=F(n*j)+log(4+y(j+1)^2);
    end
    for i=1:n%顶边
        F(n*(n-1)+i)=F(n*(n-1)+i)+log(x(i+1)^2+1);
    end
    for j=1:n%左边
        F(1+n*(j-1))=F(1+n*(j-1))+log(1+y(j+1)^2);
    end
    U=A\F;%数值解
    u=zeros(n^2,1);%精确解
    for j=1:n
        for i=1:n
            u(i+n*(j-1))=log(x(i+1)^2+y(j+1)^2);
        end
    end
    figure(k)
    x0=1:h:2;
    y0=0:h:1;
    subplot(1,3,1)
    [x0,y0]=meshgrid(x0,y0);
    z=log(x0.^2+y0.^2);
    surf(x0,y0,z)%精确解图像
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('精确解')

    subplot(1,3,2)
    Z=zeros(n+2,n+2);
    for i=1:n+2
        Z(1,i)=2*log(x(i));
        Z(n+2,i)=log(x(i)^2+1);
    end
    for j=2:n+1
        Z(j,1)=log(1+y(j)^2);
        Z(j,n+2)=log(4+y(j)^2);
    end
    Z(2:n+1,2:n+1)=reshape(U,[n,n])';
    [x,y]=meshgrid(x,y);
    surf(x,y,Z)
    xlabel('x');
    ylabel('y');
    zlabel('Z');
    title('数值解')
    e(k)=norm(U-u,'Inf');
    subplot(1,3,3)
    surf(x0,y0,abs(Z-z))%误差图
    title('误差图')
end
H=log(H);
e=log(e);
P=polyfit(H,e,1);
figure(6)
plot(H,e,'*',H,polyval(P,H),'b-')
xlabel('步长')
ylabel('误差')
title('误差与步长的关系')
p=P(1);%收敛阶




