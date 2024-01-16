%question1
H=[0.25;0.2;0.1;0.05;0.01];
e=zeros(5,1);
for k=1:5
    h=H(k);
    n1=1/h-1;
    n2=2/h-1;
    x=0:h:1;
    y=0:h:2;
    B=diag(4*ones(n1,1))+diag(-ones(n1-1,1),-1)+diag(-ones(n1-1,1),1);
    A=cell(n2,n2);
    for i=1:n2
        for j=1:n2
            if i==j
                A(i,j)={B};
            elseif abs(i-j)==1
                A(i,j)={-eye(n1)};
            else
                A(i,j)={zeros(n1)};
            end
        end
    end
    A=cell2mat(A);
    F=-4*ones(n1*n2,1);
    F=h*h*F;
    for i=1:n1%底边
        F(i)=F(i)+x(i+1)^2;
    end
    for j=1:n2%右边
        F(n1*j)=F(n1*j)+(y(j+1)-1)^2;
    end
    for i=1:n1%顶边
        F(n1*(n2-1)+i)=F(n1*(n2-1)+i)+(x(i+1)-2)^2;
    end
    for j=1:n2%左边
        F(1+n1*(j-1))=F(1+n1*(j-1))+y(j+1)^2;
    end
    U=A\F;%数值解
    u=zeros(n1*n2,1);%精确解
    for j=1:n2
        for i=1:n1
            u(i+n1*(j-1))=(x(i+1)-y(j+1))^2;
        end
    end
    figure(k)
    x0=0:h:1;
    y0=0:h:2;
    subplot(1,3,1)
    [x0,y0]=meshgrid(x0,y0);
    z=(x0-y0).^2;
    surf(x0,y0,z)%精确解图像
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('精确解')

    subplot(1,3,2)
    Z=zeros(n2+2,n1+2);
    for i=1:n1+2
        Z(1,i)=x(i)^2;
        Z(n2+2,i)=(x(i)-2)^2;
    end
    for j=2:n2+1
        Z(j,1)=y(j)^2;
        Z(j,n1+2)=(y(j)-1)^2;
    end
    Z(2:n2+1,2:n1+1)=reshape(U,[n1,n2])';
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
figure(k+1)
plot(H,e,'*',H,polyval(P,H),'b-')
xlabel('步长')
ylabel('误差')
title('误差与步长的关系')
p=P(1);%收敛阶



