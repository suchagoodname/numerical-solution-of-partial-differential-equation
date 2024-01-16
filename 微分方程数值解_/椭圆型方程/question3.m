%question1
H=[pi/8;pi/16;pi/32;pi/40;pi/48;pi/56];
e=zeros(6,1);
for k=1:6
    h=H(k);
    n1=pi/h-1;
    n2=(pi/2)/h-1;
    x=0:h:pi;
    y=0:h:pi/2;
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
    F=zeros(n1*n2,1);
   for j=1:n2
        for i=1:n1
            F(i+n1*(j-1))=h*h*(cos(x(i+1)+y(j+1))+cos(x(i+1)-y(j+1)));
        end
    end
    for i=1:n1%底边
        F(i)=F(i)+cos(x(i+1));
    end
    for j=1:n2%右边
        F(n1*j)=F(n1*j)-cos(y(j+1));
    end
    for i=1:n1%顶边
        F(n1*(n2-1)+i)=F(n1*(n2-1)+i);
    end
    for j=1:n2%左边
        F(1+n1*(j-1))=F(1+n1*(j-1))+cos(y(j+1));
    end
    U=A\F;%数值解

    u=zeros(n1*n2,1);%精确解
    for j=1:n2
        for i=1:n1
            u(i+n1*(j-1))=cos(x(i+1))*cos(y(j+1));
        end
    end
    figure(k)
    x0=0:h:pi;
    y0=0:h:pi/2;
    subplot(1,3,1)
    [x0,y0]=meshgrid(x0,y0);
    z=cos(x0).*cos(y0);
    surf(x0,y0,z)%精确解图像
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('精确解')

    subplot(1,3,2)
    Z=zeros(n2+2,n1+2);
    for i=1:n1+2
        Z(1,i)=cos(x(i));
        Z(n2+2,i)=0;
    end
    for j=2:n2+1
        Z(j,1)=cos(y(j));
        Z(j,n1+2)=-cos(y(j));
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
plot(H,e,'*',H,polyval(P,H),'b-')%步长与误差的关系
xlabel('步长')
ylabel('误差')
title('误差与步长的关系')
p=P(1);%收敛阶





