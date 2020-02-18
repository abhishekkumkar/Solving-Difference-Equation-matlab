clc
clear all
syms x y
a=input('Enter the value of a: ');
h=input('Enter the value of h: ');
f=input('Enter the f(x,y): ');
x0=input('Enter the boundary value on X axis: ');
y0=input('Enter the boundary value on Y axis: ');
x1=input('Enter the boundary value on X=a: ');
y1=input('Enter the boundary value on Y=a: ');
p=input('Enter the coeficient of Uxx: ');
q=input('Enter the coeficient of Uyy: ');
r=input('Enter the coeficient of Ux: ');
s=input('Enter the coeficient of Uy: ');
n= a/h - 1;
X=(-2*p/h^2  -2*q/h^2 + r/h + s/h)*eye(n^2);
A=linspace(0,a,n+2);
B=linspace(0,a,n+2);
for i=1:n^2-n
    X(i,i+n)=p/h^2;
    X(i+n,i)=p/h^2-r/h;
end
for i=1:n^2-1
    if mod(i,n)==0
        X(i,i+1)=0;
        X(i+1,i)=0;
    else
        X(i,i+1)=q/h^2;
        X(i+1,i)=q/h^2-s/h;
    end
end
disp(X)
Y=ones(n^2,1);
for i=2:n+1 %y
for j=2:n+1 %x
    if i==2 
        if j==2
            Y((j-1)+(i-2)*n,1)=subs(subs(f,x,A(2)),y,B(2)) -x0*(q/h^2-s/h)-y0*(p/h^2-r/h);
        elseif j==n+1
             Y((j-1)+(i-2)*n,1)=subs(subs(f,x,A(n+1)),y,B(2))-x0*(q/h^2-s/h)-x1*(p/h^2);
        else 
             Y((j-1)+(i-2)*n,1)=subs(subs(f,x,A(j)),y,B(2))-x0*(q/h^2-s/h);
        end
    elseif i~=2 && i~=n+1
        if j==2
             Y((j-1)+(i-2)*n,1)=subs(subs(f,x,A(2)),y,B(i))-y0*(p/h^2-r/h);
        elseif j==n+1
             Y((j-1)+(i-2)*n,1)=subs(subs(f,x,A(n+1)),y,B(i))-x1*(p/h^2);
        else
             Y((j-1)+(i-2)*n,1)=subs(subs(f,x,A(j)),y,B(i));
        end
    else
        if j==2
             Y((j-1)+(i-2)*n,1)=subs(subs(f,x,A(2)),y,B(n+1))-y0*(p/h^2-r/h)-y1*(q/h^2);
        elseif j==n+1
             Y((j-1)+(i-2)*n,1)=subs(subs(f,x,A(n+1)),y,B(n+1))-x1*(p/h^2)-y1*(q/h^2);
        else 
             Y((j-1)+(i-2)*n,1)=subs(subs(f,x,A(j)),y,B(n+1))-y1*(q/h^2);
        end
    end
end
end
%disp(Y)
Z=inv(X)*Y;
%disp(Z)
U=zeros(n^2,1);
for i=2:n+1
for j=2:n+1
    U((i-2)*n+j-1)=exp(A(i)+B(j));
end
end
disp(U)
for i=1:n^2
    error(i,1)=abs(Z(i,1)-U(i,1));

end
%disp(max(error))
