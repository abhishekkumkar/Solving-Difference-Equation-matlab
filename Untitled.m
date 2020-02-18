clc 
clear all
syms x
a=input('Enter the value of a: ');
b=input('Enter the value of b: ');
fa=input('Enter the value of f(a): ');
fb=input('Enter the value of f(b): ');
p=input('Enter coeficient of y2: ');
q=input('Enter coeficient of y1: ');
r=input('Enter coeficient of y: ');
f=input('Enter the function f(x): ');               
h=input('Enter the value of h: ');
n=(b-a)/h;
X1=linspace(a,b,n+1);
X=(-2*p/h^2 + q/h + r)*eye(n-1);
for i=1:n-2
    X(i,i+1)= p/h^2 - q/h ;
    X(i+1,i)= p/h^2;
end
Y=ones(n-1,1);
Y(1,1)=subs(f,x,X1(2))-fa*(p/h^2-q/h);
Y(n-1,1)=subs(f,x,X1(n))-fb*p/h^2;
for i=3:length(X1)-2
    Y(i-1,1)=subs(f,x,X1(i));
end
Z=inv(X)*Y;
disp(X)
disp(Y)
disp(Z)
 


