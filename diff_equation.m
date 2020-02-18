clc
clear all
syms x y

a=input('Enter the value of a: ');
h=input('Enter the value of h: ');
f=input('Enter the f(x,y): ');
fx0=input('Enter the boundary function on X axis: ');
fy0=input('Enter the boundary function on Y axis: ');
fx1=input('Enter the boundary function on X=a: ');
fy1=input('Enter the boundary function on Y=a: ');
p=input('Enter the coeficient of Uxx: ');
q=input('Enter the coeficient of Uyy: ');
r=input('Enter the coeficient of Ux: ');
s=input('Enter the coeficient of Uy: ');

n= (a/h) - 1;
X=(-2*p/(h^2) - 2*q/(h^2))*eye(n^2);
A=linspace(0,a,n+2);
B=linspace(0,a,n+2);  %verified

for i=1:(n^2) - n
    X(i,i+n)=q/(h^2) + r/(2*h);
    X(i+n,i)=q/(h^2) - r/(2*h);   %verified
end

for i=1:(n^2) - 1
    if mod(i,n)==0
        X(i,i+1)=0;
        X(i+1,i)=0;
    else
        X(i,i+1)=p/(h^2) + r/(2*h);
        X(i+1,i)=p/(h^2) - r/(2*h);   %verified
    end
end

%disp(X)                                 

Y=ones(n^2,1);
for j=2:n+1     %y
for i=2:n+1     %x
    
    if j==2
        
        if i==2 
            Y(i-1 + (j-2)*n,1)=subs(subs(f,x,A(i)),y,B(j)) - (p/(h^2)-r/(2*h))*subs(subs(fy0,x,A(i-1)),y,B(j)) - (q/(h^2)-r/(2*h))*subs(subs(fx0,x,A(i)),y,B(j-1));
        elseif i==n+1
            Y(i-1 + (j-2)*n,1)=subs(subs(f,x,A(i)),y,B(j)) - (p/(h^2)+r/(2*h))*subs(subs(fx1,x,A(i+1)),y,B(j)) - (q/(h^2)-r/(2*h))*subs(subs(fx0,x,A(i)),y,B(j-1));
        else
            Y(i-1 + (j-2)*n,1)=subs(subs(f,x,A(i)),y,B(j)) -  (q/(h^2)-r/(2*h))*subs(subs(fx0,x,A(i)),y,B(j-1));
        end
        
    elseif j==n+1
        
        if i==2 
            Y(i-1 + (j-2)*n,1)=subs(subs(f,x,A(i)),y,B(j)) - (p/(h^2)-r/(2*h))*subs(subs(fy0,x,A(i-1)),y,B(j)) - (q/(h^2)+r/(2*h))*subs(subs(fy1,x,A(i)),y,B(j+1));
        elseif i==n+1
            Y(i-1 + (j-2)*n,1)=subs(subs(f,x,A(i)),y,B(j)) - (p/(h^2)+r/(2*h))*subs(subs(fx1,x,A(i+1)),y,B(j)) - (q/(h^2)+r/(2*h))*subs(subs(fy1,x,A(i)),y,B(j+1));
        else
            Y(i-1 + (j-2)*n,1)=subs(subs(f,x,A(i)),y,B(j)) -  (q/(h^2)+r/(2*h))*subs(subs(fy1,x,A(i)),y,B(j+1));
        end
        
    else
        
        if i==2 
            Y(i-1 + (j-2)*n,1)=subs(subs(f,x,A(i)),y,B(j)) - (p/(h^2)-r/(2*h))*subs(subs(fy0,x,A(i-1)),y,B(j)) ;
        elseif i==n+1
            Y(i-1 + (j-2)*n,1)=subs(subs(f,x,A(i)),y,B(j)) - (p/(h^2)+r/(2*h))*subs(subs(fx1,x,A(i+1)),y,B(j)) ;
        else
            Y(i-1 + (j-2)*n,1)=subs(subs(f,x,A(i)),y,B(j)) ;
        end
        
    end
    
end
end

%disp(X)
%disp(Y)

Z=inv(X)*Y;
%disp(Z)

U=zeros(n^2,1);

for j=2:n+1
for i=2:n+1
   U(i-1+(j-2)*n,1) = exp(A(i)+B(j));
end 
end

%disp(U)

for i=1:n^2
    error(i,1)=abs(Z(i,1)-U(i,1));

end

%disp(error)

disp(max(error))
   
    
 %x0 === subs(subs(fx0,x,A(i)),y,B(j-1))
 %y0 === subs(subs(fy0,x,A(i-1)),y,B(j))
 %x1 === subs(subs(fx1,x,A(i+1)),y,B(j))
 %y1 === subs(subs(fy1,x,A(i)),y,B(j+1))
    
    
    
    
    
    
    
    
