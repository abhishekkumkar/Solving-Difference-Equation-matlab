clc
clear all
syms x y

a = input('Enter the value of a: ');
c = input('Enter the value of c: ');
d = input('Enter the value of d: ');
h = input('Enter the value of h: ');
f = input('Enter the f(x,y): ');
fx0 = input('Enter the boundary function on X axis: ');
fy0 = input('Enter the boundary function on Y axis: ');
fx1 = input('Enter the boundary function on X=a: ');
fy1 = input('Enter the boundary function on Y=a: ');

n = (a/h) - 1;

p = 1+ ((c*h)^2)/12;
q = 1+ ((d*h^2))/12;

a1 = -1/(6*(h^2)) + (c+d)/(12*h) - (c*d)/24 ;
a2 = (-p+(1/3))/(h^2) + (c/2-c/6)/h ;
a3 = -1/(6*(h^2)) + (c-d)/(12*h) + (c*d)/24 ;
a4 = (-q+(1/3))/(h^2) + (d/2-d/6)/h ;
a5 = (2*p+2*q-(2/3))/(h^2) -(c+d)/(h^2)  ;
a6 = (-q+(1/3))/(h^2)  + (d/6)/h ;
a7 = -1/(6*(h^2)) + (d-c)/(12*h) + (c*d)/(24*h)  ;
a8 = (-p+(1/3))/(h^2) + (c/6)/h ;
a9 = -1/(6*(h^2)) - (c+d)/(12*h) - (c*d)/24 ;

fx = diff(f,x);
fxx = diff(fx,x);
fy = diff(f,y);
fyy = diff(fy,y);
F = f + ((h^2)/12)*( fxx + fyy - c*fx - d*fy);

A = linspace(0,a,n+2);
B = linspace(0,a,n+2);

X = a5*eye(n^2);

for i = 1:(n^2) - n
    X(i,i+n) = a4 ;
    X(i+n,i) = a6 ;   
end

for i = 1 : (n^2)-n-1
    if (mod(i,n)==0)
        X(i,i+n+1) = 0;
        X(i+n+1,i) = 0;
    else
        X(i,i+n+1) = a1;
        X(i+n+1,i) = a9;
    end
end

for i = 1 : (n^2) - n + 1
    if (mod(i+n-1,n)==0)
        X(i,i+n-1) = 0;
        X(i+n-1,i) = 0;
    else
        X(i,i+n-1) = a7;
        X(i+n-1,i) = a3;
    end
end
        

for i = 1:(n^2) - 1
    if (mod(i,n)==0)
        X(i,i+1) = 0 ;
        X(i+1,i) = 0 ;
    else
        X(i,i+1) = a2 ;
        X(i+1,i) = a8 ;   
    end
end

%disp(X) % X is verified

Y=ones(n^2,1);

for j=2:n+1     %y
for i=2:n+1     %x
    
    fxy = subs(subs(F,x,A(i)),y,B(j)) ;
    
    if (j==2)
        
        if (i==2)
            Y(i-1 + (j-2)*n,1) = fxy - a3*subs(subs(fx0,x,A(i+1)),y,B(j-1)) - a6*subs(subs(fx0,x,A(i)),y,B(j-1)) - a9*subs(subs(fx0,x,A(i-1)),y,B(j-1)) - a8*subs(subs(fy0,x,A(i-1)),y,B(j)) - a7*subs(subs(fy0,x,A(i-1)),y,B(j+1));
        elseif (i==n+1)
            Y(i-1 + (j-2)*n,1) = fxy - a3*subs(subs(fx0,x,A(i+1)),y,B(j-1)) - a6*subs(subs(fx0,x,A(i)),y,B(j-1)) - a9*subs(subs(fx0,x,A(i-1)),y,B(j-1)) - a1*subs(subs(fx1,x,A(i+1)),y,B(j+1)) - a2*subs(subs(fx1,x,A(i+1)),y,B(j));
        else
            Y(i-1 + (j-2)*n,1) = fxy - a3*subs(subs(fx0,x,A(i+1)),y,B(j-1)) - a6*subs(subs(fx0,x,A(i)),y,B(j-1)) - a9*subs(subs(fx0,x,A(i-1)),y,B(j-1));
        end
        
    elseif(j==n+1)
        
        if (i==2)
            Y(i-1 + (j-2)*n,1) = fxy - a1*subs(subs(fy1,x,A(i+1)),y,B(j+1)) - a4*subs(subs(fy1,x,A(i)),y,B(j+1)) - a7*subs(subs(fy1,x,A(i-1)),y,B(j+1)) - a9*subs(subs(fy0,x,A(i-1)),y,B(j-1)) - a8*subs(subs(fy0,x,A(i-1)),y,B(j));
        elseif (i==n+1)
            Y(i-1 + (j-2)*n,1) = fxy - a1*subs(subs(fy1,x,A(i+1)),y,B(j+1)) - a2*subs(subs(fx1,x,A(i+1)),y,B(j)) - a3*subs(subs(fx1,x,A(i+1)),y,B(j-1)) - a4*subs(subs(fy1,x,A(i)),y,B(j+1)) - a7*subs(subs(fy1,x,A(i-1)),y,B(j+1));
        else
            Y(i-1 + (j-2)*n,1) = fxy - a1*subs(subs(fy1,x,A(i+1)),y,B(j+1)) - a4*subs(subs(fy1,x,A(i)),y,B(j+1)) - a7*subs(subs(fy1,x,A(i-1)),y,B(j+1));
        end
    
    else
        
        if (i==2)
            Y(i-1 + (j-2)*n,1) = fxy - a7*subs(subs(fy0,x,A(i-1)),y,B(j+1)) - a8*subs(subs(fy0,x,A(i-1)),y,B(j)) - a9*subs(subs(fy0,x,A(i-1)),y,B(j-1));
        elseif (i==n+1)
            Y(i-1 + (j-2)*n,1) = fxy - a1*subs(subs(fx1,x,A(i+1)),y,B(j+1)) - a2*subs(subs(fx1,x,A(i+1)),y,B(j)) - a3*subs(subs(fx0,x,A(i+1)),y,B(j-1)) ;
        else
            Y(i-1 + (j-2)*n,1) = fxy ;
        end
   
    end
end
end

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