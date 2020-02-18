%clc
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
W=input('Enter the exact solution: ');
n = (a/h) - 1;

a1=-(c*d)/24-1/6/h/h+c/12/h+d/12/h;
a2=-1/h/h-d^2/12+d/2/h+1/3/h/h-d/6/h;
a3= (c*d)/24-1/6/h/h-c/12/h+d/12/h;
a4=-1/h/h-c^2/12+c/2/h+1/3/h/h-c/6/h;
a5= 4/h/h+c^2/6+d^2/6-2/3/h/h;
a6=-1/h/h-c^2/12-c/2/h+1/3/h/h+c/6/h;           
a7= (c*d)/24-1/6/h/h+c/12/h-d/12/h;                   
a8=-1/h/h-d^2/12-d/2/h+1/3/h/h+d/6/h;           
a9=-(c*d)/24-1/6/h/h-c/12/h-d/12/h;             

f1=1/12 - (c*h)/24;
f2=2/3;
f3=1/12 + (c*h)/24;
f4=1/12 - (d*h)/24;
f5=1/12 + (d*h)/24;

A = linspace(0,a,n+2);
B = linspace(0,a,n+2);

X = a5*eye(n^2);
Ua=zeros(n^2,1);
error=zeros(n^2,1);
Y=zeros(n^2,1);

for i = 1:(n^2)
    if i< n^2-n +1
        X(i,i+n) = a4 ;
        X(i+n,i) = a6 ;  
    end
    if i< (n^2)-n
        if (mod(i,n)==0)
            X(i,i+n+1) = 0;
            X(i+n+1,i) = 0;
        else
            X(i,i+n+1) = a1;
            X(i+n+1,i) = a9;
        end
    end
    if i< (n^2) - n + 2
        if mod(i+n-1,n)==0
            X(i,i+n-1) = 0;
            X(i+n-1,i) = 0;
        else
            X(i,i+n-1) = a7;
            X(i+n-1,i) = a3;
        end
    end    
    if i< (n^2)
        if (mod(i,n)==0)
            X(i,i+1) = 0 ;
            X(i+1,i) = 0 ;
        else
            X(i,i+1) = a2 ;
            X(i+1,i) = a8 ;   
        end
    end
end

%disp(x)

for j=2:n+1     %y
for i=2:n+1     %x
    
    fxy = f1*subs(subs(f,x,A(i+1)),y,B(j)) +  f2*subs(subs(f,x,A(i)),y,B(j)) +  f3*subs(subs(f,x,A(i-1)),y,B(j)) +  f4*subs(subs(f,x,A(i)),y,B(j+1)) +  f5*subs(subs(f,x,A(i)),y,B(j-1)) ;
    Ua(i-1+(j-2)*n,1) = subs(subs(W,x,A(i)),y,B(j));
    
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
            Y(i-1 + (j-2)*n,1) = fxy - a1*subs(subs(fx1,x,A(i+1)),y,B(j+1)) - a2*subs(subs(fx1,x,A(i+1)),y,B(j)) - a3*subs(subs(fx1,x,A(i+1)),y,B(j-1)) ;
        else
            Y(i-1 + (j-2)*n,1) = fxy ;
        end
   end
end
end

Z=X\Y;
%disp(Z)

for i=1:n^2
    error(i,1)=abs(Z(i,1)-Ua(i,1));
end

disp(max(error))