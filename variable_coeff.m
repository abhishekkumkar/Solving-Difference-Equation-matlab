%clc
syms x y

a = input('Enter the value of a: ');
h = input('Enter the value of h: ');
c = input('Enter the C(x,y): ');
d = input('Enter the D(x,y): ');
f = input('Enter the f(x,y): ');
fx0 = input('Enter the boundary function on X axis: ');
fy0 = input('Enter the boundary function on Y axis: ');
fx1 = input('Enter the boundary function on X=a: ');
fy1 = input('Enter the boundary function on Y=a: ');
W=input('Enter the exact solution: ');
n = (a/h) - 1;

c1 = 1/(h*h) - c/(2*h) ;
c2 = -4/(h*h) ;
c3 = 1/(h*h) + c/(2*h) ;
c4 = 1/(h*h) - d/(2*h) ;
c5 = 1/(h*h) + d/(2*h) ;

M = 1 + (h*h/12)*(c*c - subs(c,x,x+h)/h + subs(c,x,x-h)/h);
N = 1 + (h*h/12)*(d*d - subs(d,y,y+h)/h + subs(d,y,y-h)/h);
C = c + (h*h/12)*( c1*subs(c,x,x+h) + c2*c + c3*subs(c,x,x-h) + c4*subs(c,y,y+h) + c5*subs(c,y,y-h));
D = d + (h*h/12)*( c1*subs(d,x,x+h) + c2*d + c3*subs(d,x,x-h) + c4*subs(d,y,y+h) + c5*subs(d,y,y-h));
F = f + (h*h/12)*( c1*subs(f,x,x+h) + c2*f + c3*subs(f,x,x-h) + c4*subs(f,y,y+h) + c5*subs(f,y,y-h));
G = subs(d,x,x+h)/(2*h) - subs(d,x,x-h)/(2*h) + subs(c,y,y+h)/(2*h) - subs(c,y,y-h)/(2*h) - c*d;

a1 = -1/(6*h*h) + c/(12*h) + d/(12*h) + G/24;           %i+1,j+1
a2 = -N/(h*h) + D/(2*h) + 1/(3*h*h) - d/(6*h);          %i,j+1
a3 = -1/(6*h*h) - c/(12*h) + d/(12*h) - G/24;           %i-1,j+1
a4 = -M/(h*h) + C/(2*h) + 1/(3*h*h) - c/(6*h);          %i+1,j
a5 = (2*M)/(h*h) + (2*N)/(h*h) - 2/(3*h*h);             %i,j
a6 = -M/(h*h) - C/(2*h) + 1/(3*h*h) + c/(6*h);          %i-1,j
a7 = -1/(6*h*h) + c/(12*h) - d/(12*h) - G/24;           %i+1,j-1    
a8 = -N/(h*h) - D/(2*h) + 1/(3*h*h) + d/(6*h);          %i,j-1
a9 = -1/(6*h*h) - c/(12*h) - d/(12*h) + G/24;           %i-1,j-1;   

A = linspace(0,a,n+2);
B = linspace(0,a,n+2);

X = zeros(n^2);
Ua = zeros(n^2,1);
error = zeros(n^2,1);
Y = zeros(n^2,1);
for k = 1:(n^2)
    i = floor((k-1)/n) + 1;
    if (rem(k,n)==0)
        j = n;
    else 
        j = rem(k,n);
    end
    p = j + 1;
    q = i + 1;
    
    X(k,k) = subs(subs(a5,x,A(p)),y,B(q));
    if k< n^2-n+1
        X(k,k+n) = subs(subs(a2,x,A(p)),y,B(q)) ;
        X(k+n,k) = subs(subs(a8,x,A(p)),y,B(q)) ;  
    end
    if k< (n^2)-n
        if (mod(k,n)==0)
            X(k,k+n+1) = 0;
            X(k+n+1,k) = 0;
        else
            X(k,k+n+1) = subs(subs(a1,x,A(p)),y,B(q)) ;
            X(k+n+1,k) = subs(subs(a9,x,A(p)),y,B(q)) ;
        end
    end
    if k< (n^2) - n + 2
        if mod(k+n-1,n)==0
            X(k,k+n-1) = 0;
            X(k+n-1,k) = 0;
        else
            X(k,k+n-1) = subs(subs(a3,x,A(p)),y,B(q)) ;
            X(k+n-1,k) = subs(subs(a7,x,A(p)),y,B(q)) ;
        end
    end    
    if k< (n^2)
        if (mod(k,n)==0)
            X(k,k+1) = 0 ;
            X(k+1,k) = 0 ;
        else
            X(k,k+1) = subs(subs(a4,x,A(p)),y,B(q)) ;
            X(k+1,k) = subs(subs(a6,x,A(p)),y,B(q)) ;   
        end
    end
end

%disp(x)

for j=2:n+1     %y
for i=2:n+1     %x
    
    a11=subs(subs(a1,x,A(i)),y,B(j));
    a22=subs(subs(a2,x,A(i)),y,B(j));
    a33=subs(subs(a3,x,A(i)),y,B(j));
    a44=subs(subs(a4,x,A(i)),y,B(j));
    a55=subs(subs(a5,x,A(i)),y,B(j));
    a66=subs(subs(a6,x,A(i)),y,B(j));
    a77=subs(subs(a7,x,A(i)),y,B(j));
    a88=subs(subs(a8,x,A(i)),y,B(j));
    a99=subs(subs(a9,x,A(i)),y,B(j));
    
    fxy=subs(subs(F,x,A(i)),y,B(j));
    Ua(j-1+(i-2)*n,1) = subs(subs(W,x,A(i)),y,B(j));
    
    if (j==2)
        if (i==2)
            Y(i-1 + (j-2)*n,1) = fxy - a33*subs(subs(fx0,x,A(i+1)),y,B(j-1)) - a66*subs(subs(fx0,x,A(i)),y,B(j-1)) - a99*subs(subs(fx0,x,A(i-1)),y,B(j-1)) - a88*subs(subs(fy0,x,A(i-1)),y,B(j)) - a77*subs(subs(fy0,x,A(i-1)),y,B(j+1));
        elseif (i==n+1)
            Y(i-1 + (j-2)*n,1) = fxy - a33*subs(subs(fx0,x,A(i+1)),y,B(j-1)) - a66*subs(subs(fx0,x,A(i)),y,B(j-1)) - a99*subs(subs(fx0,x,A(i-1)),y,B(j-1)) - a11*subs(subs(fx1,x,A(i+1)),y,B(j+1)) - a22*subs(subs(fx1,x,A(i+1)),y,B(j));
        else
            Y(i-1 + (j-2)*n,1) = fxy - a33*subs(subs(fx0,x,A(i+1)),y,B(j-1)) - a66*subs(subs(fx0,x,A(i)),y,B(j-1)) - a99*subs(subs(fx0,x,A(i-1)),y,B(j-1));
        end
    elseif(j==n+1)
        if (i==2)
            Y(i-1 + (j-2)*n,1) = fxy - a11*subs(subs(fy1,x,A(i+1)),y,B(j+1)) - a44*subs(subs(fy1,x,A(i)),y,B(j+1)) - a77*subs(subs(fy1,x,A(i-1)),y,B(j+1)) - a99*subs(subs(fy0,x,A(i-1)),y,B(j-1)) - a88*subs(subs(fy0,x,A(i-1)),y,B(j));
        elseif (i==n+1)
            Y(i-1 + (j-2)*n,1) = fxy - a11*subs(subs(fy1,x,A(i+1)),y,B(j+1)) - a22*subs(subs(fx1,x,A(i+1)),y,B(j)) - a33*subs(subs(fx1,x,A(i+1)),y,B(j-1)) - a44*subs(subs(fy1,x,A(i)),y,B(j+1)) - a77*subs(subs(fy1,x,A(i-1)),y,B(j+1));
        else
            Y(i-1 + (j-2)*n,1) = fxy - a11*subs(subs(fy1,x,A(i+1)),y,B(j+1)) - a44*subs(subs(fy1,x,A(i)),y,B(j+1)) - a77*subs(subs(fy1,x,A(i-1)),y,B(j+1));
        end
    else
        if (i==2)
            Y(i-1 + (j-2)*n,1) = fxy - a77*subs(subs(fy0,x,A(i-1)),y,B(j+1)) - a88*subs(subs(fy0,x,A(i-1)),y,B(j)) - a99*subs(subs(fy0,x,A(i-1)),y,B(j-1));
        elseif (i==n+1)
            Y(i-1 + (j-2)*n,1) = fxy - a11*subs(subs(fx1,x,A(i+1)),y,B(j+1)) - a22*subs(subs(fx1,x,A(i+1)),y,B(j)) - a33*subs(subs(fx1,x,A(i+1)),y,B(j-1)) ;
        else
            Y(i-1 + (j-2)*n,1) = fxy ;
        end
   end
end
end

Z=inv(X)*Y;
%disp(Z)

for i=1:n^2
    error(i,1)=abs(Z(i,1)-Ua(i,1));
end

disp(max(error))














