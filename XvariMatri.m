%clc
syms x y

a = input('Enter the value of a: ');
h = input('Enter the value of h: ');
f = input('Enter the f(x,y): ');
c = input('Enter the C(x,y): ');
d = input('Enter the D(x,y): ');
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

A = 1 + (h*h/12)*(c*c - (subs(c,x,x+h) - subs(c,x,x-h))/h  );
B = 1 + (h*h/12)*(d*d - (subs(c,y,y+h) - subs(c,y,y-h))/h  );
C = c + (h*h/12)*( c1*subs(c,x,x+h) + c2*c + c3*subs(c,x,x-h) + c4*subs(c,y,y+h) + c5*subs(c,y,y-h));
D = d + (h*h/12)*( c1*subs(d,x,x+h) + c2*d + c3*subs(d,x,x-h) + c4*subs(d,y,y+h) + c5*subs(d,y,y-h));
F = f + (h*h/12)*( c1*subs(f,x,x+h) + c2*f + c3*subs(f,x,x-h) + c4*subs(f,y,y+h) + c5*subs(f,y,y-h));
G = (subs(d,x,x+h) - subs(d,x,x-h) + subs(c,y,y+h) - subs(c,y,y-h))/(2*h) - c*d;

a1 = -1/(6*h*h) + c/(12*h) + d/(12*h) + G/24;           %i+1,j+1
a2 = -B/(h*h) + D/(2*h) + 1/(3*h*h) - d/(6*h);           %i,j+1
a3 = -1/(6*h*h) - c/(12*h) + d/(12*h) - G/24;            %i-1,j+1
a4 = -A/(h*h) + C/(2*h) + 1/(3*h*h) - c/(6*h);            %i+1,j
a5 = (2*A)/(h*h) + (2*B)/(h*h) - 2/(3*h*h);             %i,j
a6 = -A/(h*h) - C/(2*h) + 1/(3*h*h) + c/(6*h);          %i-1,j
a7 = -1/(6*h*h) + c/(12*h) - d/(12*h) - G/24;            %i+1,j-1    
a8 = -B/(h*h) - D/(2*h) + 1/(3*h*h) + d/(6*h);           %i,j-1
a9 = -1/(6*h*h) - c/(12*h) -d/(12*h) + G/24;                %i-1,j-1;             

f1=1/12 - (c*h)/24;
f2=2/3;
f3=1/12 + (c*h)/24;
f4=1/12 - (d*h)/24;
f5=1/12 + (d*h)/24;

A = linspace(0,a,n+2);
B = linspace(0,a,n+2);

X = zeros(n^2);
Ua=zeros(n^2,1);
error=zeros(n^2,1);
Y=zeros(n^2,1);

for k = 1:(n^2)
    i1=floor((k-1)/n) + 1;
    if (rem(k,n)==0)
        j1=n;
    else 
        j1=rem(k,n);
    end
    
    p = i1 + 1;
    q = j1 + 1;
    
    X(k,k)= subs(subs(a5,x,A(p)),y,B(q));
    
    if k< n^2-n +1
        X(k,k+n) = subs(subs(a4,x,A(p)),y,B(q)) ;
        X(k+n,k) = subs(subs(a6,x,A(p)),y,B(q)) ;  
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
            X(k,k+n-1) = subs(subs(a7,x,A(p)),y,B(q)) ;
            X(k+n-1,k) = subs(subs(a3,x,A(p)),y,B(q)) ;
        end
    end    
    if k< (n^2)
        if (mod(k,n)==0)
            X(k,k+1) = 0 ;
            X(k+1,k) = 0 ;
        else
            X(k,k+1) = subs(subs(a2,x,A(p)),y,B(q)) ;
            X(k+1,k) = subs(subs(a8,x,A(p)),y,B(q)) ;   
        end
    end
end

%disp(x)