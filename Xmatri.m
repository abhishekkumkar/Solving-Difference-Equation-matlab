%clc
syms x y

a = input('Enter the value of a: ');
c = input('Enter the value of c: ');
d = input('Enter the value of d: ');
h = input('Enter the value of h: ');
n = (a/h) - 1;

a1=-(c*d)/24-1/6/h/h+c/12/h+d/12/h;               %i+1,j+1
a2=-1/h/h-d^2/12+d/2/h+1/3/h/h-d/6/h;             %i,j+1
a3= (c*d)/24-1/6/h/h-c/12/h+d/12/h;               %i-1,j+1
a4=-1/h/h-c^2/12+c/2/h+1/3/h/h-c/6/h;             %i+1,j
a5= 4/h/h+c^2/6+d^2/6-2/3/h/h;                    %i,j
a6=-1/h/h-c^2/12-c/2/h+1/3/h/h+c/6/h;             %i-1,j
a7= (c*d)/24-1/6/h/h+c/12/h-d/12/h;               %i+1,j-1    
a8=-1/h/h-d^2/12-d/2/h+1/3/h/h+d/6/h;             %i,j-1
a9=-(c*d)/24-1/6/h/h-c/12/h-d/12/h;               %i-1,j-1
        
A = linspace(0,a,n+2);
B = linspace(0,a,n+2);

X = a5*eye(n^2);

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


X