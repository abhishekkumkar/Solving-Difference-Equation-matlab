% The code is for the HOC schem for convection diffusion model
% The example we have considered have diffucion coefficients (-1, -1)
%% The code begines
clc
clear all
close all
%% Convection coefficients
c = 1;
d = 1;
%% Boundaries and grid size 
a= 0;
b= 1;
h=1/4;
N=(b-a)/h;

c1=-1/h/h-c^2/12+c/2/h+1/3/h/h-c/6/h;           % coefficient of (i+1,j)up
c2= 4/h/h+c^2/6+d^2/6-2/3/h/h;                  % coefficient of (i,j)centre
c3=-1/h/h-c^2/12-c/2/h+1/3/h/h+c/6/h;           % coefficient of (i-1,j)
c4=-1/h/h-d^2/12-d/2/h+1/3/h/h+d/6/h;           % coefficient of (i,j-1)left
c5=-1/h/h-d^2/12+d/2/h+1/3/h/h-d/6/h;           % coefficient of (i,j+1)right
c6=-(c*d)/24-1/6/h/h+c/12/h+d/12/h;             % coefficient of (i+1,j+1)up-right
c7= (c*d)/24-1/6/h/h-c/12/h+d/12/h;             % coefficient of (i-1,j+1)
c8= (c*d)/24-1/6/h/h+c/12/h-d/12/h;             % coefficient of (i+1,j-1)up-left
c9=-(c*d)/24-1/6/h/h-c/12/h-d/12/h;             % coefficient of (i-1,j-1)

a0=eye(N-1)*c2;
a0=a0+diag(c4*ones(N-2,1),-1);
a0=a0+diag(c5*ones(N-2,1),1);
aoff=eye(N-1);
b0=sparse((N-1)*(N-1),(N-1)*(N-1));
    
for i=1:N-1
    b0((i-1)*(N-1)+1:(i-1)*(N-1)+(N-1),(i-1)*(N-1)+1:(i-1)*(N-1)+(N-1)) = a0;
end
full(b0);


a10=eye(N-1)*c1;
a10=a10+diag(c8*ones(N-2,1),-1);
a10=a10+diag(c6*ones(N-2,1),1);
a1off=eye(N-1);
b10=sparse((N-1)*(N-1),(N-1)*(N-1));
    
 for i=2:N-1
        b10((i-2)*(N-1)+1:(i-2)*(N-1)+(N-1),(i-1)*(N-1)+1:(i-1)*(N-1)+(N-1))=a10;
 end
full(b10);

a20=eye(N-1)*c3;
a20=a20+diag(c7*ones(N-2,1),1);
a20=a20+diag(c9*ones(N-2,1),-1);
a20off=eye(N-1);
b20=sparse((N-1)*(N-1),(N-1)*(N-1));
    
for i=2:N-1
    b20((i-1)*(N-1)+1:(i-1)*(N-1)+(N-1),(i-2)*(N-1)+1:(i-2)*(N-1)+(N-1))=a20;
end
full(b20);
C=full(b0+b10+b20);

%% coefficient Matrix %............................................................  
   
x=linspace (a, b, N+1);
y=linspace (a, b, N+1);

%% RHS & Exact Soln

F = zeros(N+1);
for i=1:N+1
    for j=1:N+1
        f1 = 6*(exp(-x(i)^6 - y(j)^6));
        f2 = x(i)^5 + 6*x(i)^10 + y(j)^5+6*y(j)^10-5*y(j)^4-5*x(i)^4;
        H(i,j) =-f1*f2;
        U(i,j) = exp(-x(i)^6 - y(j)^6);
    end
end

%H = F              % Now we are storring F in H variable
 
%% Boundary Conditions 

for i=1:N+1
    for j = 1:N+1 
        u(i, 1)     = exp(-x(i)^6   - y(1)^6);      % y = a (j = 1)
        u(N+1, j)   = exp(-x(N+1)^6 - y(j)^6);      % x = b (i = N+1)
        u(1, j)     = exp(-x(1)^6   - y(j)^6);      % x = a (i = 1)
        u(i, N+1)   = exp(-x(i)^6   - y(N+1)^6);    % y = b (j = N+1)
    end
end
 FD = norm(H-F,2);   % Before using F-values check the difference must be zero
 for i= 2: N
     for j = 2:N
         F(i, j) = (2/3)*H(i,j)+(1/12-c*h/24)*H(i+1,j)+(1/12+c*h/24)*H(i-1,j)+(1/12-d*h/24)*H(i,j+1)+(1/12+d*h/24)*H(i,j-1);
     end
 end
 FD1 = norm(H-F,2);   % Before using F-values check the difference must be zero

 %.....................
 






 %.............................
for i = 2:N-1
    for j = 2:N
        if j==N 
             F(i, j) = F(i,j)-c3*u(i-1,j)-c5*u(i,j+1)-c9*u(i-1,j-1)-c6*u(i+1,j+1)-c7*u(i-1,j+1);
        else
             F(i, j) = F(i, j)-c3*u(i-1,j)-c4*u(i,j-1)-c7*u(i-1,j+1)-c8*u(i+1,j-1)-c9*u(i-1,j-1);
        end
    end 
end

 FD2 = norm(H-F,2);   % Before using F-values check the difference must be zero

for i=N
   for j = 2:N
        if j==N
             F(i, j) = F(i, j)-c1*u(i+1,j)-c5*u(i,j+1)-c8*u(i+1,j-1)-c7*u(i-1,j+1)-c6*u(i+1,j+1);
        else
             F(i, j) = F(i, j)-c1*u(i+1,j)-c4*u(i,j-1)-c6*u(i+1,j+1)-c9*u(i-1,j-1)-c8*u(i+1,j-1);
        end
   end
end

F(:,1) = [];
F(:,end) = [];
F(1,:) = [];
F(end, :) = [];
U(:,1) = [];
U(:, end) = [];
U(1,:) = [];
U(end, :) = [];



for i = 1:N-1
    for j = 1:N-1
        k = (N-1)*(j-1)+i;
        U1(k, 1) = U(i, j);
    end
end

G=F;

%...........................

for i = 1:N-1
    for j = 1:N-1
        k = (N-1)*(j-1)+i;
        b(k, 1) = G(i, j);
    end
end

x11 = C\b;

for i = 1 : (N-1)*(N-1)
    err(i) = x11(i)-U1(i);
end
Lin1 = norm(err, inf)






