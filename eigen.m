clc 
clear all
A = input('Enter the Matrix A: ');
[P,D] = eig(A);
AP=A*P;
n=size(A);
O=zeros(n);
if(n==2)
    quiver(O(:,1),O(:,2),D(1,1)*P(:,1),D(2,2)*P(:,2));
    hold on
   quiver(O(:,1),O(:,2),AP(:,1),AP(:,2));
else
    quiver3(O(:,1),O(:,2),O(:,3),D(1,1)*P(:,1),D(2,2)*P(:,2),D(3,3)*P(:,3));
    hold on
    quiver3(O(:,1),O(:,2),O(:,3),AP(:,1),AP(:,2),AP(:,3));
end
    

%quiver3(O(:,1),O(:,2),O(:,3),D(1,1)*P(:,1),D(2,2)*P(:,2),D(3,3)*P(:,3));