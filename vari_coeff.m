%clc
clear all
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

