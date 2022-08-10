clc
clear all
close all
format shortE
syms u1 u2 u3 f1 f2 f3;

E = 2e4;
A = 250;
L = 150;
d = 1.2;
P = 6e4;
k = (E*A)/L;

A=SpringElementStiffness(k);
B=SpringElementStiffness(k);
K=zeros(3,3);

K= SpringAssemble(K,A,1,2);
K= SpringAssemble(K,B,2,3);

A
B
K
U1=P/A(2,2)

%%
u1=0;
u3=1.2;
u=[u1;u2;u3];
F=[f1;P;f3];
[u2,f1,f3]=solve(K*u==F,u2,f1,f3);
vpa(u2,5)
vpa(f1,5)
vpa(f3,5)

st1=E*[-1/L 1/L]*[u1;u2];
ST1=vpa(st1,5)

st2=E*[-1/L 1/L]*[u2;u3];
ST2=vpa(st2,5)

