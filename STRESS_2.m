%%
%Elemento      i    j

%    1         1    2
%    2         2    3

clc
clear all
close all
format shortE
syms u1 u2 u3 f1 f2 f3 p;

E = 210e9;
A = 0.003;
L1 = 1.5;
L2 = 1;
d = 0.002;
p = 10e3;

k1=LinearBarElementStiffness(E,A,L1);
k2=LinearBarElementStiffness(E,A,L2);
K=zeros(3,3);
K=LinearBarAssemble(K,k1,1,2);
K=LinearBarAssemble(K,k2,2,3);
K

u1=0;
u3=0;
f2=-p;

U=[u1;u2;u3];
F=[f1;f2;f3];

[u2,f1,f3]=solve(K*U==F,u2,f1,f3);

U2=vpa(u2,5)
F1=vpa(f1,5)
F2=vpa(f3,5)

%%

ua=[u1;u2];
st1=LinearBarElementStresses(k1,ua,A);
ST1=vpa(st1,5)

ub=[u2;u3];
st2=LinearBarElementStresses(k2,ub,A);
ST2=vpa(st2,5)

