%%
%Elemento      i    j

%    1         1    2
%    2         2    3
%    3         2    3
%    4         3    4


clc
clear all
close all
format shortE
syms u1 u2 u3 u4 f1 f2 f3 f4;

k1=1200e3;
k2=2350e3;
E=200e9;
A=0.013;
L=0.65;
p=38e3;

K1=SpringElementStiffness(k1);
K2=LinearBarElementStiffness(E,A,L);
K3=LinearBarElementStiffness(E,2*A,L);
K4=SpringElementStiffness(k2);

K=zeros(4,4);

K= SpringAssemble(K,K1,1,2);
K= LinearBarAssemble(K,K2,2,3);
K= LinearBarAssemble(K,K3,2,3);
K= SpringAssemble(K,K4,3,4);

K

u1=0;
f2=0;
f3=0;
f4=p;


U=[u1;u2;u3;u4];
F=[f1;f2;f3;f4];

[u2,u3,u4,f1]=solve(K*U==F,u2,u3,u4,f1);

U2=vpa(u2,5);
U3=vpa(u3,5);
U4=vpa(u4,5);
F1=vpa(f1,5)

U1=u1;
u=[U1;U2;U3;U4];

u
F

ua=[u2;u3];
st1=LinearBarElementStresses(K2,ua,A);
ST1=vpa(st1,5)

ub=[u2;u3];
st2=LinearBarElementStresses(K3,ub,2*A);
ST2=vpa(st2,5)


R1=(K1*[u1;u2])
R2=(K4*[u3;u4])


