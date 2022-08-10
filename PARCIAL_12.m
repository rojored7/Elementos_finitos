%%
%Elemento      i    j

%    1         1    2
%    2         2    3
%    3         2    3
%    4         2    3
%    5         3    4
%    6         3    4
%    7         4    5


clc
clear all
close all
format shortE
syms u1 u2 u3 u4 u5 f1 f2 f3 f4 f5;

E=200e9;
A1 = pi*(0.009/2)^2;
A2 = pi*(0.008/2)^2;
A3 = pi*(0.006/2)^2;
A4 = pi*(0.004/2)^2;
L1 = 0.90;
L2 = 1.5;
L3 = 0.75;
L4 = 1;
P1 = 15e3;
P2 = 25e3;


K1=LinearBarElementStiffness(E,A1,L1);
K2=LinearBarElementStiffness(E,A2,L2);
K3=LinearBarElementStiffness(E,A2,L2);
K4=LinearBarElementStiffness(E,A2,L2);
K5=LinearBarElementStiffness(E,A3,L3);
K6=LinearBarElementStiffness(E,A3,L3);
K7=LinearBarElementStiffness(E,A4,L4);

K=zeros(5,5);

K= LinearBarAssemble(K,K1,1,2);
K= LinearBarAssemble(K,K2,2,3);
K= LinearBarAssemble(K,K3,2,3);
K= LinearBarAssemble(K,K4,2,3);
K= LinearBarAssemble(K,K5,3,4);
K= LinearBarAssemble(K,K6,3,4);
K= LinearBarAssemble(K,K7,4,5);

K

U1=P1/K2(2,2)

%%
u1=0;
f2=P1;
f5=P2;
f3=0;
f4=0;




u=[u1;u2;u3;u4;u5];
F=[f1;f2;f3;f4;f5];


[u2,u3,u4,u5,f1]=solve(K*u==F,u2,u3,u4,u5,f1);

U2=vpa(u2,5);
U3=vpa(u4,5);
U4=vpa(u4,5);
U5=vpa(u5,5);
F1=vpa(f1,5);
F4=vpa(f4,5);


U=[u1;U2;U3;U4;U5]
F=[F1;f2;f3;F4;f5]

ua=[u1;u2];
st1=LinearBarElementStresses(K1,ua,A1);
ST1=vpa(st1,5)

ub=[u2;u3];
st2=LinearBarElementStresses(K2,ub,A2);
ST2=vpa(st2,5)

ub=[u2;u3];
st2=LinearBarElementStresses(K3,ub,A2);
ST2=vpa(st2,5)

ub=[u2;u3];
st2=LinearBarElementStresses(K4,ub,A2);
ST2=vpa(st2,5)

uc=[u3;u4];
st3=LinearBarElementStresses(K5,uc,A3);
ST3=vpa(st3,5)

uc=[u3;u4];
st3=LinearBarElementStresses(K6,uc,A3);
ST3=vpa(st3,5)

ud=[u4;u5];
st4=LinearBarElementStresses(K7,ud,A4);
ST4=vpa(st4,5)