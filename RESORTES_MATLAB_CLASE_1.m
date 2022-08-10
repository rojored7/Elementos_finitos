clc
clear all
close all

format shortE
syms u1 u2 u3 f1 f2 f3;
k1= 100;
k2=200;
P=15;
A=SpringElementStiffness(k1);
B=SpringElementStiffness(k2);
K=zeros(3,3);
K= SpringAssemble(K,A,1,2);
K= SpringAssemble(K,B,2,3);
U=[u1;u2;u3];
F=[f1;0;f3];

f=[0;P];
u=[u2;u3];
K
KS=K(2:3,2:3);

D=KS\f

U=[0;D(1);D(2)]
S=K*U;
