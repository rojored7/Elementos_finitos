clc
clear all
format shortE
syms u1 u2 u3 u4 u5 fa fb f3 f4 f5;
k=120;
P=20;
R1=SpringElementStiffness(k);
R2=SpringElementStiffness(k);
R3=SpringElementStiffness(k);
R4=SpringElementStiffness(k);
R5=SpringElementStiffness(k);
R6=SpringElementStiffness(k);
K=zeros(5,5);
K=SpringAssemble(K,R1,1,3);
K=SpringAssemble(K,R2,3,4);
K=SpringAssemble(K,R3,3,5);
K=SpringAssemble(K,R4,3,5);
K=SpringAssemble(K,R5,5,4);
K=SpringAssemble(K,R6,4,2);
U=[u1;u2;u3;u4;u5];
F=[fa;fb;0;0;P];
K
KS=K(3:5,3:5)
f=[0;0;P];
D=KS\f

U=[0;0;D(1);D(2);D(3)]
S=K*U

