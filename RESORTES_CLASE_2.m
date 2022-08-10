clc
close all
clear all

k1=200;
k2=200;
k3=200;
k4=200;

syms P;

E1=SpringElementStiffness(k1);
E2=SpringElementStiffness(k2);
E3=SpringElementStiffness(k3);
E4=SpringElementStiffness(k4);

K=zeros(5,5); 
K=SpringAssemble(K,E1,1,2);
K=SpringAssemble(K,E2,2,3);
K=SpringAssemble(K,E3,3,4);
K=SpringAssemble(K,E4,4,5);

K


U1=0;
U5=0.02;
syms F1 P U2 U3 U4;
U=[U1;U2;U3;U4;U5];
Ue1=[U1;U2];
Ue2=[U2;U3];
Ue3=[U3;U4];
Ue4=[U4;U5];
F=[F1;0;0;0;P];
A=200*U1-200*U2==F1;
B=-200*U1+400*U2-200*U3==0;
C=-200*U2+400*U3-200*U4==0;
D=-200*U3+400*U4-200*U5==0;
E=-200*U4+200*U5==P;

[F1,P,U2,U3,U4]=solve(A,B,C,D,E,F1,P,U2,U3,U4)
Fe1=eval(E1*Ue1)
Fe2=eval(E2*Ue2)
Fe3=eval(E3*Ue3)
Fe4=eval(E4*Ue4)
U