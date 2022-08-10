clc
clear all
close all
format shortE
syms u1 u2 u3 f1 f2 f3 p;

E = 210e9;
A = 0.003;
L = 0.6;
h1 =0.025;
h2 =0.0323;
h3 =0.04;
h4 =0.047;
h5 =0.0543;
A1 =pi*(h1/2)^2;
A2 =pi*(h2/2)^2;
A3 =pi*(h3/2)^2;
A4 =pi*(h4/2)^2;
A5 =pi*(h5/2)^2;
p = 18e3;

k1=LinearBarElementStiffness(E,A1,L);
k2=LinearBarElementStiffness(E,A2,L);
k3=LinearBarElementStiffness(E,A3,L);
k4=LinearBarElementStiffness(E,A4,L);
k5=LinearBarElementStiffness(E,A5,L);

K=zeros(6,6);

K=LinearBarAssemble(K,k1,1,2);
K=LinearBarAssemble(K,k2,2,3);
K=LinearBarAssemble(K,k3,3,4);
K=LinearBarAssemble(K,k4,4,5);
K=LinearBarAssemble(K,k5,5,6);

K

