% ELEMENTO	I	J	TETHA
% 1         1	4	-36,8699
% 2         2	4	0
% 3         3	4	45
% 4     	4	5	0
clc;
clear all;
syms u1 u2 u3 u4 u5 fx1 fx2 fx3 fx4 fx5 fy1 fy2 fy3 fy4 fy5 
A=0.01;
E=70e9;
P=10e3;
K1=3e6;
t1=-36.8699;
t2=0;
t3=45;
t4=0;


k1=PlaneTrussElementStiffness(E,A,(sqrt(3^2+4^2)),t1);
k2=PlaneTrussElementStiffness(E,A,4,t2);
k3=PlaneTrussElementStiffness(E,A,(sqrt(4^2+4^2)),t3);
k4=K1*[     cosd(t4)^2          cosd(t4)*sind(t4)       -cosd(t4)^2         -cosd(t4)*sind(t4);   % 1  X1
              cosd(t4)*sind(t4)          sind(t4)^2        -cosd(t4)*sind(t4)       -sind(t4)^2   ;   % 2  Y1
                -cosd(t4)^2        -cosd(t4)*sind(t4)         cosd(t4)^2         cosd(t4)*sind(t4);   % 3  X2
             -cosd(t4)*sind(t4)        -sind(t4)^2          cosd(t4)*sind(t4)        sind(t4)^2   ] % 4  Y2

         

K=zeros(10,10);

K=PlaneTrussAssemble(K,k1,1,4);
K=PlaneTrussAssemble(K,k2,2,4);
K=PlaneTrussAssemble(K,k3,3,4);
K=PlaneTrussAssemble(K,k4,4,5);
K
%%

U=[0;0;0;0;0;0;u1;0;u2;0];
F=[fx1;fy1;fx2;fy2;fx3;fy3;fx4;0;P;0];

f=[0;0;P]
k=K(7:9,7:9)
U=k\f
U=[0;0;0;0;0;0;U(1);U(2);U(3);0]

F=K*U

%%

 u1=[U(1);U(2);U(7);U(8)];
 u2=[U(3);U(4);U(7);U(8)];
 u3=[U(5);U(6);U(7);U(8)];
 
 
 st1=PlaneTrussElementStress(E,(sqrt(3^2+4^2)),t1,u1);
 st2=PlaneTrussElementStress(E,4,t2,u2);
 st3=PlaneTrussElementStress(E,(sqrt(4^2+4^2)),t3,u3);
 

 St1=vpa(st1,5)
 St2=vpa(st2,5)
 St3=vpa(st3,5)
 
 %%
 u4=[U(7);U(8);U(9);U(10)];
 F1=SpringElementForces(K1,u4)