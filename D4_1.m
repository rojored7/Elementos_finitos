%%
%PlaneTrussElementLength (x1,y1,x2,y2)
%PlaneTrussElementStiffness (E, A, L1, theta1)
%PlaneTrussAssemble(K,k1,1,2)
%PlaneTrussInclinedSupport(T,i,Theta)
%PlaneTrussElementStress(E,L1,theta1,u1)

%%
%Elemento               i               j                theta   
%   1                   1               4                 0                  
%   2                   1               3                 41.18     
%   3                   1               2                 90  
%   4                   2               3                 0                  
%   5                   2               4                 -41.18     
%   6                   3               4                 270
%%
E=70e9;
A=0.004;
L1=3.5;
L2=4;
L3=sqrt((3.5^2)+(4^2));
Th1=0;
Th2=41.18;
Th3=90;
Th4=0;
Th5=-41.18;
Th6=270;
i=4;

%%
k1=PlaneTrussElementStiffness(E,A,L2,Th1);
k2=PlaneTrussElementStiffness(E,A,L3,Th2);
k3=PlaneTrussElementStiffness(E,A,L1,Th3);
k4=PlaneTrussElementStiffness(E,A,L2,Th4);
k5=PlaneTrussElementStiffness(E,A,L3,Th5);
k6=PlaneTrussElementStiffness(E,A,L1,Th6);

K=zeros(8,8);

K=PlaneTrussAssemble(K,k1,1,4);
K=PlaneTrussAssemble(K,k2,1,3);
K=PlaneTrussAssemble(K,k3,1,2);
K=PlaneTrussAssemble(K,k4,2,3);
K=PlaneTrussAssemble(K,k5,2,4);
K=PlaneTrussAssemble(K,k6,3,4);
K


%%

T=eye(8,8);
t1=45;
T=PlaneTrussInclinedSupport(T,i,t1);

t=T';

Tf=T*K*t;
Tf

%%
 dx1=0;
 dy1=0;
 dy4=0;
 l=Tf(3:7,3:7);
 f=[0;0;30000;0;0];
 A=l\f;
 u=[0;0;A(1);A(2);A(3);A(4);A(5);0];
 f=Tf*u
 u
 
 %%
 
 u1=[u(1);u(2);u(7);u(8)];
 u2=[u(1);u(2);u(5);u(6)];
 u3=[u(1);u(2);u(3);u(4)];
 u4=[u(3);u(4);u(5);u(6)];
 u5=[u(3);u(4);u(7);u(8)];
 u6=[u(5);u(6);u(7);u(8)];
 
 st1=PlaneTrussElementStress(E,L2,Th1,u1);
 st2=PlaneTrussElementStress(E,L3,Th2,u2);
 st3=PlaneTrussElementStress(E,L1,Th3,u3);
 st4=PlaneTrussElementStress(E,L2,Th4,u4);
 st5=PlaneTrussElementStress(E,L3,Th5,u5);
 st6=PlaneTrussElementStress(E,L1,Th6,u6);

 St1=vpa(st1,5)
 St2=vpa(st2,5)
 St3=vpa(st3,5)
 St4=vpa(st4,5)
 St5=vpa(st5,5)
 St6=vpa(st6,5)
