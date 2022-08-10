%%
%PlaneTrussElementLength (x1,y1,x2,y2)
%PlaneTrussElementStiffness (E, A, L1, theta1)
%PlaneTrussAssemble(K,k1,1,2)
%PlaneTrussInclinedSupport(T,i,Theta)
%PlaneTrussElementStress(E,L1,theta1,u1)

%%
%Elemento               i               j                theta   
%   1                   1               2                 45                  
%   2                   2               3                 180     
%   3                   1               3                 90       

E=30e6;
A=2;
L=30;
L1=sqrt(L^2+L^2);
Th1=45;
Th2=180;
Th3=90;
i=1;

k1=PlaneTrussElementStiffness(E,A,L1,Th1);
k2=PlaneTrussElementStiffness(E,A,L ,Th2);
k3=PlaneTrussElementStiffness(E,A,L ,Th3);

K=zeros(6,6);

K=PlaneTrussAssemble(K,k1,1,2);
K=PlaneTrussAssemble(K,k2,2,3);
K=PlaneTrussAssemble(K,k3,1,3);

K
%%

T=eye(6,6);
t1=-45;
T=PlaneTrussInclinedSupport(T,i,t1);

t=T';

Tf=T*K*t;
Tf

%%

syms dy1 dy2;

dx1=0;
dx2=0;
dx3=0;
dy3=0;
fy1=0;
fy2=-2000;

U=[dy1,dy2];
F=[fy1;fy2];
l=[Tf(2,2) Tf(2,4);Tf(4,2) Tf(4,4)];
a=l\F;
u=[0;a(1);0;a(2);0;0];
f=Tf*u
u