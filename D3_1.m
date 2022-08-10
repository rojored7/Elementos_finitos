%%
%Elemento               i               j                theta   DIREC
%   1                   1               2                 0       1-2           
%   2                   2               3                 135     2-3
%   3                   1               3                 90      1-3  


% ELEMENTO              REACCIONES          DESPLA
%    1                   Y=Y1 X=0             D1
%    2                   Y=Y2 X=0             D2
%    3                   Y=Y3 X=X3            0


clc
clear all
close all

syms dx1 dx2 fx1;

th1 =0;       %BARRA
th2 =135;     %BARRA
th3 =90;      %BARRA
alpha =315;    %EJE
E=30e6;
A=2;
L1=30;
L2=30;
L3=30;

%%
%DESPLAZAR BARRAS

%                    1                         2                    3                          4 
%                   X1                        Y1                    X2                         Y2
                   
K1=(E*A/L1)*[     cosd(th1)^2          cosd(th1)*sind(th1)       -cosd(th1)^2         -cosd(th1)*sind(th1);   % 1  X1
              cosd(th1)*sind(th1)          sind(th1)^2        -cosd(th1)*sind(th1)       -sind(th1)^2   ;   % 2  Y1
                -cosd(th1)^2        -cosd(th1)*sind(th1)         cosd(th1)^2         cosd(th1)*sind(th1);   % 3  X2
             -cosd(th1)*sind(th1)        -sind(th1)^2          cosd(th1)*sind(th1)        sind(th1)^2   ;]; % 4  Y2
         
%                   3                         4                     5                          6 
%                   X2                        Y2                    X3                         Y3                    

K2=(E*A/L2)*[     cosd(th2)^2          cosd(th2)*sind(th2)       -cosd(th2)^2         -cosd(th2)*sind(th2);   % 3  X1
              cosd(th2)*sind(th2)          sind(th2)^2        -cosd(th2)*sind(th2)       -sind(th2)^2   ;   % 4  Y1
                -cosd(th2)^2        -cosd(th2)*sind(th2)         cosd(th2)^2         cosd(th2)*sind(th2);   % 5  X3
             -cosd(th2)*sind(th2)        -sind(th2)^2          cosd(th2)*sind(th2)        sind(th2)^2   ;]; % 6  Y3

         
%                    1                         2                    5                          6 
%                   X1                        Y1                    X3                         Y3

K3=(E*A/L3)*[     cosd(th3)^2          cosd(th3)*sind(th3)       -cosd(th3)^2         -cosd(th3)*sind(th3);   % 1  X1
              cosd(th3)*sind(th3)          sind(th3)^2        -cosd(th3)*sind(th3)       -sind(th3)^2   ;   % 2  Y1
                -cosd(th3)^2        -cosd(th3)*sind(th3)         cosd(th3)^2         cosd(th3)*sind(th3);   % 5  X4
             -cosd(th3)*sind(th3)        -sind(th3)^2          cosd(th3)*sind(th3)        sind(th3)^2   ;]; % 6  Y4  
         
K=zeros(6,6);


pos1 = [1,2,3,4];
pos2 = [3,4,5,6];
pos3 = [1,2,5,6];

for z = 1:1:length(pos1)
         K(pos1(z),1) = K(pos1(z),1) + K1(z,1);
         K(pos1(z),2) = K(pos1(z),2) + K1(z,2);
         K(pos1(z),3) = K(pos1(z),3) + K1(z,3);
         K(pos1(z),4) = K(pos1(z),4) + K1(z,4);
end

 for z = 1:1:length(pos2)
         K(pos2(z),3) = K(pos2(z),3) + K2(z,1);
         K(pos2(z),4) = K(pos2(z),4) + K2(z,2);
         K(pos2(z),5) = K(pos2(z),5) + K2(z,3);
         K(pos2(z),6) = K(pos2(z),6) + K2(z,4);
 end
 
 
 for z = 1:1:length(pos3)
         K(pos3(z),1) = K(pos3(z),1) + K3(z,1);
         K(pos3(z),2) = K(pos3(z),2) + K3(z,2);
         K(pos3(z),5) = K(pos3(z),5) + K3(z,3);
         K(pos3(z),6) = K(pos3(z),6) + K3(z,4);
 end


%%
% DESPLAZAR EJES

T=eye(6,6);

P1 = [1,2];
P2 = [5,6];
P3 = [3,4];


 
 t1= [ cosd(alpha) sind(alpha);
      -sind(alpha) cosd(alpha);];
 


for z = 1:1:length(P1)
    
         T(P1(z),1) = t1(z,1);% T(P2(z),5) +
         T(P1(z),2) = t1(z,2);% T(P2(z),6) +

end
t=T';
Kf=T*K*T';
K
T
t
Kf

%%
dy1=0;
dy2=0;
dx3=0;
dy3=0;
U=[dx1;dx2];
fx2=-2000;
fx1=0;
F=[fx1;fx2];
l=[K(1,1) K(1,3);K(3,1) K(3,3)]
[dx1,dx2]=solve(l*U==F,dx1,dx2);
vpa(dx1,5)
vpa(dx2,5)


U=[dx1; dx2;0;0;0;0];
um=vpa(U,5)
Ua=[U(1);U(2);U(3);U(4)]
Ub=[U(3);U(4);U(5);U(6)]
Uc=[U(1);U(2);U(5);U(6)]

ft=K*U
U
ST1=vpa((E/L1)*[-cosd(th1) -sind(th1) cosd(th1) sind(th1)]*Ua,5)
ST2=vpa((E/L2)*[-cosd(th2) -sind(th2) cosd(th2) sind(th2)]*Ub,5)
ST3=vpa((E/L3)*[-cosd(th3) -sind(th3) cosd(th3) sind(th3)]*Uc,5)



