%%
%Elemento               i               j                theta
%   1                   1               2                  90                  
%   2                   1               3                 45   
%   3                   1               4                 0  



clc
%clear all
%close all
syms dx1 dx2
th1=135;
th2=180;
th3=270;
E=210e9;
A=5e-4;
L1=5;
L2=10;
k3=2000000;
    
%                    1                         2                    3                          4 
%                   X1                        Y1                    X2                         Y2
                   
K1=(E*A/L1)*[     cosd(th1)^2          cosd(th1)*sind(th1)       -cosd(th1)^2         -cosd(th1)*sind(th1);   % 1  X1
              cosd(th1)*sind(th1)          sind(th1)^2        -cosd(th1)*sind(th1)       -sind(th1)^2   ;   % 2  Y1
                -cosd(th1)^2        -cosd(th1)*sind(th1)         cosd(th1)^2         cosd(th1)*sind(th1);   % 3  X2
             -cosd(th1)*sind(th1)        -sind(th1)^2          cosd(th1)*sind(th1)        sind(th1)^2   ;] % 4  Y2
         

%                    1                         2                    5                          6 
%                   X1                        Y1                    X3                         Y3                    

K2=(E*A/L2)*[     cosd(th2)^2          cosd(th2)*sind(th2)       -cosd(th2)^2         -cosd(th2)*sind(th2);   % 1  X1
              cosd(th2)*sind(th2)          sind(th2)^2        -cosd(th2)*sind(th2)       -sind(th2)^2   ;   % 2  Y1
                -cosd(th2)^2        -cosd(th2)*sind(th2)         cosd(th2)^2         cosd(th2)*sind(th2);   % 5  X3
             -cosd(th2)*sind(th2)        -sind(th2)^2          cosd(th2)*sind(th2)        sind(th2)^2   ;] % 6  Y3

         
%                    1                         2                    7                          8 
%                   X1                        Y1                    X4                         Y4

K3=k3*[     cosd(th3)^2          cosd(th3)*sind(th3)       -cosd(th3)^2         -cosd(th3)*sind(th3);   % 1  X1
              cosd(th3)*sind(th3)          sind(th3)^2        -cosd(th3)*sind(th3)       -sind(th3)^2   ;   % 2  Y1
                -cosd(th3)^2        -cosd(th3)*sind(th3)         cosd(th3)^2         cosd(th3)*sind(th3);   % 7  X4
             -cosd(th3)*sind(th3)        -sind(th3)^2          cosd(th3)*sind(th3)        sind(th3)^2   ;] % 8  Y4        
         
k1=[1,2,3,4];
k2=[1,2,5,6];
k3=[1,2,7,8];

%%
K=[K1(1,1)+K2(1,1)+K3(1,1)	K1(1,2)+K2(1,2)+K3(1,2)	K1(1,3),K1(1,4),K2(1,3),K2(1,4),K3(1,3),K3(1,4);
K1(2,1)+K2(2,1)+K3(2,1)	K1(2,2)+K2(2,2)+K3(2,2)	K1(2,3)	K1(2,4)	K2(2,3)	K2(2,4)	K3(2,3)	K3(2,4)
K1(3,1)	K1(3,2)	K1(3,3)	K1(3,4)	0	0	0	0
K1(4,1)	K1(4,2)	K1(4,3)	K1(4,4)	0	0	0	0
K2(3,1)	K2(3,2)	0	0	K2(3,3)	K2(3,4)	0	0
K2(4,1)	K2(4,2)	0	0	K2(4,3)	K2(4,4)	0	0
K3(3,1)	K3(3,2)	0	0	0	0	K3(3,3)	K3(3,4)
K3(4,1)	K3(4,2)	0	0	0	0	K3(4,3)	K3(4,4)]

F=[0;-25000];
U=[dx1;dx2];
[dx1,dx2]=solve(K(1:2,1:2)*U==F,dx1,dx2);
vpa(dx1,5)
vpa(dx2,5)

U=[dx1; dx2;0;0;0;0;0;0];
um=vpa(U,5)
Ua=[U(1);U(2);U(3);U(4)]
Ub=[U(1);U(2);U(5);U(6)]
%Uc=[U(1);U(2);U(7);U(8)]


ST1=vpa((E/L1)*[-cosd(th1) -sind(th1) cosd(th1) sind(th1)]*Ua,5)
ST2=vpa((E/L2)*[-cosd(th2) -sind(th2) cosd(th2) sind(th2)]*Ub,5)
%ST3=vpa((E/L3)*[-cosd(th3) -sind(th3) cosd(th3) sind(th3)]*Uc,5)
 

