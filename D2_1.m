%%
%Elemento               i               j                theta
%   1                   1               2                  0                  
%   2                   1               3                 45   
%   3                   1               4                 90  



clc
clear all
close all


th1=90;
th2=45;
th3=0;
E=30e6;
A=2;
L1=120;
L2=120*sqrt(200);
L3=120;
    
%                    1                         2                    3                          4 
%                   X1                        Y1                    X2                         Y2
                   
K1=(E*A/L1)*[     cosd(th1)          cosd(th1)*sind(th1)       -cosd(th1)^2         -cosd(th1)*sind(th1);   % 1  X1
              cosd(th1)*sind(th1)          sind(th1)^2        -cosd(th1)*sind(th1)       -sind(th1)^2   ;   % 2  Y1
                -cosd(th1)^2        -cosd(th1)*sind(th1)         cosd(th1)^2         cosd(th1)*sind(th1);   % 3  X2
             -cosd(th1)*sind(th1)        -sind(th1)^2          cosd(th1)*sind(th1)        sind(th1)^2   ;]; % 4  Y2
         

%                    1                         2                    5                          6 
%                   X1                        Y1                    X3                         Y3                    

K2=(E*A/L2)*[     cosd(th2)          cosd(th2)*sind(th2)       -cosd(th2)^2         -cosd(th2)*sind(th2);   % 1  X1
              cosd(th2)*sind(th2)          sind(th2)^2        -cosd(th2)*sind(th2)       -sind(th2)^2   ;   % 2  Y1
                -cosd(th2)^2        -cosd(th2)*sind(th2)         cosd(th2)^2         cosd(th2)*sind(th2);   % 5  X3
             -cosd(th2)*sind(th2)        -sind(th2)^2          cosd(th2)*sind(th2)        sind(th2)^2   ;]; % 6  Y3

         
%                    1                         2                    7                          8 
%                   X1                        Y1                    X4                         Y4

K3=(E*A/L3)*[     cosd(th3)          cosd(th3)*sind(th3)       -cosd(th3)^2         -cosd(th3)*sind(th3);   % 1  X1
              cosd(th3)*sind(th3)          sind(th3)^2        -cosd(th3)*sind(th3)       -sind(th3)^2   ;   % 2  Y1
                -cosd(th3)^2        -cosd(th3)*sind(th3)         cosd(th3)^2         cosd(th3)*sind(th3);   % 7  X4
             -cosd(th3)*sind(th3)        -sind(th3)^2          cosd(th3)*sind(th3)        sind(th3)^2   ;]; % 8  Y4        
         
         
%%
K=zeros(8,8);

pos1 = [1,2,3,4];
pos2 = [1,2,5,6];
pos3 = [1,2,7,8];
K1
K2
K3
for z = 1:1:length(pos1)
         K(pos1(z),1) = K(pos1(z),1) + K1(z,1);
         K(pos1(z),2) = K(pos1(z),2) + K1(z,2);
         K(pos1(z),3) = K(pos1(z),3) + K1(z,3);
         K(pos1(z),4) = K(pos1(z),4) + K1(z,4);
end

for z = 1:1:length(pos2)
        K(pos2(z),1) = K(pos2(z),1) + K2(z,1);
        K(pos2(z),2) = K(pos2(z),2) + K2(z,2);
        K(pos2(z),5) = K(pos2(z),5) + K2(z,3);
        K(pos2(z),6) = K(pos2(z),6) + K2(z,4);
end


for z = 1:1:length(pos3)
        K(pos3(z),1) = K(pos3(z),1) + K3(z,1);
        K(pos3(z),2) = K(pos3(z),2) + K3(z,2);
        K(pos3(z),7) = K(pos3(z),7) + K3(z,3);
        K(pos3(z),8) = K(pos3(z),8) + K3(z,4);
end


K