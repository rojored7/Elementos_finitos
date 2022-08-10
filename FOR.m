clc
clear all
close all

global K;
global kl;
global i;
global j;
global k;
global l;

i=4;
j=5;
k=6;
l=7;
k1=kl;

K=zeros(10,10);
pos = [i,j,k,l];
for z = 1:1:length(pos)
        K(pos(z),i) = K(pos(z),i) + k1(z,1);
        K(pos(z),j) = K(pos(z),j) + k1(z,2);
        K(pos(z),k) = K(pos(z),k) + k1(z,3);
        K(pos(z),l) = K(pos(z),l) + k1(z,2);
end


  K
