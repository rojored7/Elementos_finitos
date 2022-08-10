syms k1 k2 k3 k4;
syms u1 u2 u3;
syms f1 f2 f4 f5;

K=[k4 -k4 0 0 0; -k4 k1+k2+k4 -k2 -k1 0; 0 -k2 k2+k3 0 -k3; 0 -k1 0 k1 0; 0 0 -k3 0 k3];

f=[f1; 0; f2; f4; f5];

l=K(1:3,1:3);
n=f(1:3);

r=simplify(l\n);

u=[r;0;0];

t=simplify(K*u)