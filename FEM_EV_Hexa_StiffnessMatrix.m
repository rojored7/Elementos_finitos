function K = FEM_EV_Hexa_StiffnessMatrix(E,v,ui,uj,uk,ul,um,un,uo,up)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento con Volumen (Hexaedro) - Matriz de Rigidez

%% Variables

% Entrada:

% E                 ::>  Módulo de Young (elasticidad) del elemento.
% v                 ::>  Coeficiente de poisson.
% ui                ::>  Posición nodo i (X,Y,Z)
% uj                ::>  Posición nodo j (X,Y,Z)
% uk                ::>  Posición nodo k (X,Y,Z)
% ul                ::>  Posición nodo l (X,Y,Z)
% um                ::>  Posición nodo m (X,Y,Z)
% un                ::>  Posición nodo n (X,Y,Z)
% uo                ::>  Posición nodo o (X,Y,Z)
% up                ::>  Posición nodo p (X,Y,Z)

% Salida:

% K                 ::>  Matriz de Rigidez


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 10
    if length(ui) == 3 && length(uj) == 3 && length(uk) == 3 && length(ul) == 3 &&...
        length(um) == 3 && length(un) == 3 && length(uo) == 3 && length(up) == 3
        ok = 1;
    end
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    Xi = ui(1);
    Yi = ui(2);
    Zi = ui(3);
    Xj = uj(1);
    Yj = uj(2);
    Zj = uj(3);
    Xk = uk(1);
    Yk = uk(2);
    Zk = uk(3);
    Xl = ul(1);
    Yl = ul(2);
    Zl = ul(3);
    Xm = um(1);
    Ym = um(2);
    Zm = um(3);
    Xn = un(1);
    Yn = un(2);
    Zn = un(3);
    Xo = uo(1);
    Yo = uo(2);
    Zo = uo(3);
    Xp = up(1);
    Yp = up(2);
    Zp = up(3);
    syms s t u;    
    N1 = (1-s)*(1-t)*(1+u)/8;
    N2 = (1-s)*(1-t)*(1-u)/8;
    N3 = (1-s)*(1+t)*(1-u)/8;
    N4 = (1-s)*(1+t)*(1+u)/8;
    N5 = (1+s)*(1-t)*(1+u)/8;
    N6 = (1+s)*(1-t)*(1-u)/8;
    N7 = (1+s)*(1+t)*(1-u)/8;
    N8 = (1+s)*(1+t)*(1+u)/8;
    X = N1*Xi + N2*Xj + N3*Xk + N4*Xl + N5*Xm + N6*Xn + N7*Xo + N8*Xp;
    Y = N1*Yi + N2*Yj + N3*Yk + N4*Yl + N5*Ym + N6*Yn + N7*Yo + N8*Yp;
    Z = N1*Zi + N2*Zj + N3*Zk + N4*Zl + N5*Zm + N6*Zn + N7*Zo + N8*Zp;
    Xs = diff(X,s);
    Xt = diff(X,t);
    Xu = diff(X,u);
    Ys = diff(Y,s);
    Yt = diff(Y,t);
    Yu = diff(Y,u);
    Zs = diff(Z,s);
    Zt = diff(Z,t);
    Zu = diff(Z,u);
    J = Xs*(Yt*Zu - Zt*Yu) - Ys*(Xt*Zu - Zt*Xu) + Zs*(Xt*Yu - Yt*Xu);
    N1s = diff(N1,s);
    N2s = diff(N2,s);
    N3s = diff(N3,s);
    N4s = diff(N4,s);
    N5s = diff(N5,s);
    N6s = diff(N6,s);
    N7s = diff(N7,s);
    N8s = diff(N8,s);
    N1t = diff(N1,t);
    N2t = diff(N2,t);
    N3t = diff(N3,t);
    N4t = diff(N4,t);
    N5t = diff(N5,t);
    N6t = diff(N6,t);
    N7t = diff(N7,t);
    N8t = diff(N8,t);
    N1u = diff(N1,u);
    N2u = diff(N2,u);
    N3u = diff(N3,u);
    N4u = diff(N4,u);
    N5u = diff(N5,u);
    N6u = diff(N6,u);
    N7u = diff(N7,u);
    N8u = diff(N8,u);
    N1x = N1s*(Yt*Zu - Zt*Yu) - Ys*(N1t*Zu - Zt*N1u) + Zs*(N1t*Yu - Yt*N1u);
    N2x = N2s*(Yt*Zu - Zt*Yu) - Ys*(N2t*Zu - Zt*N2u) + Zs*(N2t*Yu - Yt*N2u);
    N3x = N3s*(Yt*Zu - Zt*Yu) - Ys*(N3t*Zu - Zt*N3u) + Zs*(N3t*Yu - Yt*N3u);
    N4x = N4s*(Yt*Zu - Zt*Yu) - Ys*(N4t*Zu - Zt*N4u) + Zs*(N4t*Yu - Yt*N4u);
    N5x = N5s*(Yt*Zu - Zt*Yu) - Ys*(N5t*Zu - Zt*N5u) + Zs*(N5t*Yu - Yt*N5u);
    N6x = N6s*(Yt*Zu - Zt*Yu) - Ys*(N6t*Zu - Zt*N6u) + Zs*(N6t*Yu - Yt*N6u);
    N7x = N7s*(Yt*Zu - Zt*Yu) - Ys*(N7t*Zu - Zt*N7u) + Zs*(N7t*Yu - Yt*N7u);
    N8x = N8s*(Yt*Zu - Zt*Yu) - Ys*(N8t*Zu - Zt*N8u) + Zs*(N8t*Yu - Yt*N8u);
    N1y = Xs*(N1t*Zu - Zt*N1u) - N1s*(Xt*Zu - Zt*Xu) + Zs*(Xt*N1u - N1t*Xu);
    N2y = Xs*(N2t*Zu - Zt*N2u) - N2s*(Xt*Zu - Zt*Xu) + Zs*(Xt*N2u - N2t*Xu);
    N3y = Xs*(N3t*Zu - Zt*N3u) - N3s*(Xt*Zu - Zt*Xu) + Zs*(Xt*N3u - N3t*Xu);
    N4y = Xs*(N4t*Zu - Zt*N4u) - N4s*(Xt*Zu - Zt*Xu) + Zs*(Xt*N4u - N4t*Xu);
    N5y = Xs*(N5t*Zu - Zt*N5u) - N5s*(Xt*Zu - Zt*Xu) + Zs*(Xt*N5u - N5t*Xu);
    N6y = Xs*(N6t*Zu - Zt*N6u) - N6s*(Xt*Zu - Zt*Xu) + Zs*(Xt*N6u - N6t*Xu);
    N7y = Xs*(N7t*Zu - Zt*N7u) - N7s*(Xt*Zu - Zt*Xu) + Zs*(Xt*N7u - N7t*Xu);
    N8y = Xs*(N8t*Zu - Zt*N8u) - N8s*(Xt*Zu - Zt*Xu) + Zs*(Xt*N8u - N8t*Xu);
    N1z = Xs*(Yt*N1u - N1t*Yu) - Ys*(Xt*N1u - N1t*Xu) + N1s*(Xt*Yu - Yt*Xu);
    N2z = Xs*(Yt*N2u - N2t*Yu) - Ys*(Xt*N2u - N2t*Xu) + N2s*(Xt*Yu - Yt*Xu);
    N3z = Xs*(Yt*N3u - N3t*Yu) - Ys*(Xt*N3u - N3t*Xu) + N3s*(Xt*Yu - Yt*Xu);
    N4z = Xs*(Yt*N4u - N4t*Yu) - Ys*(Xt*N4u - N4t*Xu) + N4s*(Xt*Yu - Yt*Xu);
    N5z = Xs*(Yt*N5u - N5t*Yu) - Ys*(Xt*N5u - N5t*Xu) + N5s*(Xt*Yu - Yt*Xu);
    N6z = Xs*(Yt*N6u - N6t*Yu) - Ys*(Xt*N6u - N6t*Xu) + N6s*(Xt*Yu - Yt*Xu);
    N7z = Xs*(Yt*N7u - N7t*Yu) - Ys*(Xt*N7u - N7t*Xu) + N7s*(Xt*Yu - Yt*Xu);
    N8z = Xs*(Yt*N8u - N8t*Yu) - Ys*(Xt*N8u - N8t*Xu) + N8s*(Xt*Yu - Yt*Xu);
    B = [N1x N2x N3x N4x N5x N6x N7x N8x   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
           0   0   0   0   0   0   0   0 N1y N2y N3y N4y N5y N6y N7y N8y   0   0   0   0   0   0   0   0;
           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 N1z N2z N3z N4z N5z N6z N7z N8z;
         N1y N2y N3y N4y N5y N6y N7y N8y N1x N2x N3x N4x N5x N6x N7x N8x   0   0   0   0   0   0   0   0;
           0   0   0   0   0   0   0   0 N1z N2z N3z N4z N5z N6z N7z N8z N1y N2y N3y N4y N5y N6y N7y N8y;
         N1z N2z N3z N4z N5z N6z N7z N8z   0   0   0   0   0   0   0   0 N1x N2x N3x N4x N5x N6x N7x N8x];
    Bnew = simplify(B);
    Jnew = simplify(J);
    D = (E/((1+v)*(1-2*v)))*[1-v   v   v         0         0         0; 
                               v 1-v   v         0         0         0;
                               v   v 1-v         0         0         0;
                               0   0   0 (1-2*v)/2         0         0;
                               0   0   0         0 (1-2*v)/2         0;
                               0   0   0         0         0 (1-2*v)/2];
    BD = transpose(Bnew)*D*Bnew/Jnew;
    r = int(int(int(BD, u, -1, 1), t, -1, 1), s, -1, 1);
    K = double(r);  
end