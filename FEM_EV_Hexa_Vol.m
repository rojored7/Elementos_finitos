function V = FEM_EV_Hexa_Vol(ui,uj,uk,ul,um,un,uo,up)

%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento con Volumen (Hexaedro) - Volumen

%% Variables

% Entrada:

% ui                ::>  Posición nodo i (X,Y,Z)
% uj                ::>  Posición nodo j (X,Y,Z)
% uk                ::>  Posición nodo k (X,Y,Z)
% ul                ::>  Posición nodo l (X,Y,Z)
% um                ::>  Posición nodo m (X,Y,Z)
% un                ::>  Posición nodo n (X,Y,Z)
% uo                ::>  Posición nodo o (X,Y,Z)
% up                ::>  Posición nodo p (X,Y,Z)

% Salida:

% V                 ::>  Volumen


%% Programa

cant_in = nargin;
ok = 0;

if cant_in == 8
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
    Jnew = simplify(J);
    r = int(int(int(Jnew, u, -1, 1), t, -1, 1), s, -1, 1);
    V = double(r);
end