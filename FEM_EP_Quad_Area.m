function A = FEM_EP_Quad_Area(ui,uj,uk,ul)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento Plano Rectangular - Área

%% Variables

% Entrada:

% ui                ::>  Posición nodo i (X,Y)
% uj                ::>  Posición nodo j (X,Y)
% uk                ::>  Posición nodo k (X,Y)
% ul                ::>  Posición nodo l (X,Y)

% Salida:

% A                 ::>  Área


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 4
    ok = 1;
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    Xi = ui(1);
    Yi = ui(2);
    Xj = uj(1);
    Yj = uj(2);
    Xk = uk(1);
    Yk = uk(2);
    Xl = ul(1);
    Yl = ul(2);
    A = (Xi*(Yj-Yk) + Xj*(Yk-Yi) + Xk*(Yi-Yj))/2+(Xi*(Yk-Yl) + Xk*(Yl-Yi) + Xl*(Yi-Yk))/2;
end

