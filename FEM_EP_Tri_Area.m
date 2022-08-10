function A = FEM_EP_Tri_Area(ui,uj,uk)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento Plano Triangular - Área

%% Variables

% Entrada:

% ui                ::>  Posición nodo i (X,Y)
% uj                ::>  Posición nodo j (X,Y)
% uk                ::>  Posición nodo k (X,Y)

% Salida:

% A                 ::>  Área


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 3
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
    A = (Xi*(Yj-Yk) + Xj*(Yk-Yi) + Xk*(Yi-Yj))/2;
end

