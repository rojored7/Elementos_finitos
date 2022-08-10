function A = FEM_EP_Tri_Area(ui,uj,uk)


%% M�TODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento Plano Triangular - �rea

%% Variables

% Entrada:

% ui                ::>  Posici�n nodo i (X,Y)
% uj                ::>  Posici�n nodo j (X,Y)
% uk                ::>  Posici�n nodo k (X,Y)

% Salida:

% A                 ::>  �rea


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 3
    ok = 1;
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la funci�n\n\n')
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

