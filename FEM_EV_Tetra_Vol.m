function V = FEM_EV_Tetra_Vol(ui,uj,uk,ul)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento con Volumen (Tetraedro) - Volumen

%% Variables

% Entrada:

% ui                ::>  Posición nodo i (X,Y,Z)
% uj                ::>  Posición nodo j (X,Y,Z)
% uk                ::>  Posición nodo k (X,Y,Z)
% ul                ::>  Posición nodo l (X,Y,Z)

% Salida:

% V                 ::>  Volumen


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 4
    if length(ui) == 3 && length(uj) == 3 && length(uk) == 3 && length(ul) == 3
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
    V = det([1 Xi Yi Zi; 
             1 Xj Yj Zj; 
             1 Xk Yk Zk; 
             1 Xl Yl Zl])/6;
end

