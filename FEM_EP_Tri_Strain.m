function Strain = FEM_EP_Tri_Strain(ui,uj,uk,Des)


%% M�TODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento Plano Triangular - Deformaci�n (Strain)

%% Variables

% Entrada:

% ui                ::>  Posici�n nodo i (X,Y)
% uj                ::>  Posici�n nodo j (X,Y)
% uk                ::>  Posici�n nodo k (X,Y)
% Des               ::>  Desplazamientos (6x1)


% Salida:

% Strain            ::>  Deformaci�n = [Def_x;Def_y;Def_xy]
%                           Def_x  : Deformaci�n normal en x
%                           Def_y  : Deformaci�n normal en y
%                           Def_xy : Deformaci�n cortante en x-y


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 4
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
    beta_i = Yj-Yk;
    beta_j = Yk-Yi;
    beta_k = Yi-Yj;
    gamma_i = Xk-Xj;
    gamma_j = Xi-Xk;
    gamma_k = Xj-Xi;
    B = [ beta_i       0  beta_j       0  beta_k       0; 
               0 gamma_i       0 gamma_j       0 gamma_k;
         gamma_i  beta_i gamma_j  beta_j gamma_k  beta_k]/(2*A);
    Strain = B*Des;
end

