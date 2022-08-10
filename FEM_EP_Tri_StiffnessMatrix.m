function K = FEM_EP_Tri_StiffnessMatrix(E,v,t,ui,uj,uk,P_Stress_Strain)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento Plano Triangular - Matriz de Rigidez

%% Variables

% Entrada:

% E                 ::>  Módulo de Young (elasticidad) del elemento.
% v                 ::>  Coeficiente de poisson.
% t                 ::>  Espesor.
% ui                ::>  Posición nodo i (X,Y)
% uj                ::>  Posición nodo j (X,Y)
% uk                ::>  Posición nodo k (X,Y)
% P_Stress_Strain   ::>  (default = 'P_Stress') 
%                        P_Stress_Strain = 'P_Stress' - Esfuerzo en eje Z = 0
%                        P_Stress_Strain = 'P_Strain' - Deformación en eje Z = 0

% Salida:

% K                 ::>  Matriz de Rigidez


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 6
    ok = 1;
    P_Stress_Strain = 'P_Stress';
elseif cant_in == 7
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
    beta_i = Yj-Yk;
    beta_j = Yk-Yi;
    beta_k = Yi-Yj;
    gamma_i = Xk-Xj;
    gamma_j = Xi-Xk;
    gamma_k = Xj-Xi;
    B = [ beta_i       0  beta_j       0  beta_k       0; 
               0 gamma_i       0 gamma_j       0 gamma_k;
         gamma_i  beta_i gamma_j  beta_j gamma_k  beta_k]/(2*A);
    if strcmp(P_Stress_Strain,'P_Stress') == 1
        D = (E/(1-v^2))*[1 v       0; 
                         v 1       0; 
                         0 0 (1-v)/2];
    end
    if strcmp(P_Stress_Strain,'P_Strain') == 1
        D = (E/((1+v)*(1-2*v)))*[1-v   v         0; 
                                 v   1-v         0; 
                                 0     0 (1-2*v)/2];
    end
    K = t*A*B'*D*B;
end

