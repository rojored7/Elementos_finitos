function Strain = FEM_EP_Strain_xEsf(E,v,Esf,P_Stress_Strain)


%% M�TODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento Plano - Deformaci�n debida al esfuerzo (Strain)

%% Variables

% Entrada:

% E                 ::>  M�dulo de Young (elasticidad) del elemento.
% v                 ::>  Coeficiente de poisson.
% DEs               ::>  Esfuerzo (3x1)
% P_Stress_Strain   ::>  (default = 'P_Stress') 
%                        P_Stress_Strain = 'P_Stress' - Esfuerzo en eje Z = 0
%                        P_Stress_Strain = 'P_Strain' - Deformaci�n en eje Z = 0


% Salida:

% Strain            ::>  Deformaci�n = [Def_x;Def_y;Def_xy]
%                           Def_x  : Deformaci�n normal en x
%                           Def_y  : Deformaci�n normal en y
%                           Def_xy : Deformaci�n cortante en x-y


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 3
    ok = 1;
    P_Stress_Strain = 'P_Stress';
elseif cant_in == 4
    ok = 1;
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la funci�n\n\n')
end

if ok == 1
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
    Strain = D\Esf;
end

