function Stress = FEM_EP_Stress_xDef(E,v,Def,P_Stress_Strain)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento Plano - Esfuerzo debido a la deformación (Stress)

%% Variables

% Entrada:

% E                 ::>  Módulo de Young (elasticidad) del elemento.
% v                 ::>  Coeficiente de poisson.
% Def               ::>  Defrmación (3x1)
% P_Stress_Strain   ::>  (default = 'P_Stress') 
%                        P_Stress_Strain = 'P_Stress' - Esfuerzo en eje Z = 0
%                        P_Stress_Strain = 'P_Strain' - Deformación en eje Z = 0


% Salida:

% Stress            ::>  Esfuerzo = [Esf_x;Esf_y;Esf_xy]
%                           Esf_x  : Esfuerzo normal en x
%                           Esf_y  : Esfuerzo normal en y
%                           Esf_xy : Esfuerzo cortante en x-y


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 3
    ok = 1;
    P_Stress_Strain = 'P_Stress';
elseif cant_in == 4
    ok = 1;
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
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
    Stress = D*Def;
end

