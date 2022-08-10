function [S1,S2,S3] = FEM_EV_PrinStress(Stress)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento con Volumen - Esfuerzos Principales (Stress)

%% Variables

% Entrada:

% Stress       ::> Estrés.

% Salida:

% S1           ::>  Esfuerzo principal X
% S2           ::>  Esfuerzo principal Y
% S3           ::>  Esfuerzo principal Z


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 1
    if length(Stress) == 6
        ok = 1;
    end
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    S1 = Stress(1) + Stress(2) + Stress(3);
    S2 = Stress(1)*Stress(2) + Stress(1)*Stress(3) + Stress(2)*Stress(3) - Stress(4)^2 -Stress(5)^2 -Stress(6)^2;
    ms3 = [Stress(1) Stress(4) Stress(6); 
           Stress(4) Stress(2) Stress(5); 
           Stress(6) Stress(5) Stress(3)];
    S3 = det(ms3);
end

