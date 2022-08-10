function S_VM = FEM_Stress_VM(SPrin)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento Plano o con Volumen - Esfuerzo de Von-Mises

%% Variables

% Entrada:

% SPrin             ::>  Esfuerzos principales [S1,S2,S3]
%                           si S1 o S2 o S3 igual a cero, se calcula
%                           Esfuerzo de Von-Mises para elemento plano

% Salida:

% S_VM              ::>  Esfuerzo de Von-Mises


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 1
    if length(SPrin) == 3 || length(SPrin) == 2
        ok = 1;
    end
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    if length(SPrin) == 2
        SPrin(length(SPrin)+1) = 0;
    end
%     S_VM = (SPrin(1)^2+SPrin(2)^2+SPrin(3)^2-SPrin(1)*SPrin(2)-SPrin(2)*SPrin(3)-SPrin(1)*SPrin(3))^0.5;
    S_VM = (((SPrin(1)-SPrin(2))^2+(SPrin(2)-SPrin(3))^2+(SPrin(1)-SPrin(3))^2)/2)^0.5;
end

