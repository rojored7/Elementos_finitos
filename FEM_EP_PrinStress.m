function [S1,S2,t] = FEM_EP_PrinStress(Stress)


%% M�TODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento Plano - Esfuerzo Principal (Stress)

%% Variables

% Entrada:

% Stress       ::> Estr�s.

% Salida:

% S1           ::>  Esfuerzo principal m�ximo
% S2           ::>  Esfuerzo principal minimo
% t            ::>  Angulo principal (inclinacion) [DEG]


%% Programa


cant_in = nargin;
ok = 0;

if cant_in == 1
    ok = 1;
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la funci�n\n\n')
end

if ok == 1
    R = (Stress(1)+Stress(2))/2;
    Q = ((Stress(1)-Stress(2))/2)^2+Stress(3)^2;
    M = 2*Stress(3)/(Stress(1)-Stress(2));
    S1 = R+sqrt(Q);
    S2 = R-sqrt(Q);
    t = (atan(M)/2)*180/pi;
end

