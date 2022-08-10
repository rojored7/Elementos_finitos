function [Stress_cent,Stress] = FEM_EP_Quad_Stress(E,v,ui,uj,uk,ul,Des,P_Stress_Strain)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento Plano Rectangular- Esfuerzo (Stress)

%% Variables

% Entrada:

% E                 ::>  Módulo de Young (elasticidad) del elemento.
% v                 ::>  Coeficiente de poisson.
% ui                ::>  Posición nodo i (X,Y)
% uj                ::>  Posición nodo j (X,Y)
% uk                ::>  Posición nodo k (X,Y)
% ul                ::>  Posición nodo l (X,Y)
% Des               ::>  Desplazamientos (8x1)
% P_Stress_Strain   ::>  (default = 'P_Stress') 
%                        P_Stress_Strain = 'P_Stress' - Esfuerzo en eje Z = 0
%                        P_Stress_Strain = 'P_Strain' - Deformación en eje Z = 0


% Salida:

% Stress_cent       ::>  Esfuerzo en el centroide - coordenadas (0,0)
% Stress            ::>  Esfuerzo = [Esf_x;Esf_y;Esf_xy] - en cualquier coordenada (s,h)
%                           Esf_x  : Esfuerzo normal en x
%                           Esf_y  : Esfuerzo normal en y
%                           Esf_xy : Esfuerzo cortante en x-y


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 7
    ok = 1;
    P_Stress_Strain = 'P_Stress';
elseif cant_in == 8
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
    Xl = ul(1);
    Yl = ul(2);
    syms s h
    a = (Yi*(s-1)+Yj*(-1-s)+Yk*(1+s)+Yl*(1-s))/4;
    b = (Yi*(h-1)+Yj*(1-h)+Yk*(1+h)+Yl*(-1-h))/4;
    c = (Xi*(h-1)+Xj*(1-h)+Xk*(1+h)+Xl*(-1-h))/4;
    d = (Xi*(s-1)+Xj*(-1-s)+Xk*(1+s)+Xl*(1-s))/4;
    B1 = [a*(h-1)/4-b*(s-1)/4                   0; 
                            0 c*(s-1)/4-d*(h-1)/4;
          c*(s-1)/4-d*(h-1)/4 a*(h-1)/4-b*(s-1)/4];
    B2 = [a*(1-h)/4-b*(-1-s)/4                    0; 
                             0 c*(-1-s)/4-d*(1-h)/4;
          c*(-1-s)/4-d*(1-h)/4 a*(1-h)/4-b*(-1-s)/4];
    B3 = [a*(h+1)/4-b*(s+1)/4                   0;
                            0 c*(s+1)/4-d*(h+1)/4;
          c*(s+1)/4-d*(h+1)/4 a*(h+1)/4-b*(s+1)/4];
    B4 = [a*(-1-h)/4-b*(1-s)/4                    0;
                             0 c*(1-s)/4-d*(-1-h)/4;
          c*(1-s)/4-d*(-1-h)/4 a*(-1-h)/4-b*(1-s)/4];
    Bfirst = [B1 B2 B3 B4];
    Jfirst = [  0  1-h  h-s  s-1; 
              h-1    0  s+1 -s-h;
              s-h -s-1    0  h+1;
              1-s  s+h -h-1    0];
    J = [Xi Xj Xk Xl]*Jfirst*[Yi;Yj;Yk;Yl]/8;
    B = Bfirst/J;    
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
    Stress = D*B*Des;
    % Esfuerzo en el centroide del elemento
    Stress_cent = subs(Stress,{s,h},{0,0});
    Stress_cent = double(Stress_cent);
end