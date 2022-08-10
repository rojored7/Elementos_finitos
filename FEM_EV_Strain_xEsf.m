function Strain = FEM_EV_Strain_xEsf(E,v,Esf)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento con Volumen - Deformación debida al esfuerzo (Strain)

%% Variables

% Entrada:

% E                 ::>  Módulo de Young (elasticidad) del elemento.
% v                 ::>  Coeficiente de poisson.
% DEs               ::>  Esfuerzo (6x1)


% Salida:

% Strain            ::>  Deformación = [Def_x;Def_y;Def_z;Def_xy;Def_yz;Def_xz]
%                           Def_x  : Deformación normal en x
%                           Def_y  : Deformación normal en y
%                           Def_z  : Deformación normal en z
%                           Def_xy : Deformación cortante en x-y
%                           Def_yz : Deformación cortante en y-z
%                           Def_xz : Deformación cortante en x-z


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 3
    ok = 1;
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    D = (E/((1+v)*(1-2*v)))*[1-v   v   v         0         0         0; 
                               v 1-v   v         0         0         0;
                               v   v 1-v         0         0         0;
                               0   0   0 (1-2*v)/2         0         0;
                               0   0   0         0 (1-2*v)/2         0;
                               0   0   0         0         0 (1-2*v)/2];
    Strain = D\Esf;
end

