function Strain = FEM_EV_Strain_xEsf(E,v,Esf)


%% M�TODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento con Volumen - Deformaci�n debida al esfuerzo (Strain)

%% Variables

% Entrada:

% E                 ::>  M�dulo de Young (elasticidad) del elemento.
% v                 ::>  Coeficiente de poisson.
% DEs               ::>  Esfuerzo (6x1)


% Salida:

% Strain            ::>  Deformaci�n = [Def_x;Def_y;Def_z;Def_xy;Def_yz;Def_xz]
%                           Def_x  : Deformaci�n normal en x
%                           Def_y  : Deformaci�n normal en y
%                           Def_z  : Deformaci�n normal en z
%                           Def_xy : Deformaci�n cortante en x-y
%                           Def_yz : Deformaci�n cortante en y-z
%                           Def_xz : Deformaci�n cortante en x-z


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 3
    ok = 1;
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la funci�n\n\n')
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

