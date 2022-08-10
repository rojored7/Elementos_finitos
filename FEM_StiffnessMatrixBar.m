function K_bar = FEM_StiffnessMatrixBar (E,A,L,DEG_X,DEG_YoDoR,DEG_Z,DoR)

%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Matriz de rigidez para barras

%% Variables

% Entrada:

% E        ::>  Módulo de Young (elasticidad) del elemento.
% A        ::>  Área transversal del elemento.
% L        ::>  Longitud del elemento.
% DEG_X    ::>  Ángulo de dirección del 'elemento barra' según el eje cordenado global (eje X en el espacio).
% DEG_Y    ::>  Ángulo de dirección del 'elemento barra' según el eje cordenado global en el eje Y.
% DEG_Z    ::>  Ángulo de dirección del 'elemento barra' según el eje cordenado global en el eje Z.
% DoR      ::>  'DEG' o 'RAD'

%   DATOS DE ENTRADA:
%       Para Barras 1D las entradas son: E,A,L
%       Para Barras 2D las entradas son: E,A,L,DEG_X
%       Para Barras 3D las entradas son: E,A,L,DEG_X,DEG_Y,DEG_Z

% Salida:

% K_bar2D  ::>  Matriz de rigidez del 'elemento barra'


%% Programa


cant_in = nargin;
ok = 0;

if cant_in == 3
    ok = 1;
    Dim = '1D';  
elseif cant_in == 4
    ok = 1;
    Dim = '2D';    
    DEG_YoDoR = 'DEG';
elseif cant_in == 5
    ok = 1;
    Dim = '2D'; 
elseif cant_in == 6
    ok = 1;
    Dim = '3D';
    DoR = 'DEG';
elseif cant_in == 7
    ok = 1;
    Dim = '3D';
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    if strcmp(Dim,'1D') == 1
        K_bar = E*A/L*[ 1 -1;
                       -1  1];
    end
    if strcmp(Dim,'2D') == 1
        if strcmp(DEG_YoDoR,'DEG') == 1
            C = cosd(DEG_X);
            S = sind(DEG_X);
        elseif strcmp(DEG_YoDoR,'RAD') == 1
            C = cos(DEG_X);
            S = sin(DEG_X);
        end
        C2 = C^2;
        S2 = S^2;
        CS = C*S;
        Lamda = [C2 CS;
                 CS S2];
        K_bar = E*A/L*[ Lamda -Lamda;
                       -Lamda  Lamda;];
    end
    if strcmp(Dim,'3D') == 1
        % Tabla de sin-cos de los ángulos
        if strcmp(DoR,'DEG') == 1
            Cx = cosd(DEG_X);
            Cy = cosd(DEG_YoDoR);
            Cz = cosd(DEG_Z);
        elseif strcmp(DoR,'RAD') == 1
            Cx = cos(DEG_X);
            Cy = cos(DEG_YoDoR);
            Cz = cos(DEG_Z);
        end        
        Cx2 = Cx^2;
        Cy2 = Cy^2;
        Cz2 = Cz^2;
        Cxy = Cx*Cy;
        Cxz = Cx*Cz;
        Cyz = Cy*Cz;
        Lamda = [Cx2 Cxy Cxz;
                 Cxy Cy2 Cyz;
                 Cxz Cyz Cz2];
        K_bar = E*A/L*[ Lamda -Lamda;
                       -Lamda  Lamda;];
    end
end