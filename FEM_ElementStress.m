function Stress = FEM_ElementStress (EoK,LoU,UioA,Uj,DEG_X,DEG_Y_DoR,DEG_Z,DoR)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Estres del elemento

%% Variables

% Entrada:

% EoK       ::>  Módulo de Young (elasticidad) del elemento | K - matriz rigidez elemento.
% LoU       ::>  Longitud del elemento | U - desplazamientos nodos i y j del elemento.
% UioA      ::>  Vector de desplazamiento en el nodo de origen del elemento | A - área del elemento.
% Uj        ::>  Vector de desplazamiento en el nodo de llegada del elemento.
% DEG_X     ::>  Ángulo de dirección del 'elemento barra' según el eje cordenado global (eje X en el espacio).
% DEG_Y_DoR ::>  DEG_Y - Ángulo de dirección del 'elemento barra' según el eje cordenado global en el eje Y.
%                | DoR - (default DEG)'DEG' o 'RAD' (para 2D)
% DEG_Z     ::>  Ángulo de dirección del 'elemento barra' según el eje cordenado global en el eje Z.
% DoR       ::>  (default DEG)'DEG' o 'RAD'

% Salida:

% K_bar2D  ::>  Matriz de rigidez del 'elemento barra'


%% Programa

cant_in = nargin;
ok = 0;

if cant_in >= 3 && cant_in <= 8
    ok = 1;
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    [f_Ui,c_Ui]=size(UioA);
    [f_LoU,c_LoU] = size(LoU);
    [f_EoK,c_EoK] = size(EoK);
    if f_Ui >= 1 && c_Ui == 1
        if f_LoU == 2 && f_EoK == 2
            Dim = '1D_caso2';
        end
        if f_Ui == 1 && (cant_in == 4)
            Dim = '1D';
        elseif f_Ui == 2 && (cant_in == 5 || cant_in == 6)
            Dim = '2D';
        elseif f_Ui == 3 && (cant_in == 7 || cant_in == 8)
            Dim = '3D';
        end
        if strcmp(Dim,'1D_caso2') == 1
            Stress = EoK*LoU/UioA;
        end
        if strcmp(Dim,'1D') == 1
            Stress = EoK/LoU*[-1 1]*[UioA;Uj];
        end
        if strcmp(Dim,'2D') == 1
            if cant_in == 5 || strcmp(DEG_Y_DoR,'DEG') == 1
                C = cosd(DEG_X);
                S = sind(DEG_X);
            elseif strcmp(DEG_Y_DoR,'RAD') == 1
                C = cos(DEG_X);
                S = sin(DEG_X);
            end
            Stress = EoK/LoU*[-C -S C S]*[UioA;Uj];
        end
        if strcmp(Dim,'3D') == 1
            if cant_in == 7 || strcmp(DoR,'DEG') == 1
                Cx = cosd(DEG_X);
                Cy = cosd(DEG_Y_DoR);
                Cz = cosd(DEG_Z);
            elseif strcmp(DoR,'RAD') == 1
                Cx = cos(DEG_X);
                Cy = cos(DEG_Y_DoR);
                Cz = cos(DEG_Z);
            end
            Stress = EoK/LoU*[-Cx -Cy -Cz Cx Cy Cz]*[UioA;Uj];
        end
    else
        fprintf('\n\nVerifique las dimensiones de los datos de entrada en la función\n\n')    
    end
end