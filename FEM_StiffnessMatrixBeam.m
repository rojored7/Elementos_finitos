function K_beam = FEM_StiffnessMatrixBeam (E,I,L_A,A_G,DEG_J,DEG_COOR_DoR_ui,DoR_uj)

%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Matriz de rigidez para vigas

%% Variables

% Entrada:

% E            ::>  Módulo de Young (elasticidad) del elemento.
% I            ::>  Momento de inercia del elemento.
%                   | I = [Iy,Iz] cuando se trabaja con 3D
% L            ::>  L - Longitud del elemento.
%                   | A - área (para 3d_3a)
% AoG          ::>  A - área | G - módulo de torsión del material.
% DEGoJ        ::>  DEG - Ángulo de dirección del 'elemento barra' según el eje cordenado global
%                   | J - constante torsional del material.
% DEGoCOOR_DoR ::>  DEG - Ángulo de dirección del 'elemento barra' según el eje cordenado global
%                   | COOR - matriz 2x2 con coordenadas del elemento [Xj,Xi;Zj,Zi].
%                   | DoR - (default DEG)'DEG' o 'RAD' (para 2D_1A)
%                   | ui - Coordenadas X,Y,Z del nodo de origen i (para 3d_3a)
% DoR          ::>  DoR - (default DEG)'DEG' o 'RAD'
%                   | uj - Coordenadas X,Y,Z del nodo de llegada j (para 3d_3a)

%   DATOS DE ENTRADA:
%       Para Vigas con 1d y 1a las entradas son: E,I,L
%       Para Vigas con 2d y 1a las entradas son: E,I,L,A,DEG
%       Para Vigas con 1d y 2a las entradas son: E,I,L,G,J,DEG
%       Para Vigas con 3d y 3a las entradas son: E,[Iy,Iz],A,G,J,ui(x,y,x),uj(x,y,z)


% Salida:

% K_bar2D      ::>  Matriz de rigidez del 'elemento barra'


% Nota: d - Desplazamiento a los largo de los ejes X, Y o Z
%       a - Ángulo en el eje coordenado cartesiano



%% Programa


cant_in = nargin;
ok = 0;

if cant_in == 3
    ok = 1;
    Dim = '1d_1a';
elseif cant_in == 5
    ok = 1;
    Dim = '2d_1a';
    DEG_COOR_DoR_ui = 'DEG';
elseif cant_in == 6
    ok = 1;
    Dim = '1d_2a';
    DoR_uj = 'DEG';
    if strcmp(DEG_COOR_DoR_ui,'DEG') == 1 || strcmp(DEG_COOR_DoR_ui,'RAD') == 1
        Dim = '2D_1A';
    end
elseif cant_in == 7
    ok = 1;
    if length(I) == 2
        Dim = '3d_3a';
    else
    	Dim = '1d_2a';
    end
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    if strcmp(Dim,'1d_1a') == 1
        L = L_A;
        K_beam = (E*I/L^3)*[ 12   6*L  -12   6*L;
                            6*L 4*L^2 -6*L 2*L^2;
                            -12  -6*L   12  -6*L;
                            6*L 2*L^2 -6*L 4*L^2];
    end
    if strcmp(Dim,'2d_1a') == 1
        I = I(1);
        L = L_A;
        A = A_G;
        DEG = DEG_J;
        if strcmp(DEG_COOR_DoR_ui,'DEG') == 1
            C = cosd(DEG);
            S = sind(DEG);
        elseif strcmp(DEG_COOR_DoR_ui,'RAD') == 1
            C = cos(DEG);
            S = sin(DEG);
        end
        ka1 = A*C^2+12*I*S^2/(L^2);
        ka2 = A*S^2+12*I*C^2/(L^2);
        ka3 = (A-12*I/(L^2))*C*S;
        ka4 = 6*I*S/L;
        ka5 = 6*I*C/L;
        ka6 = 4*I;
        ka7 = 2*I;
        K_s1 = [ ka1  ka3 -ka4;
                 ka3  ka2  ka5;
                -ka4  ka5  ka6];
        K_s2 = [-ka1 -ka3 -ka4;
                -ka3 -ka2  ka5;
                 ka4 -ka5  ka7];
        K_s3 = [ ka1  ka3  ka4;
                 ka3  ka2 -ka5;
                 ka4 -ka5  ka6];
        K_beam = E/L*[K_s1  K_s2;
                      K_s2' K_s3];
    end
    if strcmp(Dim,'1d_2a') == 1     
        [f_DEGoCOOR,c_DEGoCOOR] = size(DEG_COOR_DoR_ui);
        if f_DEGoCOOR == c_DEGoCOOR && f_DEGoCOOR == 2
            COOR = DEG_COOR_DoR_ui;
            L = sqrt((COOR(1,1)-COOR(1,2))^2+(COOR(2,1)-COOR(2,2))^2);
            C = (COOR(1,1)-COOR(1,2))/L;
            S = (COOR(2,1)-COOR(2,2))/L;
        elseif f_DEGoCOOR == c_DEGoCOOR && f_DEGoCOOR == 1
            DEG = DEG_COOR_DoR_ui;
            if strcmp(DoR_uj,'DEG') == 1
                C = cosd(DEG);
                S = sind(DEG);
            end
            if strcmp(DoR_uj,'RAD') == 1
                C = cos(DEG);
                S = sin(DEG);
            end
        end
        L = L_A;
        I = I(1);
        G = A_G;
        J = DEG_J;
        Tg = [ 1  0  0  0  0  0;
               0  C  S  0  0  0;
               0 -S  C  0  0  0;
               0  0  0  1  0  0;
               0  0  0  0  C  S;
               0  0  0  0 -S  C];
        Kg = [ 12*E*I/L^3      0  6*E*I/L^2 -12*E*I/L^3      0  6*E*I/L^2;
                        0  G*J/L          0           0 -G*J/L          0;
                6*E*I/L^2      0    4*E*I/L  -6*E*I/L^2      0    2*E*I/L;
              -12*E*I/L^3      0 -6*E*I/L^2  12*E*I/L^3      0 -6*E*I/L^2;
                        0 -G*J/L          0           0  G*J/L          0;
                6*E*I/L^2      0    2*E*I/L  -6*E*I/L^2      0    4*E*I/L];
        K_beam = Tg'*Kg*Tg;
    end
    if strcmp(Dim,'3d_3a') == 1        
        ui = DEG_COOR_DoR_ui;
        uj = DoR_uj;
        Xi = ui(1);
        Yi = ui(2);
        Zi = ui(3);
        Xj = uj(1);
        Yj = uj(2);        
        Zj = uj(3);
        Iy = I(1);
        Iz = I(2);
        A = L_A;
        G = A_G;
        J = DEG_J;
        L = sqrt((Xj-Xi)*(Xj-Xi) + (Yj-Yi)*(Yj-Yi) + (Zj-Zi)*(Zj-Zi));
        w1 = E*A/L;
        w2 = 12*E*Iz/(L*L*L);
        w3 = 6*E*Iz/(L*L);
        w4 = 4*E*Iz/L;
        w5 = 2*E*Iz/L;
        w6 = 12*E*Iy/(L*L*L);
        w7 = 6*E*Iy/(L*L);
        w8 = 4*E*Iy/L;
        w9 = 2*E*Iy/L;
        w10 = G*J/L;
        ke = [ w1   0   0    0   0   0 -w1   0   0    0   0   0;
                0  w2   0    0   0  w3   0 -w2   0    0   0  w3;
                0   0  w6    0 -w7   0   0   0 -w6    0 -w7   0;
                0   0   0  w10   0   0   0   0   0 -w10   0   0;
                0   0 -w7    0  w8   0   0   0  w7    0  w9   0;
                0  w3   0    0   0  w4   0 -w3   0    0   0  w5;
              -w1   0   0    0   0   0  w1   0   0    0   0   0;
                0 -w2   0    0   0 -w3   0  w2   0    0   0 -w3;
                0   0 -w6    0  w7   0   0   0  w6    0  w7   0;
                0   0   0 -w10   0   0   0   0   0  w10   0   0;
                0   0 -w7    0  w9   0   0   0  w7    0  w8   0;
                0  w3   0    0   0  w5   0 -w3   0    0   0  w4];
        if Xi == Xj && Yi == Yj
           if Zj > Zi
              Lambda = [ 0 0 1; 
                         0 1 0;
                        -1 0 0];
           else
              Lambda = [0 0 -1; 
                        0 1  0; 
                        1 0  0];
           end
        else
            CXx = (Xj-Xi)/L;
            CYx = (Yj-Yi)/L;
            CZx = (Zj-Zi)/L;
            D = sqrt(CXx*CXx + CYx*CYx);
            CXy = -CYx/D;
            CYy = CXx/D;
            CZy = 0;
            CXz = -CXx*CZx/D;
            CYz = -CYx*CZx/D;
            CZz = D;
            Lambda = [CXx CYx CZx; 
                      CXy CYy CZy;
                      CXz CYz CZz];
        end
        R = [  Lambda zeros(3) zeros(3) zeros(3); 
             zeros(3)   Lambda zeros(3) zeros(3);
             zeros(3) zeros(3)   Lambda zeros(3);
             zeros(3) zeros(3) zeros(3)   Lambda];
        K_beam = R'*ke*R;
    end
end