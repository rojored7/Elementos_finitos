function [K_si,T] = FEM_SystemBarInclinedSupport (K,i,tetha,DoR,DOF)

%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Matriz del sistema en 2D y 3D 

%% Variables

% Entrada:

% K        ::>  Matriz de rigidez del elemento.
% i        ::>  Vector de nodos con soportes inclinados (nx1).
% tetha    ::>  Vector de angulos de los soportes inclinados (nx1).
% DoR      ::>  Sistema de coordenasa (ángulos): 'DEG' o 'RAD' (defualt DEG).
% DOF      ::>  Grados de libertad: '2D' o '3D' (defualt 2D).


% Salida:

% K_si     ::>  Matriz de rigidez del 'elemento barra' con transformación para soporte inclinado
% T        ::>  Matriz de transformación para soporte inclinado


%% Programa

cant_in = nargin;
ok = 0;

if cant_in == 3    
    [f_i,c_i]=size(i);
    [f_tetha,c_tetha]=size(tetha);
    if f_i == f_tetha && c_i == c_tetha
        ok = 1;
    else
        fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
    end
    Dim = '2D';
    CoorAng = 'DEG';
elseif cant_in == 4
    [f_i,c_i]=size(i);
    [f_tetha,c_tetha]=size(tetha);
    if f_i == f_tetha && c_i == c_tetha
        ok = 1;
    else
        fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
    end
    Dim = '2D';
    CoorAng = DoR;
elseif cant_in == 5
    [f_i,c_i]=size(i);
    [f_tetha,c_tetha]=size(tetha);
    if f_i == f_tetha && c_i == c_tetha
        ok = 1;
    else
        fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
    end
    Dim = DOF;
    CoorAng = DoR;
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    [f_K,c_K] = size(K); 
    T = eye(f_K);
    for t = 1:length(i)               
        if strcmp(Dim,'2D') == 1
            if strcmp(CoorAng,'DEG') == 1
                x = tetha(t)*pi/180;
            end
            if strcmp(CoorAng,'RAD') == 1
                x = tetha(t);
            end        
            T(2*i(t)-1,2*i(t)-1) =  cos(x);
            T(2*i(t)-1,2*i(t))   =  sin(x);
            T(2*i(t),2*i(t)-1)   = -sin(x) ;
            T(2*i(t),2*i(t))     =  cos(x);
        end
        if strcmp(Dim,'3D') == 1
            if strcmp(CoorAng,'DEG') == 1
                x = tetha(t)*pi/180;
            end
            if strcmp(CoorAng,'RAD') == 1
                x = tetha(t);
            end
            T(3*i(t)-2,3*i(t)-2) =  cos(x);
            T(3*i(t)-2,3*i(t)-1) =  sin(x);
            T(3*i(t)-2,3*i(t))   =  0;
            T(3*i(t)-1,3*i(t)-2) = -sin(x);
            T(3*i(t)-1,3*i(t)-1) =  cos(x);
            T(3*i(t)-1,3*i(t))   =  0;
            T(3*i(t),3*i(t)-2)   =  0;
            T(3*i(t),3*i(t)-1)   =  0;
            T(3*i(t),3*i(t))     =  1;
        end  
    end
    K_si = T*K*T';
end
