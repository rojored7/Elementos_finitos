function [Cx,Cy,Cz,L] = FEM_cxcycz (Vec_O,Vec_E)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Calculo de los cosenos en X, Y y Z entre 2 puntos en el espacio

%% Variables

% Entrada:

% Vec_O    ::>  Vector de posiciones X, Y y Z en el espacio del nodo de origen.
% Vec_E    ::>  Vector de posiciones X, Y y Z en el espacio del nodo de llegada.

% Salida:

% Cx       ::>  Matriz de rigidez del 'elemento barra'
% Cy       ::>  Matriz de rigidez del 'elemento barra'
% Cz       ::>  Matriz de rigidez del 'elemento barra'
% L        ::>  Longitud del elemento.


%% Programa


cant_in = nargin;
ok = 0;

if cant_in == 2
    ok = 1;
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    L = ((Vec_E(1)-Vec_O(1))^2+(Vec_E(2)-Vec_O(2))^2+(Vec_E(3)-Vec_O(3))^2)^(1/2);
    Cx = (Vec_E(1)-Vec_O(1))/L;
    Cy = (Vec_E(2)-Vec_O(2))/L;
    Cz = (Vec_E(3)-Vec_O(3))/L;
end

