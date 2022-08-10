function K_global = FEM_AssembleElementMatrix (K,ke,i,j,k,l,m,n,o,p)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Ensamblar Matriz de rigidez de los elementos

%% Variables

% Entrada:

% K        ::>  Matriz de rigidez global.
% ke       ::>  Matriz de rigidez a ensamblar.
% nodos    ::>  Para elementos 1D = i - j
%          ::>  Para elementos 2D (Elemento plano Triangular)= i - j - k
%          ::>  Para elementos 2D (Elemento plano Rectangular)= i - j - k - l
%          ::>  Para elementos 3D (Elemento Tetraedro)= i - j - k - l
%          ::>  Para elementos 3D (Elemento Hexaedro)= i - j - k - l - m - n - o - p

% Salida:

% K_global ::>  Matriz de rigidez del elemento ensamblada


%% Programa

cant_in = nargin;
ok = 0;

if cant_in == 4
    ok = 1;
    [f_k,c_k]=size(ke);
    if f_k == 2
        Dim = '1D_1DOF';
    elseif f_k == 4
        Dim = '1D_2DOF';
    elseif f_k == 6
        Dim = '1D_3DOF';
    elseif f_k == 12
        Dim = '1D_6DOF';
    end
elseif cant_in == 5
    ok = 1;
    Dim = '2D_eTrian';
elseif cant_in == 6
    ok = 1;
    [f_k,c_k]=size(ke);
    if f_k == 8
        Dim = '2D_eQuad';
    elseif f_k == 12
        Dim = '3D_Tetra';
    end
elseif cant_in == 10
    ok = 1;
    Dim = '3D_Hexa';
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    [f_k,c_k]=size(ke);
    if f_k == c_k
        if strcmp(Dim,'1D_1DOF') == 1
            pos = [i,j];
            for z = 1:1:length(pos)
                K(pos(z),i) = K(pos(z),i) + ke(z,1);
                K(pos(z),j) = K(pos(z),j) + ke(z,2);
            end
        end
        if strcmp(Dim,'1D_2DOF') == 1
            pos = [i,j];
            dof = 2;
            for z = 1:1:length(pos)
                for d = dof-1:-1:0
                    K(2*pos(z)-d,2*i-1) = K(2*pos(z)-d,2*i-1) + ke(2*z-d,1);
                    K(2*pos(z)-d,2*i)   = K(2*pos(z)-d,2*i)   + ke(2*z-d,2);
                    K(2*pos(z)-d,2*j-1) = K(2*pos(z)-d,2*j-1) + ke(2*z-d,3);
                    K(2*pos(z)-d,2*j)   = K(2*pos(z)-d,2*j)   + ke(2*z-d,4);
                end
            end
        end
        if strcmp(Dim,'1D_3DOF') == 1
            pos = [i,j];
            dof = 3;
            for z = 1:1:length(pos)
                for d = dof-1:-1:0
                    K(3*pos(z)-d,3*i-2) = K(3*pos(z)-d,3*i-2) + ke(3*z-d,1);
                    K(3*pos(z)-d,3*i-1) = K(3*pos(z)-d,3*i-1) + ke(3*z-d,2);
                    K(3*pos(z)-d,3*i)   = K(3*pos(z)-d,3*i)   + ke(3*z-d,3);
                    K(3*pos(z)-d,3*j-2) = K(3*pos(z)-d,3*j-2) + ke(3*z-d,4);
                    K(3*pos(z)-d,3*j-1) = K(3*pos(z)-d,3*j-1) + ke(3*z-d,5);
                    K(3*pos(z)-d,3*j)   = K(3*pos(z)-d,3*j)   + ke(3*z-d,6);
                end
            end
        end
        if strcmp(Dim,'1D_6DOF') == 1
            pos = [i,j];
            dof = 6;
            for z = 1:1:length(pos)                
                for d = dof-1:-1:0
                    K(6*pos(z)-d,6*i-5) = K(6*pos(z)-d,6*i-5) + ke(6*z-d,1);
                    K(6*pos(z)-d,6*i-4) = K(6*pos(z)-d,6*i-4) + ke(6*z-d,2);
                    K(6*pos(z)-d,6*i-3) = K(6*pos(z)-d,6*i-3) + ke(6*z-d,3);
                    K(6*pos(z)-d,6*i-2) = K(6*pos(z)-d,6*i-2) + ke(6*z-d,4);
                    K(6*pos(z)-d,6*i-1) = K(6*pos(z)-d,6*i-1) + ke(6*z-d,5);
                    K(6*pos(z)-d,6*i)   = K(6*pos(z)-d,6*i)   + ke(6*z-d,6);
                    K(6*pos(z)-d,6*j-5) = K(6*pos(z)-d,6*j-5) + ke(6*z-d,7);
                    K(6*pos(z)-d,6*j-4) = K(6*pos(z)-d,6*j-4) + ke(6*z-d,8);
                    K(6*pos(z)-d,6*j-3) = K(6*pos(z)-d,6*j-3) + ke(6*z-d,9);
                    K(6*pos(z)-d,6*j-2) = K(6*pos(z)-d,6*j-2) + ke(6*z-d,10);
                    K(6*pos(z)-d,6*j-1) = K(6*pos(z)-d,6*j-1) + ke(6*z-d,11);
                    K(6*pos(z)-d,6*j)   = K(6*pos(z)-d,6*j)   + ke(6*z-d,12);
                end
            end
        end
        if strcmp(Dim,'2D_eTrian') == 1
            pos = [i,j,k];
            xyz = 2;
            for z = 1:1:length(pos)                
                for eje = xyz-1:-1:0
                    K(2*pos(z)-eje,2*i-1) = K(2*pos(z)-eje,2*i-1) + ke(2*z-eje,1);
                    K(2*pos(z)-eje,2*i)   = K(2*pos(z)-eje,2*i)   + ke(2*z-eje,2);
                    K(2*pos(z)-eje,2*j-1) = K(2*pos(z)-eje,2*j-1) + ke(2*z-eje,3);
                    K(2*pos(z)-eje,2*j)   = K(2*pos(z)-eje,2*j)   + ke(2*z-eje,4);
                    K(2*pos(z)-eje,2*k-1) = K(2*pos(z)-eje,2*k-1) + ke(2*z-eje,5);
                    K(2*pos(z)-eje,2*k)   = K(2*pos(z)-eje,2*k)   + ke(2*z-eje,6);
                end
            end
        end
        if strcmp(Dim,'2D_eQuad') == 1
            pos = [i,j,k,l];
            xyz = 2;
            for z = 1:1:length(pos)
                for eje = xyz-1:-1:0
                    K(2*pos(z)-eje,2*i-1) = K(2*pos(z)-eje,2*i-1) + ke(2*z-eje,1);
                    K(2*pos(z)-eje,2*i)   = K(2*pos(z)-eje,2*i)   + ke(2*z-eje,2);
                    K(2*pos(z)-eje,2*j-1) = K(2*pos(z)-eje,2*j-1) + ke(2*z-eje,3);
                    K(2*pos(z)-eje,2*j)   = K(2*pos(z)-eje,2*j)   + ke(2*z-eje,4);
                    K(2*pos(z)-eje,2*k-1) = K(2*pos(z)-eje,2*k-1) + ke(2*z-eje,5);
                    K(2*pos(z)-eje,2*k)   = K(2*pos(z)-eje,2*k)   + ke(2*z-eje,6);
                    K(2*pos(z)-eje,2*l-1) = K(2*pos(z)-eje,2*l-1) + ke(2*z-eje,7);
                    K(2*pos(z)-eje,2*l)   = K(2*pos(z)-eje,2*l)   + ke(2*z-eje,8);
                end
            end
        end
        if strcmp(Dim,'3D_Tetra') == 1
            pos = [i,j,k,l];
            xyz = 3;
            for z = 1:1:length(pos)
                for eje = xyz-1:-1:0
                    K(3*pos(z)-eje,3*i-2) = K(3*pos(z)-eje,3*i-2) + ke(3*z-eje,1);
                    K(3*pos(z)-eje,3*i-1) = K(3*pos(z)-eje,3*i-1) + ke(3*z-eje,2);
                    K(3*pos(z)-eje,3*i)   = K(3*pos(z)-eje,3*i)   + ke(3*z-eje,3);
                    K(3*pos(z)-eje,3*j-2) = K(3*pos(z)-eje,3*j-2) + ke(3*z-eje,4);
                    K(3*pos(z)-eje,3*j-1) = K(3*pos(z)-eje,3*j-1) + ke(3*z-eje,5);
                    K(3*pos(z)-eje,3*j)   = K(3*pos(z)-eje,3*j)   + ke(3*z-eje,6);                    
                    K(3*pos(z)-eje,3*k-2) = K(3*pos(z)-eje,3*k-2) + ke(3*z-eje,7);
                    K(3*pos(z)-eje,3*k-1) = K(3*pos(z)-eje,3*k-1) + ke(3*z-eje,8);
                    K(3*pos(z)-eje,3*k)   = K(3*pos(z)-eje,3*k)   + ke(3*z-eje,9);
                    K(3*pos(z)-eje,3*l-2) = K(3*pos(z)-eje,3*l-2) + ke(3*z-eje,10);
                    K(3*pos(z)-eje,3*l-1) = K(3*pos(z)-eje,3*l-1) + ke(3*z-eje,11);
                    K(3*pos(z)-eje,3*l)   = K(3*pos(z)-eje,3*l)   + ke(3*z-eje,12);
                end
            end
        end
        if strcmp(Dim,'3D_Hexa') == 1
            pos = [i,j,k,l,m,n,o,p];
            xyz = 3;
            for z = 1:1:length(pos)
                for eje = xyz-1:-1:0
                    K(3*pos(z)-eje,3*i-2) = K(3*pos(z)-eje,3*i-2) + ke(3*z-eje,1);
                    K(3*pos(z)-eje,3*i-1) = K(3*pos(z)-eje,3*i-1) + ke(3*z-eje,2);
                    K(3*pos(z)-eje,3*i)   = K(3*pos(z)-eje,3*i)   + ke(3*z-eje,3);
                    K(3*pos(z)-eje,3*j-2) = K(3*pos(z)-eje,3*j-2) + ke(3*z-eje,4);
                    K(3*pos(z)-eje,3*j-1) = K(3*pos(z)-eje,3*j-1) + ke(3*z-eje,5);
                    K(3*pos(z)-eje,3*j)   = K(3*pos(z)-eje,3*j)   + ke(3*z-eje,6);
                    K(3*pos(z)-eje,3*k-2) = K(3*pos(z)-eje,3*k-2) + ke(3*z-eje,7);
                    K(3*pos(z)-eje,3*k-1) = K(3*pos(z)-eje,3*k-1) + ke(3*z-eje,8);
                    K(3*pos(z)-eje,3*k)   = K(3*pos(z)-eje,3*k)   + ke(3*z-eje,9);
                    K(3*pos(z)-eje,3*l-2) = K(3*pos(z)-eje,3*l-2) + ke(3*z-eje,10);
                    K(3*pos(z)-eje,3*l-1) = K(3*pos(z)-eje,3*l-1) + ke(3*z-eje,11);
                    K(3*pos(z)-eje,3*l)   = K(3*pos(z)-eje,3*l)   + ke(3*z-eje,12);                    
                    K(3*pos(z)-eje,3*m-2) = K(3*pos(z)-eje,3*m-2) + ke(3*z-eje,13);
                    K(3*pos(z)-eje,3*m-1) = K(3*pos(z)-eje,3*m-1) + ke(3*z-eje,14);
                    K(3*pos(z)-eje,3*m)   = K(3*pos(z)-eje,3*m)   + ke(3*z-eje,15);
                    K(3*pos(z)-eje,3*n-2) = K(3*pos(z)-eje,3*n-2) + ke(3*z-eje,16);
                    K(3*pos(z)-eje,3*n-1) = K(3*pos(z)-eje,3*n-1) + ke(3*z-eje,17);
                    K(3*pos(z)-eje,3*n)   = K(3*pos(z)-eje,3*n)   + ke(3*z-eje,18);
                    K(3*pos(z)-eje,3*o-2) = K(3*pos(z)-eje,3*o-2) + ke(3*z-eje,19);
                    K(3*pos(z)-eje,3*o-1) = K(3*pos(z)-eje,3*o-1) + ke(3*z-eje,20);
                    K(3*pos(z)-eje,3*o)   = K(3*pos(z)-eje,3*o)   + ke(3*z-eje,21);
                    K(3*pos(z)-eje,3*p-2) = K(3*pos(z)-eje,3*p-2) + ke(3*z-eje,22);
                    K(3*pos(z)-eje,3*p-1) = K(3*pos(z)-eje,3*p-1) + ke(3*z-eje,23);
                    K(3*pos(z)-eje,3*p)   = K(3*pos(z)-eje,3*p)   + ke(3*z-eje,24);
                end
            end
        end
        K_global = K;
    else
        fprintf('\n\nVerifique las dimensiones de los datos de entrada en la función\n\n')
    end
end