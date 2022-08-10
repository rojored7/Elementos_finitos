function K = FEM_EV_Tetra_StiffnessMatrix(E,v,ui,uj,uk,ul)


%% MÉTODO DE ELEMENTOS FINITOS (MEF o FEM): Elemento con Volumen (Tetraedro) - Matriz de Rigidez

%% Variables

% Entrada:

% E                 ::>  Módulo de Young (elasticidad) del elemento.
% v                 ::>  Coeficiente de poisson.
% ui                ::>  Posición nodo i (X,Y,Z)
% uj                ::>  Posición nodo j (X,Y,Z)
% uk                ::>  Posición nodo k (X,Y,Z)
% ul                ::>  Posición nodo l (X,Y,Z)

% Salida:

% K                 ::>  Matriz de Rigidez


%% Programa
cant_in = nargin;
ok = 0;

if cant_in == 6
    if length(ui) == 3 && length(uj) == 3 && length(uk) == 3 && length(ul) == 3
        ok = 1;
    end
else
    fprintf('\n\nVerifique la cantidad de datos de entrada en la función\n\n')
end

if ok == 1
    Xi = ui(1);
    Yi = ui(2);
    Zi = ui(3);
    Xj = uj(1);
    Yj = uj(2);
    Zj = uj(3);
    Xk = uk(1);
    Yk = uk(2);
    Zk = uk(3);
    Xl = ul(1);
    Yl = ul(2);
    Zl = ul(3);
    V = det([1 Xi Yi Zi; 
             1 Xj Yj Zj; 
             1 Xk Yk Zk; 
             1 Xl Yl Zl])/6;
    m_beta1 = [1 Yj Zj ; 1 Yk Zk ; 1 Yl Zl];
    m_beta2 = [1 Yi Zi ; 1 Yk Zk ; 1 Yl Zl];
    m_beta3 = [1 Yi Zi ; 1 Yj Zj ; 1 Yl Zl];
    m_beta4 = [1 Yi Zi ; 1 Yj Zj ; 1 Yk Zk];
    m_gamma1 = [1 Xj Zj ; 1 Xk Zk ; 1 Xl Zl];
    m_gamma2 = [1 Xi Zi ; 1 Xk Zk ; 1 Xl Zl];
    m_gamma3 = [1 Xi Zi ; 1 Xj Zj ; 1 Xl Zl];
    m_gamma4 = [1 Xi Zi ; 1 Xj Zj ; 1 Xk Zk];
    m_delta1 = [1 Xj Yj ; 1 Xk Yk ; 1 Xl Yl];
    m_delta2 = [1 Xi Yi ; 1 Xk Yk ; 1 Xl Yl];
    m_delta3 = [1 Xi Yi ; 1 Xj Yj ; 1 Xl Yl];
    m_delta4 = [1 Xi Yi ; 1 Xj Yj ; 1 Xk Yk];
    beta1 = -1*det(m_beta1);
    beta2 =    det(m_beta2);
    beta3 = -1*det(m_beta3);
    beta4 =    det(m_beta4);
    gamma1 =    det(m_gamma1);
    gamma2 = -1*det(m_gamma2);
    gamma3 =    det(m_gamma3);
    gamma4 = -1*det(m_gamma4);
    delta1 = -1*det(m_delta1);
    delta2 =    det(m_delta2);
    delta3 = -1*det(m_delta3);
    delta4 =    det(m_delta4);
    B1 = [ beta1      0      0;
               0 gamma1      0; 
               0      0 delta1; 
          gamma1  beta1      0; 
               0 delta1 gamma1; 
          delta1      0  beta1];
    B2 = [ beta2      0      0;
               0 gamma2      0; 
               0      0 delta2; 
          gamma2  beta2      0;
               0 delta2 gamma2;
          delta2      0  beta2];
    B3 = [ beta3      0      0;
               0 gamma3      0;
               0      0 delta3;
          gamma3  beta3      0;
               0 delta3 gamma3;
          delta3      0  beta3];
    B4 = [ beta4      0      0;
               0 gamma4      0;
               0      0 delta4;
          gamma4  beta4      0;
               0 delta4 gamma4;
          delta4      0  beta4];
    B = [B1 B2 B3 B4]/(6*V);
    D = (E/((1+v)*(1-2*v)))*[1-v   v   v         0         0         0;
                               v 1-v   v         0         0         0;
                               v   v 1-v         0         0         0;
                               0   0   0 (1-2*v)/2         0         0;
                               0   0   0         0 (1-2*v)/2         0;
                               0   0   0         0         0 (1-2*v)/2];
    K = V*B'*D*B;
end

