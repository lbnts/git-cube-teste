
function euler = MCDparaEuler (A)
% Obtém os ângulos de Euler a partir da matriz de cossenos diretores
%   . sequência 3-2-1
%   . [roll pitch yaw] [rad]

euler = [atan2(A(2,3), A(3,3)) ;
         asin(-A(1,3)) ;
         atan2(A(1,2), A(1,1))];
