
function euler = QuaternionParaEuler (q)
% Obtém os ângulos de Euler a partir do quaternion de atitude
%   . sequência 3-2-1
%   . [roll pitch yaw] [rad]

euler = MCDparaEuler ( QuaternionParaMCD (q) );