
function A = QuaternionParaMCD (q)
% Calcula a matriz de cossenos diretores a partir do quaternion

A = (q(4)^2 - q(1:3)' * q(1:3)) * eye(3) + 2 * q(1:3) * q(1:3)' ...
                                               - 2 * q(4) * Skew(q(1:3));