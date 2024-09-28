
function w = TransfVetor (v, q)
% Transforma vetor
%   w = q* x v x q

tmp = MultQuat ([q(1) ; -q(2) ; -q(3) ; -q(4)], [0 ; v(1) ; v(2) ; v(3)]);
w = MultQuat (tmp, q);
w = w(2:4);