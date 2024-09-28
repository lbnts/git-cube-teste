
function q = TransfElementar (a, e)
% Transformação elementar
%   a: ângulo [rad]
%   e: eixo

q = [cos(0.5*a) ;
     e(1) * sin(0.5*a) ;
     e(2) * sin(0.5*a) ;
     e(3) * sin(0.5*a)];