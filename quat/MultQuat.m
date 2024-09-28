
function r = MultQuat (p, q)
% Multiplica quaternions
% r = p x q

r = [p(1) -p(2) -p(3) -p(4) ;
     p(2)  p(1) -p(4)  p(3) ;
     p(3)  p(4)  p(1) -p(2) ;
     p(4) -p(3)  p(2)  p(1)] * q;