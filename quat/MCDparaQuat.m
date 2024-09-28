
function q = MCDparaQuat (C)
% Transforma matriz de cossenos diretores para quaternions

[~, i] = max([abs(C(1,1)) abs(C(2,2)) abs(C(3,3)) abs(trace(C))]);
switch i
    case 1 
        q = [1 + 2*C(1,1) - trace(C) ;
             C(1,2) + C(2,1) ;
             C(1,3) + C(3,1) ;
             C(2,3) - C(3,2)];        
    case 2
        q = [C(2,1) + C(1,2) ;
             1 + 2*C(2,2) - trace(C) ;
             C(2,3) + C(3,2) ;
             C(3,1) - C(1,3)];
    case 3
        q = [C(3,1) + C(1,3) ;
             C(3,2) + C(2,3) ; 
             1 + 2*C(3,3) - trace(C) ;
             C(1,2) - C(2,1)];
    case 4
        q = [C(2,3) - C(3,2) ;
             C(3,1) - C(1,3) ; 
             C(1,2) - C(2,1) ;
             1 + trace(C)]; 
    otherwise
        q = [1 0 0 0]';
end

q = q ./ norm(q);