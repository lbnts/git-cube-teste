
function G = gmst (J0, ts)
% 
% Tempo sideral médio de Greenwich
% 
% Entradas:
%   J0: Data J2000 às 00:00:00
%   ts: segundos do dia [s]
% Saída:
%   G: Tempo sideral médio de Greenwich [º]
% 

T = J0 / 36525;
G = mod(100.4606184 + 36000.77005361 * T + 0.00038793 * T^2 + ...
                                       0.2506844773374215 * ts / 60, 360);

end
