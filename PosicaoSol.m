
function Xsol = PosicaoSol (J2000)

M = 357.528 + 0.9856003 * J2000;
V = 280.461 + 0.9856474 * J2000 + 1.915 * sind(M) + 0.02 * sind(2*M);
e = 23.4393 + 0.0000004 * J2000;
Xsol = [cosd(V) ; cosd(e)*sind(V) ; sind(e)*sind(V)];