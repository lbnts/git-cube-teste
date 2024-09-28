
function [Xi, Vi, qoi] = ProblemaDireto (el, t)

mu = 3.986004418e14; % [m^3/s^2]

% Elementos orbitais
a = el(1);  % [m]
i = el(3);  % [rad]
W = el(4);  % [rad]
e = el(2); 
w = el(5);  % [rad]
M0 = el(6); % [rad]

% Movimento médio [rad/s]
n = sqrt(mu / el(1)^3);       

% Calcular a anomalia média
M = M0 + n*t;

% Resolver a equação de Kepler
u = kepler (mod(M, 2*pi), e);

% Calcular a distância geocêntrica
r = a * (1 - e * cos(u));

% Calcular as coordenadas no plano orbital
x = a * (cos(u) - e);
y = a * sin(u) * sqrt(1 - e^2);
z = 0;
vx = -(n * a^2 / r) * sin(u);
vy = (n * a^2 / r) * cos(u) * sqrt(1 - e^2);
vz = 0;

% Calcular o quaternion de transformação da orientação da órbita
%   plano orbital -> inercial
q1 = TransfElementar (-w, [0 0 1]);
q2 = TransfElementar (-i, [1 0 0]);
q3 = TransfElementar (-W, [0 0 1]);

qpoi = MultQuat (MultQuat (q1, q2), q3);

% Transformar posição e velocidade para o sistema inercial
Xi = TransfVetor ([x ; y ; z], qpoi);
Vi = TransfVetor ([vx ; vy ; vz], qpoi);

% Quaternion de transformação do sistema orbital para o inercial 
%   Fundamentals of Spacecraft Attitude Determination and Control, seção 2.6.4
%   https://link.springer.com/book/10.1007/978-1-4939-0802-8
o3 = -Xi ./ norm(Xi);
o2 = -Skew(Xi) * Vi ./ norm(Skew(Xi) * Vi);
o1 = Skew(o2) * o3;
C = [o1 o2 o3];
qoi = MCDparaQuat (C);

end

function u = kepler (M, e)
%%
% Resolve a equação de Kepler
% 
% A Practical Method for Solving the Kepler Equation
% Marc A. Murison
% http://www.alpheratz.net/murison/dynamics/twobody/KeplerIterations summary.pdf
% 
%% 

tol = 1.0e-14;

% Mnorm = mod(M, 2 * pi);

t34 = e * e;
t35 = e * t34;
t33 = cos(M);
u0 = M + (-0.5 * t35 + e + (t34 + 1.5 * t33 * t35) * t33) * sin(M);

u = u0;
dE = tol + 1;
count = 0; 

while dE > tol 
    
    t1 = cos(u0);
    t2 = e * t1 - 1;
    t3 = sin(u0);
    t4 = e * t3;
    t5 = -u0 + t4 + M;
    t6 = t5 / (0.5 * t5 * t4 / t2 + t2);
    eps3 = t5 / ((0.5 * t3 - 1/6 * t1 * t6) * e * t6 + t2);
    
    u = u0 - eps3;
    dE = abs(u - u0);
    u0 = u;
    count = count + 1;
    if count == 100 
        disp('Astounding! KeplerSolve failed to converge!');
        break;
    end

end

end % functon kepler

