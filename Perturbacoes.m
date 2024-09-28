
function [Ngg, Naero, Nsol] = Perturbacoes (q, n, J, dim, Xi, Vi, J2000)
% q = qob

% Torque de gradiente de gravidade
Ngg = 6*n^2 * [(J(3,3)-J(2,2)) * (q(3)*q(4) + q(2)*q(1)) * (-q(2)^2 - q(3)^2 + q(4)^2 + q(1)^2) ;
    (J(1,1)-J(3,3)) * (q(2)*q(4) - q(3)*q(1)) * (-q(2)^2 - q(3)^2 + q(4)^2 + q(1)^2) ;
    (J(2,2)-J(1,1)) * (q(3)*q(4) + q(2)*q(1)) * (q(2)*q(4) - q(3)*q(1))];

% Torque aerodinâmico
rcmp = [0.005 ; 0.004 ; 0.006]; % distância do cg para o cp [m]
rho = 4.89e-13; % densidade do ar a 500 km [kg/m^3]
Cd = 2.0; % coeficiente de arrasto
Faero = -0.5 * Cd * rho * (dim(2)*dim(3)) * norm(Vi)^2; % força aerodinâmica [N]
Naero = -abs(Faero) * [2*(q(2)*q(4) + q(3)*q(1))*rcmp(2) - 2*(q(2)*q(3) - q(4)*q(1))*rcmp(3) ;
        (q(2)^2 - q(3)^2 - q(4)^2 + q(1)^2)*rcmp(3) - 2*(q(2)*q(4) + q(3)*q(1))*rcmp(1) ;
        2*(q(2)*q(3) - q(4)*q(1))*rcmp(1) - (q(2)^2 - q(3)^2 - q(4)^2 + q(1)^2)*rcmp(2)];


% Torque de pressão de radiação
c = 2.99792458e8; % velocidade da luz [m/s]
EAM0 = 1366.9; % irradiância solar a 1 UA [W/m^2]
Cp = 1.6; % constante de absorção
Fsol = -Cp * (EAM0 / c) * (dim(2)*dim(3)); % força de pressão de radiação [N]
Nsol = Fsol * Skew(rcmp) * (Xi/norm(Xi) - PosicaoSol (J2000)) ./ norm(Xi/norm(Xi) - PosicaoSol (J2000));

