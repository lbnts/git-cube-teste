clear; close all; clc;
addpath ./mag
addpath ./quat

%%% ATUALIZACAO HEAD GIT HUB CODE %%%
%%%%%%%%% NOTE A DIFERENCA %%%%%%%%%%
%%%%%%%%%%%% DIFF %%%%%%%%%%%%%%%%%%%
%%%%%%%  NOVA ATT %%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%
%%%     Constantes     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

GpR = pi / 180; % graus para radianos
RpG = 180 / pi; % radianos para graus

% WGS-84
RT = 6378137;        % raio da Terra [m]
F = 1/298.257223563; % achatamento
mu = 3.986004418e14; % constante gravitacional [m^3/s^2]

%%% --- Constantes --- %%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Propriedades do CubeSat     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dimensões
M = 3.3;             % massa [kg]
dim = [0.1 0.1 0.3]; % lados [m]
% Matriz de inércias [kg.m^2]
J = diag([(M/12)*(dim(2)^2+dim(3)^2) 
          (M/12)*(dim(1)^2+dim(3)^2) 
          (M/12)*(dim(1)^2+dim(2)^2)]);

% Momento magnético residual
m_res = [0 0 0.0005]'; % momento de dipolo residual [Am^2]

%%% --- Propriedades do CubeSat --- %%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Definição da órbita     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data J2000 - 'sábado, 26/03/2022 14:12:20';
J2000_0 = 8120.5;              % Data J2000 às 0h00min00s
ts = (14 * 60 + 12) * 60 + 20; % Segundos do dia [s]

el(1) = 500e3 + RT;            % Semi eixo maior [m]
el(2) = 0;                     % Excentricidade
el(3) = 35 * GpR;              % Inclinação [rad]
el(4) = 30 * GpR;              % Ascensão reta do nodo ascendente [rad]
el(5) = 10 * GpR;              % Argumento do perigeu [rad]
el(6) = 0 * GpR;               % Anomalia média na época [rad]

n = sqrt(mu / el(1)^3);        % Movimento médio [rad/s]
P = 2 * pi / n;                % Período orbital [s]

%%% --- Definição da órbita --- %%%

%% Condições iniciais
t = 0; % tempo de simulação [s]
wioo = [0 -n 0]'; % w_io^o [rad/s]
wibb = [0.087 -0.192 0.052]'; % w_ib^b [rad/s]
qob = [1 0 0 0]' ./ norm([1 0 0 0]);
wobb = wibb - TransfVetor (wioo, qob); % w_ob^b [rad/s]

%% Parâmetros do PD
kp = 3e-7; % ganho proporcional
kd = 5e-5; % ganho derivativo
wref = [0 0 0]';   % velocidade angular de referência [rad/s]
qref = [1 0 0 0]'; % quaternion de referência

%% Inicializa os vetores de armazenamento
tempo = t;
gq(1,:) = qob;
gwo(1,:) = wibb - TransfVetor (wioo, qob);
[gXi(1,:), ~] = ProblemaDireto (el, t);

Ba = [0 0 0]';

%% Simulação
K = 40000;         % número de passos
h = 4*P / (K - 1); % tamanho do passo [s]

for k = 2 : K
    
    t = t + h;   % tempo de simulação [s]
    ts = ts + h; % tempo em segundos do dia [s]

    %% Órbita
    [Xi, Vi, qoi] = ProblemaDireto (el, t); 

    % ECI -> ECEF
    g = gmst (J2000_0, ts) * GpR; % tempo sideral médio de Greenwich [rad]
    qie = TransfElementar (g, [0 0 1]);
    Xe = TransfVetor (Xi, qie); % posição no sistema ECEF [m]
    L = ECEFparaGeod (Xe); % L = [lat lon alt], [rad rad m]

    %% Campo magnético IGRF13
    [Beci, Becef, Bned] = IGRF13 (L(1), L(2), L(3), J2000_0, ts);
    Bo = TransfVetor (Beci, [qoi(1) -qoi(2) -qoi(3) -qoi(4)]'); % no sistema orbital [nT]
    Bb = TransfVetor (Bo, qob) .* 1e-9; % no sistema do corpo [T]

    %% Detumbling
    % dB = Bb - Ba;
    % mc = 12000 * dB / h; % momento magnético [Am^2]
    % Nc = Skew(Bb) * mc; % torque [Nm]
    % Ba = Bb; % salva para a próxima iteração

    %% Controle PD
    dw = wref - wobb; % erro de velocidade angular
    dq = MultQuat (qref, [qob(1) -qob(2) -qob(3) -qob(4)]'); % erro de quaternions
    Nd = -(kp * dq(2:4) + kd * dw); % Torque de controle desejado [Nm]

    % Projeção do torque, eqs. (4.9) a (4.11) de 
    %   https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/260833
    mc = Skew(Bb) * Nd / norm(Bb)^2; % momento magnético [Am^2]
    Nc = Skew(Bb) * mc; % torque [Nm]

    % Torques de perturbação
    [Ngg, Naero, Nsol] = Perturbacoes (qob, n, J, dim, Xi, Vi, J2000_0+t/86400);
    N = Nc + Ngg + Naero + Nsol + Skew(m_res) * Bb; % Torque total
    
    %% Integra a dinâmica e cinemática (RK4)

    % Dinâmica 
    k1 = dinamica (wibb, N, J);
    k2 = dinamica (wibb + 0.5 * h * k1, N, J);
    k3 = dinamica (wibb + 0.5 * h * k2, N, J);
    k4 = dinamica (wibb + h * k3, N, J);

    wibb = wibb + h * (k1 + 2 * (k2 + k3) + k4) / 6;
    wobb = wibb - TransfVetor (wioo, qob);

    % Cinemática
    k1 = cinematica (qob, wobb);
    k2 = cinematica (qob + 0.5 * h * k1, wobb);
    k3 = cinematica (qob + 0.5 * h * k2, wobb);
    k4 = cinematica (qob + h * k3, wobb);

    qob = qob + h * (k1 + 2 * (k2 + k3) + k4) / 6;    

    %% Armazena valores para os gráficos
    tempo(k) = t;
    gwo(k,:) = wobb;
    gq(k,:) = qob;
    gEob(k,:) = QuaternionParaEuler (qob);
    % gXi(k,:) = Xi;
    % gB(k,:) = Bb;
    % gN(k,:) = N;
    gmc(k,:) = mc;
    
end

plot(tempo/P, gq)