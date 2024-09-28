clear; close all; clc;
% Teste de convenções:
% q = [qw qx qy qz]
%   qw: componente escalar
%   qx, qy, qz: componentes vetoriais
% 
% Conveção de Hamilton: ij = k, jk = i, ki = j

%% Transformação A -> B

vA = [1 2 3]'; % vetor em A

a = pi/6; % ângulo [rad]
T1 = [cos(a) sin(a) 0 ; -sin(a) cos(a) 0 ; 0 0 1]; % matriz de transformação

e = [0 0 1]'; % eixo
% quaternion de transformação
q1 = [cos(a/2) ; e(1)*sin(a/2) ; e(2)*sin(a/2) ; e(3)*sin(a/2)]; 
% matriz de transformação
Tq1 = eye(3) - 2*q1(1)*Skew([q1(2) q1(3) q1(4)]) + 2 * (Skew([q1(2) q1(3) q1(4)]))^2;

disp('Vetor em B:')
disp('vB = T1 * vA'); disp(T1*vA);
disp('vB = Tq1 * vA'); disp(Tq1*vA);
vB = MultQuat (MultQuat ([q1(1) -q1(2) -q1(3) -q1(4)]', [0 ; vA]), q1);
vB = vB(2:4);
disp('vB = q1^* * vA * q1'); disp(vB);

%% B -> C

a = pi/4; % ângulo [rad]
T2 = [1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a)];
e = [1 0 0]';
q2 = [cos(a/2) ; e(1)*sin(a/2) ; e(2)*sin(a/2) ; e(3)*sin(a/2)];
Tq2 = eye(3) - 2*q2(1)*Skew([q2(2) q2(3) q2(4)]) + 2 * (Skew([q2(2) q2(3) q2(4)]))^2;

disp('Vetor em C:')
disp('vC = T2 * vB'); disp(T2*vB);
disp('vC = Tq2 * vB'); disp(Tq2*vB);
vC = MultQuat (MultQuat ([q2(1) -q2(2) -q2(3) -q2(4)]', [0 ; vB]), q2);
vC = vC(2:4);
disp('vC = q2^* * vB * q2'); disp(vC);

%% A -> C
T3 = T2 * T1;
Tq3 = Tq2 * Tq1;
q3 = MultQuat (q1, q2); % q3 = q1 x q2

disp('Vetor em C:')
disp('vC = T3 * vA'); disp(T3*vA);
disp('vC = Tq3 * vA'); disp(Tq3*vA);
vC = MultQuat (MultQuat ([q3(1) -q3(2) -q3(3) -q3(4)]', [0 ; vA]), q3);
vC = vC(2:4);
disp('vC = q3^* * vA * q3'); disp(vC);
