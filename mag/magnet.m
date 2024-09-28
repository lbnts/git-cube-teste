
function [Br, Bt, Bp] = magnet (r, theta, phi, days)
% Inputs
%   r: Geocentric radius
%   theta: Latitude measured in degrees positive from equator
%   phi: Longitude measured in degrees positive east from Greenwich
%   days: Decimal days since January 1, 2025
%
% Outputs - magnetic field strength in local tangential coordinates
%   Br: B in radial direction
%   Bt: B in theta direction
%   Bp: B in phi direction

% g and h Schmidt quasi-normalized coefficients
gn = [1 1 2 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 11 11 11 11 11 11 11 11 11 11 11 11 12 12 12 12 12 12 12 12 12 12 12 12 12 13 13 13 13 13 13 13 13 13 13 13 13 13 13];
gm = [0 1 0 1 2 0 1 2 3 0 1 2 3 4 0 1 2 3 4 5 0 1 2 3 4 5 6 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 8 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 10 0 1 2 3 4 5 6 7 8 9 10 11 0 1 2 3 4 5 6 7 8 9 10 11 12 0 1 2 3 4 5 6 7 8 9 10 11 12 13];
gvali = [-29404.8 -1450.9 -3749.4 5164.98 1452.32 3408 -7290.91 2393.89 415.6 3950.63 4479.76 337.7 -647.16 35.5 -1845.11 3692.5 1443.28 -662.16 -335.44 9.47 952.88 1238.15 1089.43 -1210.48 -197.54 31.41 -43.46 2161.09 -2720.52 -237.48 1157.03 195.11 39.52 -17.44 6.34 1191.48 650.2 -987.05 -20.71 -564.13 226.91 94.05 -41.36 -0.19 474.8 1070.19 315.09 -124.47 -62.01 -444.72 19.14 66.3 -24.03 -7.25 -342.81 -1508.37 -21.07 280.98 -105.18 51.74 -37.19 38.08 11.46 -6.37 -2.26 1033.35 -652.94 -1022.62 754.33 -215.56 47.53 -65.88 -4.96 31.86 -5.29 0.54 1.8 -1320.39 -89.7 397.56 843.99 -584.3 233.81 62.49 58.53 -17.56 -12.77 0.94 -3.06 -0.17 126.96 -1557.16 773.76 898.19 -295.24 556.7 0 213.71 0 27.23 2.84 5.02 -1.42 -0.22];
gsvi = [5.7 7.4 -16.5 -12.12 -1.82 5.5 -18.06 6 -9.49 -5.25 -8.85 -23.09 10.88 -3.77 -2.36 5.08 -4.61 0.94 2.88 0.63 -7.22 -5.67 5.98 12.95 -7.64 0 0.6 -2.68 -7.09 0 14.33 1.23 -3.09 -1.94 0.52 0 6.7 -5.61 16.57 -2.67 5.93 2.06 -0.25 0.25 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

hn = [1 1 2 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 11 11 11 11 11 11 11 11 11 11 11 11 12 12 12 12 12 12 12 12 12 12 12 12 12 13 13 13 13 13 13 13 13 13 13 13 13 13 13];
hm = [0 1 0 1 2 0 1 2 3 0 1 2 3 4 0 1 2 3 4 5 0 1 2 3 4 5 6 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 8 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 10 0 1 2 3 4 5 6 7 8 9 10 11 0 1 2 3 4 5 6 7 8 9 10 11 12 0 1 2 3 4 5 6 7 8 9 10 11 12 13];
hvali = [0 4652.5 0 -5181.6 -636.18 0 -251.38 468.44 -429.6 0 1560.03 -619.84 417.7 -258.61 0 484.95 1600.83 -570.39 71.66 69.38 0 -361.05 375.1 526.04 -351.97 20.71 45.74 0 -1826.68 -489.44 45.05 290.2 -13.58 -65.87 -1.17 0 563.06 -858.06 530.17 -312.81 220.98 24.71 -17.3 1.75 0 -2981.24 1195.15 813.23 -287.52 -212.25 135.7 3.01 -3.62 5.85 0 827.17 -42.14 595.01 560.98 -635.67 -4.13 -86.19 -27.82 -0.27 -5.22 0 0 1022.62 -196.78 -95.81 95.05 -18.82 -84.33 -36.42 -26.44 -5.44 -1.51 0 -1076.43 397.56 908.91 -876.45 33.4 166.64 -23.41 35.12 5.11 -8.49 0 0.28 0 -1557.16 928.51 1796.39 -393.65 -904.64 -45.15 80.14 -14.28 34.04 14.19 -4.01 -1.14 -0.33];
hsvi = [0 -25.9 0 -52.31 -19.4 0 18.37 -2.13 0.4 0 -0.55 25.44 7.53 -3.7 0 0 19.21 -2.82 6.66 0.21 0 0 -23.91 -12.95 4.37 0 0.67 0 21.28 17.38 -16.38 -2.47 -6.79 0.24 0.19 0 -13.41 33.65 -8.28 13.37 -4.45 -2.75 1.25 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% Checks to see if located at either pole to avoid singularities
if (theta > -0.00000001 && theta < 0.00000001)
    theta = 0.00000001;
elseif (theta < 180.00000001 && theta > 179.99999999)
    theta = 179.99999999;
end

% The angles must be converted from degrees into radians
theta = (90 - theta) * pi/180;
phi = phi * pi/180;

a = 6371.2; % Reference radius used in IGRF

N = max(gn);
g = zeros(N,N+1);
h = zeros(N,N+1);
for x = 1 : length(gn)
    g(gn(x),gm(x)+1) = gvali(x) + gsvi(x) * days / 365;
    h(hn(x),hm(x)+1) = hvali(x) + hsvi(x) * days / 365;
end

% Initialize each of the variables
% Br B in the radial driection
% Bt B in the theta direction
% Bp B in the phi direction
% 
% P The associated Legendre polynomial evaluated at cos(theta)
% The nomenclature for the recursive values generally follows
% the form P10 = P(n-1,m-0)
% dP The partial derivative of P with respect to theta

Br = 0; Bt = 0; Bp = 0;
P11 = 1; P10 = P11;
dP11 = 0; dP10 = dP11;
for m = 0 : N
    for n = 1 : N
        if m <= n
            % Calculate Legendre polynomials and derivatives recursively
            if n == m
                P2 = sin(theta) * P11;
                dP2 = sin(theta) * dP11 + cos(theta) * P11;
                P11 = P2; P10=P11; P20=0;
                dP11 = dP2; dP10 = dP11; dP20 = 0;
            elseif n == 1
                P2 = cos(theta) * P10;
                dP2 = cos(theta) * dP10 - sin(theta) * P10;
                P20 = P10; P10 = P2; 
                dP20 = dP10; dP10 = dP2;
            else
                K = ((n - 1)^2 - m^2) / ((2 * n - 1) * (2 * n - 3));
                P2 = cos(theta) * P10 - K * P20;
                dP2 = cos(theta) * dP10 - sin(theta) * P10 - K * dP20;
                P20 = P10; P10 = P2;
                dP20 = dP10; dP10 = dP2;
            end

            % Calculate Br, Bt, and Bp
            Br = Br + (a / r)^(n + 2) * (n+1) *  ...
               ((g(n,m+1) * cos(m * phi) + h(n,m+1) * sin(m * phi)) * P2);
            Bt = Bt + (a / r)^(n + 2) * ...
               ((g(n,m+1) * cos(m * phi) + h(n,m+1) * sin(m * phi)) * dP2);
            Bp = Bp + (a / r)^(n + 2) * ...
               (m * (-g(n,m+1) * sin(m * phi) + h(n,m+1) * cos(m * phi)) * P2);
        end
    end
end

% Br = Br;
Bt = -Bt;
Bp = -Bp / sin(theta);