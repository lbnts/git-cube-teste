
function [Beci, Becef, Bned] = IGRF13 (lat, lon, h, J2000, ts)
% 
% Campo magnético terrestre, modelo IGRF-13
% 
% Inputs:
%   h: altitude [m]
%   lat: Latitude positive from equator [°]
%   lon: Longitude positive east from Greenwich [°]
%   days: Decimal days since January 1, 2020
%
% Outputs - magnetic field strength in local tangential coordinates
%   Beci: B in inertial system
%   Becef: B in ECEF system
%   Bned: B in NED system
% 
% http://hanspeterschaub.info/Papers/UnderGradStudents/MagneticField.pdf
% 

GpR = pi / 180; % graus -> radianos
RpG = 180 / pi; % radianos -> graus

J2000_0 = 7304.5; % Data J2000 em 1-jan-2020, 00:00:00

r = h * 1e-3 + 6371.2; % Geocentric radius [km]

days = J2000 + ts / 86400 - J2000_0; % número de dias desde 1-1-2020

% Campo magnético
[Br, Bt, Bp] = magnet (r, lat*RpG, lon*RpG, days);
[X, Y, Z] = msph2cart (Br, Bt, Bp);

Bned = [X ; Y ; Z]; % Campo magnético no sistema NED [nT]

% NED -> ECEF
R = [-cos(lon)*sin(lat) -sin(lon) -cos(lon)*cos(lat) ;
     -sin(lon)*sin(lat)  cos(lon) -sin(lon)*cos(lat) ;
     cos(lat) 0 -sin(lat)];

Becef = R * [X ; Y ; Z]; % [nT]

% ECEF -> ECI
g = gmst (J2000, ts) * GpR; % tempo sideral médio de Greenwich [rad]
T = [cos(g) -sin(g) 0 ;
     sin(g) cos(g) 0 ;
     0 0 1];

Beci = T * Becef; % [nT]
