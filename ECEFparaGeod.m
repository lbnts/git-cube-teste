function L = ECEFparaGeod (x)

a = 6378137; % [m]
f = 1 / 298.257222101;
b = a * (1 - f);
e = sqrt(1 - b^2/a^2);

tol = 1e-10;

h = 0; ha = 2 * tol;
lat = atan2( x(3) , (sqrt(x(1)^2 + x(2)^2) * (1 - e^2)) ); lata = lat + 2 * tol;
lon = atan2(x(2) , x(1));

while (abs((lat - lata) / lat) > tol && abs((h - ha) / h) > tol)
    
    ha = h;
    lata = lat;
    
    Rn = a / sqrt(1 - e^2 * sin(lat)^2);
    h = sqrt(x(1)^2 + x(2)^2) / cos(lat) - Rn;
    lat = atan( (x(3) / sqrt(x(1)^2 + x(2)^2)) * ( (Rn + h) / (Rn * (1 - e^2) + h) ) );    
     
end

L = [lat lon h];