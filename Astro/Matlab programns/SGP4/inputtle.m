% TLE input for a hypthetical satellite

% Set values%
Re = 6378.137*10^3;
ph = 1350*10^3;
ah = 9850*10^3;
mu = 3.986004415*10^14;
Cd = 2.2;
m = 250.;
A = 1.;
Cdc3 = 2.0
mc3 = 3
Ac3 = 3.265*10^(-2)
tday = 24*60*60;
rho0 = 0.157
Bstarc3 = 0.73533*10^(-4)
% Calculate other starting values
perigee = ph + Re;
apogee = ah + Re;
a = (perigee + apogee)/2;
e = 1. - perigee/a

bcoeffc3 = (mc3/(Cdc3*Ac3))
rho0true = 2*Bstarc3/(1/bcoeffc3)
% Calculate basic parameters
T = 2*pi*sqrt((a^3)/(mu))
v = sqrt(mu*((2/(a*(1-e)))-1/a));
n = tday/T                              %revolutions per day

bcoeff =(m/(Cd*A))
Bstar = rho0true*(1/bcoeff)/2
% Calculate inclination for sun-sync
imax = 108.;
imin = 98.;
incl = (180/pi)*acos(-((a/12352000)^(7./2.)))



