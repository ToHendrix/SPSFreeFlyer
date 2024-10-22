%==========================================================================
% MADE BY WILLEM VAN LYNDEN
%==========================================================================
%==========================================================================
% Modified simple sgp4 orbit propagation programn
%==========================================================================
% This programn propogates an orbit of a satellite, using a Two Line
% Element(TLE) input and a time (minutes since beginning of orbit) input.
% Outputs are three-dimensional location and velocity in space
clc
clear
format long g

%==========================================================================
% Set initial values
%==========================================================================
ge = 398600.8; % Earth gravitational constant
Re = 6371.000;
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440.;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

%==========================================================================
% TLE file name
%==========================================================================
fname = 'correctedtle.txt';

%==========================================================================
% Open the TLE file and read TLE elements
%==========================================================================
fid = fopen(fname, 'r');
% 19-32	04236.56031392	Element Set Epoch (UTC)
% 3-7	25544	Satellite Catalog Number
% 9-16	51.6335	Orbit Inclination (degrees)
% 18-25	344.7760	Right Ascension of Ascending Node (degrees)
% 27-33	0007976	Eccentricity (decimal point assumed)
% 35-42	126.2523	Argument of Perigee (degrees)
% 44-51	325.9359	Mean Anomaly (degrees)
% 53-63	15.70406856	Mean Motion (revolutions/day)
% 64-68	32890	Revolution Number at Epoch

%==========================================================================
% Main programn
%==========================================================================

% Read TLE lines
while (1)
    % read first line
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    Cnum = tline(3:7);      			        % Catalog Number (NORAD)
    SC   = tline(8);					        % Security Classification
    ID   = tline(10:17);			            % Identification Number
    epoch = str2num(tline(19:32));              % Epoch
    TD1   = str2num(tline(34:43));              % first time derivative
    TD2   = str2num(tline(45:50));              % 2nd Time Derivative
    ExTD2 = tline(51:52);                       % Exponent of 2nd Time Derivative
    BStar = str2num(tline(54:59));              % Bstar/drag Term
    ExBStar = str2num(tline(60:61));            % Exponent of Bstar/drag Term
    BStar = BStar*1e-5*10^ExBStar;
    Etype = tline(63);                          % Ephemeris Type
    Enum  = str2num(tline(65:end));             % Element Number
    
    % read second line
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    i = str2num(tline(9:16));                   % Orbit Inclination (degrees)
    raan = str2num(tline(18:25));               % Right Ascension of Ascending Node (degrees)
    e = str2num(strcat('0.',tline(27:33)));     % Eccentricity
    omega = str2num(tline(35:42));              % Argument of Perigee (degrees)
    M = str2num(tline(44:51));                  % Mean Anomaly (degrees)
    no = str2num(tline(53:63));                 % Mean Motion
    a = ( ge/(no*2*pi/86400)^2 )^(1/3);         % semi major axis (m)
    rNo = str2num(tline(64:68));                % Revolution Number at Epoch
end
fclose(fid);

% Calculate necessary values
satdata.epoch = epoch;
satdata.norad_number = Cnum;
satdata.bulletin_number = ID;
satdata.classification = SC; % almost always 'U'
satdata.revolution_number = rNo;
satdata.ephemeris_type = Etype;
satdata.xmo = M * (pi/180);
satdata.xnodeo = raan * (pi/180);
satdata.omegao = omega * (pi/180);
satdata.xincl = i * (pi/180);
satdata.eo = e;
satdata.xno = no * TWOPI / MINUTES_PER_DAY;
satdata.xndt2o = TD1 * 1e-8 * TWOPI / MINUTES_PER_DAY_SQUARED;
satdata.xndd6o = TD2 * TWOPI / MINUTES_PER_DAY_CUBED;
satdata.bstar = BStar;

%==========================================================================
% Set Time input of moment of interest IN MINUTES
%==========================================================================
tsince = 27000*100 -18;

%==========================================================================
% Calculate output
%==========================================================================
[pos, vel] = sgp4(tsince, satdata);
rnew = sqrt(pos(1)^2+pos(2)^2+pos(3)^2);
alt = rnew - Re;
years = tsince/(60*24*365.25);

%==========================================================================
% Give output
%==========================================================================
fprintf('     TSINCE            radius              X                Y                Z     [km]\n');
fprintf(' %9.1f%22.8f%18.8f%18.8f%18.8f \n',tsince,rnew,pos(1),pos(2),pos(3));
fprintf('     Years             altitude            XDOT             YDOT             ZDOT    [km/s]\n');
fprintf(' %9.1f%22.8f%18.8f%18.8f%18.8f \n',years,alt,vel(1),vel(2),vel(3));