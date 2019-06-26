%==========================================================================
% MADE BY WILLEM VAN LYNDEN
%==========================================================================
%==========================================================================
% Modified simple sgp4 orbit propagation programn
%==========================================================================
% This programn propogates an orbit of a satellite, using a Two Line
% Element(TLE) input and a time (minutes since beginning of orbit) input.
% Outputs can be: visual orbit(s) of initial (ti) and final orbit (tf), including
% years to final orbit. Also possible is to propagate multiple orbit from a
% certain time (t1) to another (t2) (time in minutes)
clc
clear
format long g

%==========================================================================
% Set initial values
%==========================================================================
ge = 398600.8; % Earth gravitational constant
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
X = [];
Y = [];
Z = [];
position = [];
position2 = [];
i = 0;
j = 0;

%==========================================================================
% Set Time input of moments of interest IN MINUTES
%==========================================================================
tfe = 18900*100 + 45;   % final orbit time end
tfs = tfe - 100;      % final orbit time start
ti = 0; % initial orbit time

%==========================================================================
% Calculate output
%==========================================================================
while tfs < tfe
    tsince = tfs;
    [pos, vel] = sgp4(tsince, satdata);
    position = [position; [transpose(pos) i]];
    tfs = tfs + 1;
    i = i + 1;
end
while ti < 100 % initial orbit time period range
    tsince = ti;
    [pos, vel] = sgp4(tsince, satdata);
    position2 = [position2; [transpose(pos) j]];
    ti = ti + 1;
    j = j + 1;
end


% colormap(hsv)
% patch(position(:,1),position(:,2),position(:,3),position(:,4),'FaceColor','none','EdgeColor','interp')
% colorbar
% view(3)

%==========================================================================
% Give and plot output
%==========================================================================
figure;hold on
Re = 6371.000;
[x y z] = sphere;
a=[0 0 0 Re];
s1=surf(x*a(1,4)+a(1,1),y*a(1,4)+a(1,2),z*a(1,4)+a(1,3));
daspect([1 1 1])
plot3(position(:,1),position(:,2),position(:,3))
plot3(position2(:,1),position2(:,2),position2(:,3))
view(30,10)
rnew = sqrt(position(end,1)^2+position(end,2)^2+position(end,3)^2);
altitude = rnew - Re
years = tfe/(60*24*365.25)
