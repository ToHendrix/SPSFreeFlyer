%Thermal Analysis for one Quadpack 272*272*391, Solar Panels Body Fixed
Quadwidth=272;
Quadlength=391;
QuadHeight=272;
%---------------------------------------------------------------------------------------------------
Area1=272*272/(1000*1000);%Radiator Area, small panel ontop
alpha1=0.248;
epsilon1=0.888;

%---------------------------------------------------------------------------------------------------
Area2=272*272/(1000*1000);%Perfectly Insulated Side
alpha2=0;
epsilon2=0;

%---------------------------------------------------------------------------------------------------
Area3=272*391/(1000*1000);%Panel facing facing the sun such that its surface is perpendicular to the Z-axis
alpha3=0.805;%From old SMAD pg 422
epsilon3=0.825;%From old SMAD pg 422
eff_sp=0.33;%panel efficiency
packing_sp=0.8;%packing density of solar panels
Kdiss=1;%Electronics Heat Dissipation factor

%---------------------------------------------------------------------------------------------------
Area4=272*391/(1000*1000);%Panel 4, Kapton Film
alpha4=0.49;
epsilon4=0.3;

%---------------------------------------------------------------------------------------------------
Area5=272*391/(1000*1000);%Panel 5,
alpha5=0.248;
epsilon5=0.888;
%epsilonopen=0.618;
%epsilonclosed=0.0998;
%---------------------------------------------------------------------------------------------------
Area6=272*391/(1000*1000);%Panel 6, Kapton Film
alpha6=0.49;
epsilon6=0.3;
%---------------------------------------------------------------------------------------------------
EffectiveAbsorption=[alpha1*Area1,alpha3*Area3,alpha4*Area4,alpha5*Area5,alpha6*Area6];
EffectiveEmmision=[epsilon1*Area1,epsilon3*Area3,epsilon4*Area4,epsilon5*Area5,epsilon6*Area6];
%---------------------------------------------------------------------------------------------------
qs=1365;
qIR=258;
albedo=0.35;
boltz=5.67*10^-8;%boltzmann cst
T_orb=91*60;
T_ecl=35*60;

%---------------------------------------------------------------------------------------------------
%Solar Panel Power Input of one quadpack
Qs=(alpha3-(eff_sp*packing_sp))*Area3*qs;%Power dissipated in by the solar panels as heat

Qelec=50;%Electrical Power generated by the solar panels while in the sun
Avg_Power=(T_orb-T_ecl)*Qelec/T_orb;% Average power usage per orbit per quadpack
Q_int=Avg_Power*Kdiss; %Power dissipated as heat internal to the quadpacks

%Earth IR and Albedo Input 
Q_earth_max=epsilon1*Area1*qIR;
Q_earth_min=epsilon4*Area4*qIR;

Q_alb_max=alpha6*Area6*albedo*qs;


%Energy Balance
Power_In_Max=Qs+81+Q_earth_max+Q_alb_max;
Power_In_Min=43+Q_earth_min;

EffectRadi_open=boltz*(epsilon1*Area1+epsilon3*Area3+epsilon4*Area4+epsilon5*Area5+epsilon6*Area6);
EffectRadi_closed=boltz*(epsilon1*Area1+epsilon3*Area3+epsilon4*Area4+epsilon5*Area5+epsilon6*Area6);

Tmax=(Power_In_Max/(EffectRadi_open))^0.25-273.5;
Tmin=((Power_In_Min)/(EffectRadi_closed))^0.25-273.5;


