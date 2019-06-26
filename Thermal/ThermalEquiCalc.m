%Thermal Preliminary Steady State Calculations

alpha=0.15;%absorp
beta=0.05;%emissivity
boltz=5.67*10^-8;%boltzmann cst

%%Ring Definition
%area that is exposed if the z-axis is directly aligned with the incoming radiation(8)
Aring1=0.02975;
alpha1_ring=0.15;
epsilon1_ring=0.05;

%inside area that is exposed for a adcs pointing offset
Rinner=0.4685;
heightIllu=0.450;
Aring2=heightIllu*2*Rinner;
alpha2_ring=0.15;
epsilon2_ring=0.05;

%outside area that is exposed for a adcs pointing offset
Router=0.4735;
PFactor=0.6718;
Aring3=heightIllu*2*Rinner*PFactor;
alpha3_ring=0.15;
epsilon3_ring=0.05;
    
%outside area that is radiating
Aring4=2*Aring3;
alpha4_ring=0.15;
epsilon4_ring=0.05;



qs=1365;
qIR=258;
albedo=0.35;
boltz=5.67*10^-8;%boltzmann cst
Re=6371000;
H=350*1000;
circ_ecli=(pi/2)+atan(Re/(Re+H));
Redmax=((Re^2))/((Re+800000)^2);
Redmin=((Re^2))/((Re+350000)^2);
adcsError=0.5;
    



F_IR=258;%Earth IR Flux
Re=6371*10^3;%Earth Radius
h=850*10^3;%Satellite Altitude
A_ir=0.28935;%Area on which infrared radiation is incident
Qir=((Re^2)*F_IR)/((Re+h)^2)*beta*A_ir;





%Q_int=43.5;%Internal Power Dissipation
Q_int=0;

A_rad=0.5787;%area to radiate

% Qinmax=qs*Aring1*alpha1_ring*cos(deg2rad(0.5))+qs*Aring2*alpha2_ring*cos(deg2rad(89.5))+qs*Aring3*alpha3_ring*cos(deg2rad(89.5))...
% +qIR*Redmax*epsilon3_ring*Aring3+...
% +qs*albedo*alpha3_ring*Aring3;
Qinmax=47.6653;
Qinmin=3.2834;
% Qinmin=Q_int+Qir;

Tmax=(Qinmax/(A_rad*(5.67*10^-8)*beta))^0.25-273.5%Equivalent Temperature Hot Case in °C
Tmin=(Qinmin/(A_rad*(5.67*10^-8)*beta))^0.25-273.5%Equivalent Temperature Cold Case in °C

