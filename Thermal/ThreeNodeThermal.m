%Three Node Thermal Analysis
clear all
numorbits=1;
Torb=90*60;
timestep=1;


%%Makrolon GV 30 Characteristics(Interface from Ring to Service Module)
%%
kSR=0.16;%Thermal Conductivity W/mKZ
LSR=0.01;%thickness of the glass fiber insulation
ASR=0.051076;%Interface area of Ring and Service Module 226mm X 226mm

%%Alu Interface(Interface from Ring to QuadPack)
%%
hQR=150;
AQR=0.051076;%Interface area of Ring and Quadpack 226mm X 226mm


%%Ring Characteristics
%%
mr=15; %mass of the ring
cpr=910;%specific heat capacity of the ring
Ar=0.5787;%Area of the ring which radiates
epsilonR=0.05;%Emissivity of the Ring


%%
%SVM Characteristics
mSVM=40;%mass of the quadpack
cpSVM=200; %specific heat of the quadpack
%%
%QuadPack Characteristics
mQP=7.5;%mass of the quadpack
cpQP=910; %specific heat of the quadpack
epsilonQP=0.6402;%thermal emissivity of the quadpack(average value from QuadThermal all sides)
%%
Area1=272*272/(1000*1000);%Radiator Area, small panel ontop
alpha1SVM=0.248;
epsilon1SVM=0.888;

alpha1QP=0.65;
epsilon1QP=0.82;
%---------------------------------------------------------------------------------------------------
Area2=272*272/(1000*1000);%Perfectly Insulated Side
alpha2SVM=0;
epsilon2SVM=0;

alpha2QP=0;
epsilon2QP=0;
%---------------------------------------------------------------------------------------------------
Area3=272*391/(1000*1000);%Panel facing facing the sun such that its surface is perpendicular to the Z-axis
epsilon3=0.825;%From old SMAD pg 422

%---------------------------------------------------------------------------------------------------
Area4=272*391/(1000*1000);%Panel 4, Kapton Film
alpha4SVM=0.49;
epsilon4SVM=0.3;

alpha4QP=0.53;
epsilon4QP=0.82;
%---------------------------------------------------------------------------------------------------
Area5=272*391/(1000*1000);%Panel 5,
alpha5SVM=0.248;
epsilon5SVM=0.888;

alpha5QP=0.53;
epsilon5QP=0.82;

%---------------------------------------------------------------------------------------------------
Area6=272*391/(1000*1000);%Panel 6, Kapton Film
alpha6SVM=0.49;
epsilon6SVM=0.3;

alpha6QP=0.53;
epsilon6QP=0.82;
%
%%
%QuadPack Characteristics
mQP=7.5;%mass of the quadpack
cpQP=910; %specific heat of the quadpack
epsilonQP=0.6402;%thermal emissivity of the quadpack(average value from QuadThermal all sides)
AQP=0.4994;%area over which the quadpack radiates, Total area of outside quadpack except Side 2


%%
%Initial Conditions, Heat Input and Constants\

Tr_icold=310;%Initial temperature of the ring
Tq_icold=310;%Temperature of the Quadpack interface with the glass fiber, fixed for now
Tr_ihot=310;%Initial temperature of the ring
TSVM_ihot=310;%Temperature of the Quadpack interface with the glass fiber, fixed for now
TQP_ihot=310;
TQP_icold=310;
QinR=RingFlux((Torb/timestep),numorbits);%Input raidation heat to the ring, connect to Flux function later
%Its different because you divide by 6, so do not care about this
Qin=QuadFlux((Torb/timestep),numorbits);%Input radiation to the cubesat


boltz=5.67*10^-8;%boltzmann constant


%%
%%Propagator
time=[];%initialize time output array

Trcold=[];%initialize Ring Temperature output array
TSVMcold=[];%initialize Quadpack Temperature output array

Trhot=[];%initialize Ring Temperature output array
TSVMhot=[];%initialize Service Module Temperature output array
TQPhot=[];%initialize QuadPack Temperature output array
uhot_out=[];%
ucold_out=[];%
rad=[];
condRingQP=[];
Input=[];
i=1;
t=0;
diff=[];
while t<Torb*numorbits

    %Matrix computation for the worst case of the service module such that it
    %is pointing at earth when triangle sps-earth-sun is a right triangle
    Ahot=[mr*cpr/timestep,0,0;0,mSVM*cpSVM/timestep,0;0,0,mQP*cpQP/timestep];
    Bhot=[1,(mr*cpr/timestep)-(kSR*ASR/LSR)-5*(hQR*AQR),kSR*ASR/LSR,-boltz*Ar*epsilonR,0,0,5*hQR*AQR,0,0;0,(kSR*ASR/LSR),(mSVM*cpSVM/timestep)-(kSR*ASR/LSR),0,-boltz*(Area1*epsilon1SVM+Area2*epsilon2SVM+Area3*epsilon3+Area4*epsilon4SVM+Area5*epsilon5SVM+Area6*epsilon6SVM),1,0,0,0;0,hQR*AQR,0,0,0,0,(mQP*cpQP/timestep)-(hQR*AQR),-boltz*(Area1*epsilon1QP+Area2*epsilon2QP+Area3*epsilon3+Area4*epsilon4QP+Area5*epsilon5QP+Area6*epsilon6QP),1];
    uhot=[QinR(i);Tr_ihot;TSVM_ihot;Tr_ihot^4;TSVM_ihot^4;Qin(1,i);TQP_ihot;TQP_ihot^4;Qin(2,i)];
    uhot_out=[uhot_out,uhot];
    xhot=Ahot\(Bhot*uhot);

    
    %update of variables for next step of the hot case
    Tr_ihot=xhot(1);
    TSVM_ihot=xhot(2);
    TQP_ihot=xhot(3);
    %Generation of the output for the hot case
    Trhot=[Trhot,xhot(1)];
    TSVMhot=[TSVMhot,xhot(2)];
    TQPhot=[TQPhot,xhot(3)];
    rad=[rad,boltz*(Area1*epsilon1QP+Area2*epsilon2QP+Area3*epsilon3+Area4*epsilon4QP+Area5*epsilon5QP+Area6*epsilon6QP)*TQP_ihot^4];
    condRingQP=[condRingQP,(hQR*AQR)*(Tr_ihot-TQP_ihot)];
    %---------------------------------------------------------------------------------------------------------------------------------------------------------
    
 
    
    
    time=[time,t];

    
    
    %Time keeping and clock
    t=t+timestep;
    i=i+1
    

    
end

% figure
% plot(time/60/90,Trhot-273.5,'b');
% hold on
% plot(time/60/90,TSVMhot-273.5,'r')
% hold on
% plot(time/60/90,TQPhot-273.5,'k')
% hold off