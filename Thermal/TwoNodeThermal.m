%Two Node Thermal Analysis
clear all
numorbits=2;
Torb=90*60;
timestep=0.1;
Tref=30+273.5;
HeaterPow=0;


%%Glass Fiber Characteristics
%%
k=0.16;%Thermal Conductivity W/mKZ
L=0.01;%thickness of the glass fiber insulation
Aint=0.051076;%Interface area of glass fiber and ring and glass fiber and quadpacks 226mm X 226mm
mgf=Aint*L*2*1.9;%mass of the glass fiber insulator plate
cpgf=800;%specific heat capacity of glass fiber


%%Ring Characteristics
%%
mr=40; %mass of the ring
cpr=910;%specific heat capacity of the ring
Ar=0.5787;%Area of the ring which radiates
epsilonR=0.05;%Emissivity of the Ring


%%
%Quadpack Characteristics
mq=40;%mass of the quadpack
cpq=200; %specific heat of the quadpack
epsilonq=0.6402;%thermal emissivity of the quadpack(average value from QuadThermal all sides)
Aq=0.4994;%area over which the quadpack radiates, Total area of outside quadpack except Side 2

%%
%Initial Conditions, Heat Input and Constants\

Tr_icold=310;%Initial temperature of the ring
Tq_icold=310;%Temperature of the Quadpack interface with the glass fiber, fixed for now
Tr_ihot=310;%Initial temperature of the ring
Tq_ihot=310;%Temperature of the Quadpack interface with the glass fiber, fixed for now

QinR=RingFlux((Torb/timestep),numorbits);%Input raidation heat to the ring, connect to Flux function later
%Its different because you divide by 6, so do not care about this
QinQ=QuadFlux((Torb/timestep),numorbits);%Input radiation to the cubesat


boltz=5.67*10^-8;%boltzmann constant


%%
%%Propagator
time=[];%initialize time output array

Trcold=[];%initialize Ring Temperature output array
Tqcold=[];%initialize Quadpack Temperature output array

Trhot=[];%initialize Ring Temperature output array
Tqhot=[];%initialize Quadpack Temperature output array

uhot_out=[];%
ucold_out=[];%
rad=[];
Input=[];
i=1;
t=0;
diff=[];
while t<Torb*numorbits

    %Matrix computation for the hot case of the service module such that it
    %is pointing at earth when triangle sps-earth-sun is a right triangle
    Ahot=[mr*cpr/timestep,0;0,mq*cpq/timestep];
    Bhot=[1,(mr*cpr/timestep)-(k*Aint/L),k*Aint/L,-boltz*Ar*epsilonR,0,0;0,(k*Aint/L),mq*cpq/timestep-k*Aint/L,0,-boltz*Aq*epsilonq,1];
    uhot=[QinR(i);Tr_ihot;Tq_ihot;Tr_ihot^4;Tq_ihot^4;QinQ(2,i)];
    uhot_out=[uhot_out,uhot];
    xhot=Ahot\(Bhot*uhot);

    
    %update of variables for next step of the hot case
    Tr_ihot=xhot(1);
    Tq_ihot=xhot(2);
    
    %Generation of the output for the hot case
    Trhot=[Trhot,xhot(1)];
    Tqhot=[Tqhot,xhot(2)];
    %---------------------------------------------------------------------------------------------------------------------------------------------------------
    
    %Matrix computation for the cold case of the service module such that it
    %is pointing away from earth when triangle sps-earth-sun is a right triangle
    Acold=[mr*cpr/timestep,0;0,mq*cpq/timestep];
    Bcold=[1,(mr*cpr/timestep)-(k*Aint/L),k*Aint/L,-boltz*Ar*epsilonR,0,0;0,(k*Aint/L),mq*cpq/timestep-k*Aint/L,0,-boltz*Aq*epsilonq,1];
    ucold=[QinR(i);Tr_icold;Tq_icold;Tr_icold^4;Tq_icold^4;QinQ(1,i)];
    ucold_out=[ucold_out,ucold];
    xcold=Acold\(Bcold*ucold);
    
    %update of variables for next step of the cold case
    Tr_icold=xcold(1);
    Tq_icold=xcold(2);
    
    %output storage
    
    %rad=[rad,boltz*Ar*epsilonR*Tr_i^4];
    %stored_ring=[stored_ring,mr*cpr*Tr_i];
    Trcold=[Trcold,xcold(1)];
    Tqcold=[Tqcold,xcold(2)];
    
    
    time=[time,t];

    
    
    %Time keeping and clock
    t=t+timestep;
    i=i+1;
    

    
end

figure
plot(time/60/90,Trcold-273.5,'b');
hold on
plot(time/60/90,Tqcold-273.5,'r')
hold off

