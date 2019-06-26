%Four Node Thermal Analysis
clear all
numorbits=6;
Torb=90*60;
timestep=1;
sensitivityTr=[];
sensitivityTSVM=[];
sensitivityTQP1=[];
sensitivityQPo=[];
senseTime=[];

%%Makrolon GV 30 Characteristics(Interface from Ring to Service Module)
%%
coun=0;
hQR=500;%Thermal conductivity over an alu interface in vacuum


kSR=0.16;%Thermal Conductivity W/mKZ
LSR=0.01;%thickness of the glass fiber insulation
ASR=0.051076;%Interface area of Ring and Service Module 226mm X 226mm

%%Alu Interface(Interface from Ring to QuadPack)
%%

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
cpSVM=50; %specific heat of the service module set to 200J/kgK
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



PinR=RingFlux((Torb/timestep),numorbits);%Input raidation heat to the ring, connect to Flux function later

PinQ=QuadFlux((Torb/timestep),numorbits);%Input radiation to the cubesat


boltz=5.67*10^-8;%boltzmann constant

while coun<10
    
    
    Tr_i=273.5;%Initial temperature of the ring
    TSVM_i=273.5;%Initial temeprature of the Service Module
    TQP1_i=273.5;%Initial temeprature of the QuapdPack opposite to SVM
    TQPo_i=273.4;%Initial temeprature of other Quadpacks


    %%
    %%Propagator
    time=[];%initialize time output array
    Tr=[];%initialize Ring Temperature output array
    TSVM=[];%initialize Service Module Temperature output array
    TQP1=[];%initialize QuadPack opposite to Service Module Temperature output array
    TQPo=[];%initialize other QuadPack Temperature output array

    P_rad=[];

    i=1;
    t=0;
    diff=[];

    
    while t<Torb*numorbits-1

        %Matrix computation for the worst case of the service module such that it
        %is pointing at earth when triangle sps-earth-sun is a right triangle
        %Note that for the first matrix the components for the quadpack
        %opposite of the SVM and the ones for the other quadpacks is identical.
        %The distinction arises from the power input
        A=[mr*cpr/timestep,0,0,0;0,mSVM*cpSVM/timestep,0,0;0,0,mQP*cpQP/timestep,0;0,0,0,mQP*cpQP/timestep];



        B=[(mr*cpr/timestep)-(kSR*ASR/LSR)-(hQR*AQR)-4*(hQR*AQR),-boltz*Ar*epsilonR,(kSR*ASR/LSR),0,(hQR*AQR),0,4*(hQR*AQR),0,1,0,0,0;...
            (kSR*ASR/LSR),0,(mSVM*cpSVM/timestep)-(kSR*ASR/LSR),-boltz*(Area1*epsilon1SVM+Area2*epsilon2SVM+Area3*epsilon3+Area4*epsilon4SVM+Area5*epsilon5SVM+Area6*epsilon6SVM),0,0,0,0,0,1,0,0;...
            (hQR*AQR),0,0,0,(mQP*cpQP/timestep)-(hQR*AQR),-boltz*(Area1*epsilon1QP+Area2*epsilon2QP+Area3*epsilon3+Area4*epsilon4QP+Area5*epsilon5QP+Area6*epsilon6QP),0,0,0,0,1,0;...
            (hQR*AQR),0,0,0,0,0,(mQP*cpQP/timestep)-(hQR*AQR),-boltz*(Area1*epsilon1QP+Area2*epsilon2QP+Area3*epsilon3+Area4*epsilon4QP+Area5*epsilon5QP+Area6*epsilon6QP),0,0,0,1];

        u=[Tr_i;Tr_i^4;TSVM_i;TSVM_i^4;TQP1_i;TQP1_i^4;TQPo_i;TQPo_i^4;PinR(i);PinQ(1,i);PinQ(2,i);PinQ(3,i)];

        x=A\(B*u);

        %update of variables for next step 
        Tr_i=x(1);

        TSVM_i=x(2);

        TQP1_i=x(3);

        TQPo_i=x(4);

        %Generation of the output 
        Tr=[Tr,x(1)];
        TSVM=[TSVM,x(2)];
        TQP1=[TQP1,x(3)];
        TQPo=[TQPo,x(4)];


        P_rad=[P_rad,boltz*(Area1*epsilon1QP+Area2*epsilon2QP+Area3*epsilon3+Area4*epsilon4QP+Area5*epsilon5QP+Area6*epsilon6QP)*TQPo_i^4];

        time=[time,t];

        %Time keeping and clock
        t=t+timestep;
        i=i+1
        
    end
%     sensitivityTr=[sensitivityTr;Tr];
%     sensitivityTSVM=[sensitivityTSVM;TSVM];
%     sensitivityTQP1=[sensitivityTQP1;TQP1];
%     sensitivityQPo=[sensitivityQPo;TQPo];
%     senseTime=[senseTime;time];
    timestep=timestep/2
    coun=coun+2;
end

% pp=1;
% figure
% while pp<6
%     plot(senseTime(pp,:),sensitivityTr(pp,:)-273.5)
%     hold on
%     pp=pp+1;
% end
% title('Mounting Ring Node Temperature Variation for Different Conductivity Values')
% xlabel('Number of Orbits')
% ylabel('Temperature °C')
% grid on
% hold off
% 
% pp=1;
% figure
% while pp<6
%     plot(senseTime(pp,:),sensitivityTSVM(pp,:)-273.5)
%     hold on
%     pp=pp+1;
% end
% title('Service Module Node Temperature Variation for Different Conductivity Values')
% xlabel('Number of Orbits')
% ylabel('Temperature °C')
% %legend('Conductivity of 1W/m^2K','Conductivity of 101W/m^2K','Conductivity of 201W/m^2K','Conductivity of 301W/m^2K','Conductivity of 401W/m^2K','Conductivity of 501W/m^2K')
% grid on
% hold off
% 
% pp=1;
% figure
% while pp<6
%     plot(senseTime(pp,:),sensitivityTQP1(pp,:)-273.5)
%     hold on
%     pp=pp+1;
% end
% title('QP1')
% grid on
% hold off
% 
% pp=1;
% figure
% while pp<6
%     plot(senseTime(pp,:),sensitivityQPo(pp,:)-273.5)
%     hold on
%     pp=pp+1;
% end
% title('QP Other')
% grid on
% figure
% plot(time/60/90,Tr-273.5,'b');
% hold on
% plot(time/60/90,TSVM-273.5,'r')
% hold on
% plot(time/60/90,TQP1-273.5,'k')
% hold on
% plot(time/60/90,TQPo-273.5,'y')
% title('Temperature Variation of Nodes for the max Power input to the Service Module')
% xlabel('Number of orbits')
% ylabel('Temperature[°C]')
% grid on
% legend(' Mounting Ring','Service Module Node','QuadPack 1 Node','Other Quadpacks')
% hold off
