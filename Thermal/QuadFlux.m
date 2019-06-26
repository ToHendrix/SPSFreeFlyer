%Flux Quadpacks

 function f=QuadFlux(numelem,orbnum)
    %%QuadPackDefinition   
    
    QintSun=81;
    QintEclipse=43;
    heaterpow=0;
    %Heat Flux input into the quadpack for a circular 350km orbit
    Quadwidth=272;
    Quadlength=391;
    QuadHeight=272;
    
    
    
    
    
    %The LABELLLLS ARE WRONG USE THE LABELS IN THE EXCEL!!!PLEASE
    
    
    
    
    
    %---------------------------------------------------------------------------------------------------
    Area1=272*272/(1000*1000);%Small Panel ontop(Service Module =Radiator Area,Quadpack=Black Anodised Alu)
    alpha1SVM=0.248;
    epsilon1SVM=0.888;
    
    alpha1QP=0.65;
    epsilon1QP=0.82;
    %---------------------------------------------------------------------------------------------------
    Area2=272*272/(1000*1000);%Bottom Panel(Service Module =insulated,Quadpack=insulated)
    alpha2SVM=0;
    epsilon2SVM=0;
    
    alpha2QP=0;
    epsilon2QP=0;
    %---------------------------------------------------------------------------------------------------
    Area3=272*391/(1000*1000);%Solar Panel Covered(Service Module =Solar Panel,Quadpack=Solar Panel)
    alpha3=0.805;%From old SMAD pg 422
    epsilon3=0.825;%From old SMAD pg 422
    eff_sp=0.33;%panel efficiency
    packing_sp=0.8;%packing density of solar panels
    Kdiss=1;%Electronics Heat Dissipation factor

    %---------------------------------------------------------------------------------------------------
    Area4=272*391/(1000*1000);%Side Panel (Service Module =Kapton Film,Quadpack=Blue Anodised Alu)
    alpha4SVM=0.49;
    epsilon4SVM=0.3;

    alpha4QP=0.53;
    epsilon4QP=0.82;
    %---------------------------------------------------------------------------------------------------
    Area5=272*391/(1000*1000);%Side Panel (Service Module =Highly Radiative Surface,Quadpack=Blue Anodised Alu)
    alpha5SVM=0.248;
    epsilon5SVM=0.888;

    alpha5QP=0.53;
    epsilon5QP=0.82;
    %---------------------------------------------------------------------------------------------------

    Area6=272*391/(1000*1000);%Side Panel (Service Module =Kapton Film,Quadpack=Blue Anodised Alu)
    alpha6SVM=0.49;
    epsilon6SVM=0.3;

    alpha6QP=0.53;
    epsilon6QP=0.82;
    
    %Constants
    qs=1365;
    qIR=258;
    albedo=0.35;
    boltz=5.67*10^-8;%boltzmann cst
    Re=6371000;
    H=350*1000;
    circ_ecli=(pi/2)+atan(Re/(Re+H));
    Redmax=((Re^2))/((Re+800000)^2);
    Redmin=((Re^2))/((Re+350000)^2);
    Torb=91;
    
    %%
    %%HeatFlux Input to the quadpack per orbit calculation
    circ1=linspace(0,pi/2,numelem/6);

    SVM1_1=(alpha3-(eff_sp*packing_sp))*Area3*qs+qIR*Redmax*epsilon5SVM*Area5*cos(circ1)+...
        qs*albedo*alpha5SVM*Area5*cos(circ1)+QintSun;
    SVM4_1=(alpha3-(eff_sp*packing_sp))*Area3*qs+qIR*Redmin*epsilon5SVM*Area5*cos(circ1)+...
        qIR*epsilon1SVM*Redmin*Area1*cos((pi/2)-circ1)+qs*albedo*alpha5SVM*Area5*cos(circ1)+...
        qs*albedo*alpha1SVM*Area1*cos((pi/2)-circ1)+QintSun;
    
    
    QuadA_1=(alpha3-(eff_sp*packing_sp))*Area3*qs+qIR*Redmax*epsilon5QP*Area5*cos(circ1)+...
        qs*albedo*alpha5QP*Area5*cos(circ1);
    
    QuadB_1=(alpha3-(eff_sp*packing_sp))*Area3*qs+qIR*Redmin*epsilon5QP*Area5*cos(circ1)+...
        qIR*epsilon1QP*Redmin*Area1*cos((pi/2)-circ1)+qs*albedo*alpha5QP*Area5*cos(circ1)+...
        qs*albedo*alpha1QP*Area1*cos((pi/2)-circ1);
    

    
    
    
    %--------------------------------------------------------------------------------------------------------
    circ2=linspace(pi/2,circ_ecli,numelem/6);

    SVM4_2=(alpha3-(eff_sp*packing_sp))*Area3*qs+qIR*Redmin*epsilon1SVM*Area1*cos((pi/2)-circ2)+...
        +qIR*Redmin*epsilon3*Area3*(-cos(circ2))+qs*albedo*alpha1SVM*Area1*cos((pi/2)-circ2)+...
        +qs*albedo*alpha1SVM*Area3*(-cos(circ2))+QintSun;
    SVM1_2=(alpha3-(eff_sp*packing_sp))*Area3*qs+qIR*Redmax*epsilon3*Area3*(-cos(circ2))+...
        qs*albedo*alpha1SVM*Area3*(-cos(circ2))+QintSun;
    
    
    
    QuadB_2=(alpha3-(eff_sp*packing_sp))*Area3*qs+qIR*Redmin*epsilon1QP*Area1*cos((pi/2)-circ2)+...
        +qIR*Redmin*epsilon3*Area3*(-cos(circ2))+qs*albedo*alpha1QP*Area1*cos((pi/2)-circ2)+...
        +qs*albedo*alpha1QP*Area3*(-cos(circ2));
    
    QuadA_2=(alpha3-(eff_sp*packing_sp))*Area3*qs+qIR*Redmax*epsilon3*Area3*(-cos(circ2))+...
        qs*albedo*alpha1QP*Area3*(-cos(circ2));



    %-----------------------------------------------------------------------------------------------------
    circ3=linspace(circ_ecli,pi,numelem/6);
    cst=ones(numelem/6,1);
    %q4best_3=qIR*Redmin*epsilon1*Area1*cos((pi/2)-circ3)+qIR*Redmin*epsilon3*Area3*(-cos(circ3))+QintEclipse;
    SVM4_3=qIR*Area4*epsilon4SVM*cst+QintEclipse;
    SVM1_3=QintEclipse*cst;
    
    QuadB_3=qIR*Area4*epsilon4QP*cst;
    QuadA_3=cst*0;

    
    

    

    %Now just mirror the plot to get the worst case, in reality it could be
    %that the SPS2 might rotate around the z-axis such that quadpack1 will have
    %the orientation to earth which quadpack4 had in the first part of the
    %orbit, meaning that the cold case will now be quadpack 4

    circ4=linspace(pi,2*pi,numelem/2);
    circ=[circ1,circ2,circ3,circ4];

%---------------------------------------------------------------------------------------------------------------------------------------------------------

   %There are basically four combinations of hot and cold cases that have
    %to be considered for the quadpacks
    
    
%     
%     %1)HotSunHotEclipse
%     QPHotHot_half=[QuadB_1,QuadB_2,QuadB_3'];
%     QPHotHot=[QPHotHot_half,flip(QPHotHot_half)];
%     figure
%     plot(rad2deg(circ),QPHotHot,'b')
%     grid on
%     title('Worst Case Power Input , Hot Case Sunlight,Cold Case Eclipse')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 40 210])
%     [0.4 0.4 0.4]
%    
%     
%     
%     %4)ColdSunColdEclipse
%     QPColdCold_half=[QuadA_1,QuadA_2,QuadA_3'];
%     QPColdCold=[QPColdCold_half,flip(QPColdCold_half)];
%     
%     figure
%     plot(rad2deg(circ),QPColdCold,'c')
%     grid on
%     title('Worst Case Power Input , Cold Case Sunlight,Cold Case Eclipse')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 40 210])
    
    
%     %Quadpack Power Input
%     QPHotHot_half=[QuadB_1,QuadB_2,QuadB_3'];
%     QPHotHot=[QPHotHot_half,flip(QPHotHot_half)];
%     
%     QPColdCold_half=[QuadA_1,QuadA_2,QuadA_3'];
% 	QPColdCold=[QPColdCold_half,flip(QPColdCold_half)];
%     
%     figure
%     plot(rad2deg(circ),(QPHotHot+QPColdCold)/2,'k')
%     hold on
%     
%     plot(rad2deg(circ),QPColdCold,'b')
%     hold on
% 
% 
%     plot(rad2deg(circ),QPHotHot,'r')
%     grid on
%     title('Average Worst Case Power Input , Hot Case Sunlight,Cold Case Eclipse')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 0 150])
%     legend('Power Input Quadpack 2,3,4,5 Averaged','Power Input Quadpack 1 Cold Case','Power Input Quadpack 1 Hot Case');
%     hold off
    
    
%---------------------------------------------------------------------------------
    


%------------------------------------------------------------------------------------------------------------------------------------------------------------

    
    
    
    
    %There are basically four combinations of hot and cold cases that have
    %to be considered for the service module 
    
%     %1)HotSunColdEclipse
%     SVMHotCold_half=[SVM4_1,SVM4_2,SVM1_3'];
%     SVMHotCold=[SVMHotCold_half,flip(SVMHotCold_half)];
%     figure
%     plot(rad2deg(circ),SVMHotCold,'b')
%     grid on
%     title('Worst Case Power Input , Hot Case Sunlight,Cold Case Eclipse')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 40 210])
%     [0.4 0.4 0.4]
%     %2)HotSunHotEclipse
%     SVMHotHot_half=[SVM4_1,SVM4_2,SVM4_3'];
%     SVMHotHot=[SVMHotHot_half,flip(SVMHotHot_half)];
%     
%     figure
%     plot(rad2deg(circ),SVMHotHot,'k')
%     grid on
%     title('Worst Case Power Input , Hot Case Sunlight,Hot Case Eclipse')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 40 210])
%     
%     %3)ColdSunHotEclipse
%     SVMColdHot_half=[SVM1_1,SVM1_2,SVM4_3'];
%     SVMColdHot=[SVMColdHot_half,flip(SVMColdHot_half)];
%     axis([0 360 40 210])
%     
%     figure
%     plot(rad2deg(circ),SVMColdHot,'r')
%     grid on
%     title('Worst Case Power Input , Cold Case Sunlight,Hot Case Eclipse')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 40 210])
%     
%     %4)ColdSunColdEclipse
%     SVMColdCold_half=[SVM1_1,SVM1_2,SVM1_3'];
%     SVMColdCold=[SVMColdCold_half,flip(SVMColdCold_half)];
%     
%     figure
%     plot(rad2deg(circ),SVMColdCold,'c')
%     grid on
%     title('Worst Case Power Input , Cold Case Sunlight,Cold Case Eclipse')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 40 210])
    
    

    
    
    %---------------------------------------------------------------------------------
    
    
    SVMworst_half=[SVM4_1,SVM4_2,SVM4_3']
    SVMworst=[SVMworst_half,flip(SVMworst_half)];
    
    Quad_half1=[QuadA_1,QuadA_2,QuadB_3'];
    Quad1=[Quad_half1,flip(Quad_half1)];
    
    Quad_halfavg=[(QuadA_1+QuadB_1)/2,(QuadA_2+QuadB_2)/2,(QuadA_3'+QuadB_3')/2];
    Quado=[Quad_halfavg,flip(Quad_halfavg)];
    
    
    o=1;
    FluxSVMworst=[];
    FluxQuad1=[];
    FluxQuado=[];
    while o<=orbnum
        FluxQuado=[FluxQuado,Quado];
        FluxSVMworst=[FluxSVMworst,SVMworst];
        FluxQuad1=[FluxQuad1,Quad1];
        o=o+1;
        
    end

    f=[FluxSVMworst;FluxQuad1;FluxQuado];
 end