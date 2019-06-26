function p=RingFlux(numelem,orbnum)
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
    adcsError=0.5;
    

    %%
    
    
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
    
    
    %%
    %%
    %%HeatFlux Input to the quadpack per orbit calculation
    %qs*Aring3*alpha2_ring*cos(deg2rad(adcsError))
   
    circ1=linspace(0,pi/2,numelem/6);

    Qring1_hot=qs*Aring1*alpha1_ring*cos(deg2rad(adcsError))+qs*Aring2*alpha2_ring*cos(deg2rad(90-adcsError))+qs*Aring3*alpha3_ring*cos(deg2rad(90-adcsError))...
        +qIR*Redmax*epsilon3_ring*Aring3*cos(circ1)+...
        qIR*Redmax*epsilon2_ring*Aring2*cos(pi/2-circ1)+qs*albedo*alpha3_ring*Aring3*cos(circ1)+...
        qs*albedo*alpha2_ring*Aring2*cos(pi/2-circ1)
    
    Qring1_cold=qs*Aring1*alpha1_ring*cos(deg2rad(adcsError))+...
        +qIR*Redmax*epsilon3_ring*Aring3*cos(circ1)+...
        qIR*Redmax*epsilon2_ring*Aring2*cos(pi/2-circ1)+qs*albedo*alpha3_ring*Aring3*cos(circ1)+...
        qs*albedo*alpha2_ring*Aring2*cos(pi/2-circ1);

    %--------------------------------------------------------------------------------------------------------
    circ2=linspace(pi/2,circ_ecli,numelem/6);

    Qring2_hot=qs*Aring1*alpha1_ring*cos(deg2rad(adcsError))+qs*Aring2*alpha2_ring*cos(deg2rad(90-adcsError))+qs*Aring3*alpha3_ring*cos(deg2rad(90-adcsError))...
        +qIR*Redmax*epsilon1_ring*Aring1*(-cos(circ2))+...
        qIR*Redmax*epsilon2_ring*Aring2*cos(circ2-pi/2)+qs*albedo*alpha1_ring*Aring1*(-cos(circ2))+...
        qs*albedo*alpha1_ring*Aring2*cos(circ2-pi/2);

        
    Qring2_cold=qs*Aring1*alpha1_ring*cos(deg2rad(adcsError))+...
        +qIR*Redmax*epsilon1_ring*Aring1*(-cos(circ2))+...
        qIR*Redmax*epsilon2_ring*Aring2*cos(circ2-pi/2)+qs*albedo*alpha1_ring*Aring1*(-cos(circ2))+...
        qs*albedo*alpha1_ring*Aring2*cos(circ2-pi/2);
    %-----------------------------------------------------------------------------------------------------
    circ3=linspace(circ_ecli,pi,numelem/6);
    cst1=ones(numelem/6,1);
    
    Qring3_hot= qIR*Redmin*epsilon3_ring*Aring3*1*cst1;
    Qring3_cold= qIR*Redmax*epsilon1_ring*Aring1*cst1;
    
    circ4=linspace(pi,2*pi,numelem/2);
    circ=[circ1,circ2,circ3,circ4];
    %Now just mirror the plot to get the worst case, in reality it could be
    %that the SPS2 might rotate around the z-axis such that quadpack1 will have
    %the orientation to earth which quadpack4 had in the first part of the
    %orbit, meaning that the cold case will now be quadpack 4
    
%     %1)HotSunColdEclipse
%     RingHotCold_half=[Qring1_hot,Qring2_hot,Qring3_cold'];
%     RingHotCold=[RingHotCold_half,flip(RingHotCold_half)];
%     figure
%     plot(rad2deg(circ),RingHotCold,'b')
%     grid on
%     title('Ring Worst Case Power Input , Hot Case Sunlight,Cold Case Eclipse ')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 0 50])
%     %2)HotSunHotEclipse
%     RingHotHot_half=[Qring1_hot,Qring2_hot,Qring3_hot'];
%     RingHotHot=[RingHotHot_half,flip(RingHotHot_half)];
%     
%     figure
%     plot(rad2deg(circ),RingHotHot,'k')
%     grid on
%     title('Ring Worst Case Power Input , Hot Case Sunlight,Hot Case Eclipse ')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 0 50])
%     
%     %3)ColdSunHotEclipse
%     RingColdHot_half=[Qring1_cold,Qring2_cold,Qring3_hot'];
%     RingColdHot=[RingColdHot_half,flip(RingColdHot_half)];
% 
%     
%     figure
%     plot(rad2deg(circ),RingColdHot,'r')
%     grid on
%     title('Ring Worst Case Power Input , Cold Case Sunlight,Hot Case Eclipse ')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 0 50])
%     
%     %4)ColdSunColdEclipse
%     RingColdCold_half=[Qring1_cold,Qring2_cold,Qring3_cold'];
%     RingColdCold=[RingColdCold_half,flip(RingColdCold_half)];
%     
%     figure
%     plot(rad2deg(circ),RingColdCold,'c')
%     grid on
%     title('Ring Worst Case Power Input , Cold Case Sunlight,Cold Case Eclipse ')
%     xlabel('Orbit Position given by \beta')
%     ylabel('Total Power Input[W]')
%     axis([0 360 0 50])
%     
%     
    
    
    %Current assumption is that the worst case is hot in sunglight and cold
    %in eclipse.
    QRworst_half=[Qring1_cold,Qring2_cold,Qring3_hot'];
    QRworst=[QRworst_half,flip(QRworst_half)];
   
    
    
    j=1;
    fluxR=[];
    Position=[];
    while j<=orbnum
        
        fluxR=[fluxR,QRworst];
        Position=[Position,circ];
        j=j+1;
        
    end
    

    
    p=fluxR;

end







