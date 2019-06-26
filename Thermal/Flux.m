function [RingFlux]=Flux(numelem)
    

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
    
    
    %%Ring Definition
    Aring1=0.02975;%area that is exposed if the z-axis is directly aligned with the incoming radiation(8)
    alpha1_ring=0.15;
    epsilon1_ring=0.05;

    Aring2=0.28935%area that is exposed if the x&yaxis are aligned with the incoming radiation(half cylinder)
    alpha2_ring=0.15;
    epsilon2_ring=0.05;

    Aring3=0.02975%area that is exposed if the x&yaxis are aligned with the incoming radiation(half cylinder)
    alpha3_ring=0.15;
    epsilon3_ring=0.05;
    
    
    %%
    %%
    %%HeatFlux Input to the quadpack per orbit calculation
    
    
    circ1=linspace(0,pi/2,numelem/6);

    qr_1=qs*Aring1*alpha1_ring+qIR*Redmax*epsilon3_ring*Aring3*cos(circ1)+...
        qIR*Redmax*epsilon2_ring*Aring2*cos(pi/2-circ1)+qs*albedo*alpha3_ring*Aring3*cos(circ1)+...
        qs*albedo*alpha2_ring*Aring2*cos(pi/2-circ1);


    %--------------------------------------------------------------------------------------------------------
    circ2=linspace(pi/2,circ_ecli,numelem/6);

    qr_2=qs*Aring1*alpha1_ring+qIR*Redmax*epsilon1_ring*Aring1*(-cos(circ2))+...
        qIR*Redmax*epsilon2_ring*Aring2*cos(circ2-pi/2)+qs*albedo*alpha1_ring*Aring1*(-cos(circ2))+...
        qs*albedo*alpha1_ring*Aring2*cos(circ2-pi/2);

    %-----------------------------------------------------------------------------------------------------
    circ3=linspace(circ_ecli,pi,numelem/6);


    qrBest_3=qIR*Redmax*epsilon1_ring*Aring1*(-cos(circ2))+...
        qIR*Redmax*epsilon2_ring*Aring2*cos(circ2-pi/2);

    %Now just mirror the plot to get the worst case, in reality it could be
    %that the SPS2 might rotate around the z-axis such that quadpack1 will have
    %the orientation to earth which quadpack4 had in the first part of the
    %orbit, meaning that the cold case will now be quadpack 4

    circ4=linspace(pi,2*pi,numelem/2);
    
    qr_half=[qr_1,qr_2,qrBest_3];
    qr=[qr_half,flip(qr_half)];
    circ=[circ1,circ2,circ3,circ4];
    
    
    
    QuadFluxHot=q4;
    QuadFluxCold=q1;
    RingFlux=qr;
end







