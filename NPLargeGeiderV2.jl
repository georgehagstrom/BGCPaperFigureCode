


function NPLargeGeiderV2(IDay,DayLength,IrrInst,T,po4,nh4,no3,q_fe,params)
    sperd = 86400.0

    kappa_phyto = params["kappa_phyto"];
   
    kappa_photo = params["kappa_photo"];
    P_C_max_lp = params["P_C_max_lp"];
    thetamax_lp = params["thetamax_lp"]
    zeta_lp = params["zeta_lp"]
    bresp_lp = params["bresp_lp"]
    alpha_lp = params["alpha_lp"]
    kno3_lp = params["kno3_lp"]
    knh4_lp = params["knh4_lp"]

    r_lp = params["r_lp"]
    kpo4_lp = params["kpo4_lp"]


	thetaFac = params["thetaFactorLg"]
    q_fe_k_lp = params["q_fe_k_lp"]

    PStorMax = params["PStorMax_lp"] # Maximum level of P storage
    q_fe_lp_max = params["q_fe_lp_max"]
    tfac_lp = exp(kappa_phyto*T);
    bresp = bresp_lp*tfac_lp
    P_C_max = P_C_max_lp*tfac_lp
    k_p_stor_lp = params["k_p_stor_lp"]
    #
    PExtra_lp = params["PExtra_lp"]

    LFac = params["LFac"]

   tfac_photolp = exp(kappa_photo*(T))
    k1 = params["k1"]*tfac_photolp
    k2 = params["k2"]*tfac_photolp
   aph = params["aph_lp"]
 
#aph = params["aph_sm"]
phiM = params["phiM"]



   tfac_lp = exp(kappa_phyto*T);
   bresp = bresp_lp*tfac_lp

   c = cOptDi(IDay*4.57,tfac_photolp)

   PCM_s = min(k1*(1-c),k2*c)/3600.0

   aPhoto = PCM_s*(1-exp(-aph*phiM*c*4.57*IDay/PCM_s))


   kE = P_C_max_lp*tfac_lp



    q_fe_lp_max = 50e-6*(106/16.0)

    alphaE = 0.5    # This parameter is the ribosome fraction of the synthesis pool
    alphaDNA = 0.1 # This is the RNA fraction of the synthesis pool. I don't know what this number really should be, theoretically it could combine with alphaE above
    kPhoto = 7.5e-2/sperd*4.57 # This is the affinity of 1-unit of photosynthetic apparatus for light
    alphaS = .2;     # Parameter relating the size of the periplasmic lpace to the radius. Fraction of cell that is membrane+periplasm is alphaS/r. 0.35 is a high value
    gammaDNA = 0.0005;   # Fraction of cell that is DNA+non-ribosomal RNA
    gamma = .2;     # Fraction of cell devoted to "structure" modeled here as the gammaDNA above plus a lipid pool. The composition of this pool is a bit arbitrary


    fLipid = 0.25;   # Fraction of inner and outermembrane composed of lipids

    fProt=1-fLipid;      # Fraction of inner and outer membrane composed of proteins
    molarN = 14.0;   # mass of Nitrogen in g/mol

    gammaLipid = gamma-gammaDNA;  # Fraction of structure pool that is composed of ordinary lipids


    # These are the percent of dry mass that is Pholphorus, Nitrogen, or Carbon of different macromolecules 

    PDNA = 0.095;    # P in DNA
    PRib = 0.047;   # P in a ribosome (roughly 50# RNA 50# protein)
    NProt = .16;    # N in Protein
    NDNA = .16;     # N in DNA
    PPhospholipid = 0.042;  # P in pholpholipids
    #PE=PRib*alphaE+PDNA*alphaDNA; # The biosynthetic apparatus is alphaE ribosomes and (1-alphaE) proteins. PE is the P content of the biosynthesis apparatus
    PE = 0.05
    S_lp = alphaS/r_lp+gamma # This constant is just how big the structure pool is
    
    #print(S_lp)
    #NS_lp = (gammaDNA*NProt+fProt*NProt*alphaS/r_lp)/S_lp  # This is the N content of the structure pool overall, in gN/g dry mass
    NS_lp = NProt * 0.8

    molarP=31.0; # Molar mass of P


    # calculate the nutrient limitation (Frost and Franzen) Standard Cobalt all
    no3lim_lp = no3/( (kno3_lp + no3)*(1 + nh4/knh4_lp) );


    nh4lim_lp = nh4/(knh4_lp + nh4);
    nlim_lp = no3lim_lp+nh4lim_lp;

    plim_lp = po4/(kpo4_lp+po4)

    plim_stor_lp = po4/(k_p_stor_lp+po4)
    felim_lp = q_fe^2/(q_fe^2+q_fe_k_lp^2)



    
    nutlim_lp = min(nlim_lp,plim_lp,felim_lp)

    # Here we insist the cell does not investment any more in E that would cause it to grow faster than nutrient limitation would imply
    
    

    ENutLim_lp = nutlim_lp
    
    #LNutLim_lp = -P_C_max/(kPhoto*IDay)*log(1 - (ENutLim_lp*(1+zeta_lp) + bresp/P_C_max)/D   )
    LNutLim_lp = (kE*ENutLim_lp*(1+zeta_lp)+bresp)/(DayLength*aPhoto/LFac)


    # Next estimate the alphaPhoto coefficient


    #LNutLim = -P_C_max/(kPhoto*IDay)*log(1 - (ENutLim_sp*(1+zeta_sp) + bresp/P_C_max)/D   )
    

    
    
    
          
                   
    # Check if the nutrient limited solution overinvests
    if (LNutLim_lp+ENutLim_lp>1-S_lp) ||  (LNutLim_lp<0)


        # If yes, find the max growth rate solution
        # We define an RHS function
	#EMaxGrowth = (1-S_lp)/(1 + (kE*(1+zeta_lp)+bresp)/(DayLength*aPhoto))
	

	EMaxGrowth = ((1-S_lp)*DayLength*aPhoto/LFac-bresp)/((kE*(1+zeta_lp))+(DayLength*aPhoto/LFac))
	LMaxGrowth =  (kE*EMaxGrowth*(1+zeta_lp)+bresp)/(DayLength*aPhoto/LFac)

	if EMaxGrowth<0
		EMaxGrowth = 0
		LMaxGrowth = 1-S_lp
	end



       		
	EOpt_lp = EMaxGrowth
	LOpt_lp = LMaxGrowth
        nutlimFlag_lp = 0

		
	else
	# If not, lock in the nutrient limited solution
        	EOpt_lp = ENutLim_lp
	        LOpt_lp = LNutLim_lp
        	nutlimFlag_lp = argmin([nlim_lp,plim_lp,felim_lp])
    end
       PStor_lp = plim_stor_lp*PStorMax # Level of P Storage

    	theta_lp = LOpt_lp * kPhoto/(alpha_lp*(S_lp+EOpt_lp+LOpt_lp)) # Chlorophyl to Carbon

    # Calculate NP_lp and structural NP. Note the abscence of storage term in the latter. Probably can disregard NP_struct_lp for Cobalt version unless you want to 
    # explicitly track storage

    NP_lp = ((S_lp*NS_lp+NProt*(EOpt_lp+LOpt_lp))/(PE*EOpt_lp+PExtra_lp+PStor_lp*(EOpt_lp+LOpt_lp+S_lp) ))*molarP/molarN

    NP_struct_lp = ((S_lp*NS_lp+NProt*(EOpt_lp+LOpt_lp))/(PE*EOpt_lp+PExtra_lp ))*molarP/molarN
    
    #print("NP_lp ",NP_lp," NP_struct_lp ",NP_struct_lp," EOpt_lp ", EOpt_lp , " LOpt_lp ", LOpt_lp,"\n")
    # Return
    IrrLimInst = (PCM_s*LOpt_lp*(1-exp(-aph*phiM*c*4.57*IrrInst/PCM_s))-bresp)/(1+zeta_lp)

    theta = thetaFac*LOpt_lp*c/(LOpt_lp+EOpt_lp+S_lp)
    lpSol = @SVector [EOpt_lp,LOpt_lp,EOpt_lp*kE*sperd,NP_lp,NP_struct_lp,nutlimFlag_lp,theta,IrrLimInst,LOpt_lp*c,nlim_lp]
end

