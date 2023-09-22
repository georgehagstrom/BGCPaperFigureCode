function NPSmallGeiderV2(IDay,DayLength,IrrInst,T,po4,nh4,no3,q_fe,params)
sperd = 86400.0

kappa_phyto = params["kappa_phyto"];
kappa_photo = params["kappa_photo"]
alpha_sp = params["alpha_sp"];
P_C_max_sp = params["P_C_max_sp"];
zeta_sp = params["zeta_sp"]
bresp_sp = params["bresp_sp"]
r_sp = params["r_sp"]
kno3_sp = params["kno3_sp"]
knh4_sp = params["knh4_sp"]


kpo4_sp = params["kpo4_sp"]


q_fe_k_sp = params["q_fe_k_sp"]

k_fed_sp = params["k_fed_sp"]

k_p_stor_sp = params["k_p_stor_sp"]
q_fe_sp_max = params["q_fe_sp_max"]


tfac_photo_sp = exp(kappa_photo*(T-15.0))
k1 = params["k1"]*tfac_photo_sp
k2 = params["k2"]*tfac_photo_sp
aph = params["aph_sm"]
phiM = params["phiM"]
LFac = params["LFac"]
thetaFac = params["thetaFactorSm"]

tfac_sp = exp(kappa_phyto*T);
bresp = bresp_sp*tfac_sp

c = cOpt(4.57*IDay,tfac_photo_sp)

PCM_s = min(k1*(1-c),k2*c)/3600.0

aPhoto = PCM_s*(1-exp(-aph*phiM*c*4.57*IDay/PCM_s))
kE = P_C_max_sp*tfac_sp


PStorMax = params["PStorMax_sp"]  # Maximum level of P storage



    q_fe_sp_max = 50e-6*(106/16.0)

    PExtra_sp = params["PExtra_sp"]
    alphaE = 0.5  # ribosome content of the synthesis pool
    alphaDNA = 0.1 # RNA content (ex-ribosomes) of the synthesis pool, sort of a fudge factor that should be merged with alphaE 
    kPhoto = 7.5e-2/sperd*4.57  # Affinity of a unit of photosynthesis investment
    alphaS = .20;     # Parameter relating the size of the periplasmic space to the radius. Fraction of cell that is membrane+periplasm is alphaS/r. 0.35 is a high value
    gammaDNA = 0.0035;   # Fraction of cell that is DNA+non-ribosomal RNA
    gamma = .2;     # Fraction of cell devoted to "structure" modeled here as the gammaDNA above plus a lipid pool. The composition of this pool is a bit arbitrary
    kST0=0.185;       # Efficiecny of the biosynthetic apparatus at 25 celsius. This also impacts the strength of the growth rate hypothesis. Might not need two parameters for this.


    fLipid = 0.25;   # Fraction of inner and outermembrane composed of lipids

    fProt=1-fLipid;      # Fraction of inner and outer membrane composed of proteins
    molarN = 14.0;   # mass of Nitrogen in g/mol


    gammaLipid = gamma-gammaDNA;  # Fraction of structure pool that is composed of lipids


        # These are the percent of dry mass that is Phosphorus, Nitrogen, or Carbon of different macromolecules 

    PDNA = 0.095;    # P in DNA
    PRib = 0.047;   # P in a ribosome (roughly 50# RNA 50# protein)
    NProt = .16;    # N in Protein
    NDNA = .16;     # N in DNA


    PPhospholipid = 0.042;  # P in phospholipids

        # Here we calculate the temperature impact on diffusivity of inorganic P and N. 
        # A Q10 of 1.5 is not a bad approximation, better ones look much more complex.

        # Here we invoke the alpha function we defined earlier. This is the efficiency of photosynthetic apparatus at irradiance I. We assume it is temperature dependent



#    PE=PRib*alphaE+PDNA*alphaDNA; # The biosynthetic apparatus is alphaE ribosomes and (1-alphaE) proteins. PE is the P content of the biosynthesis apparatus
    PE = 0.05
    S_sp = alphaS/r_sp+gamma  # Size of structure pool
    
    PS_stor_sp = alphaS/r_sp*fLipid*PPhospholipid # Maximum level of P storage due to phospholipids
#    NS_sp = (gammaDNA*NProt+fProt*NProt*alphaS/r_sp)/S_sp  # N content of structure pool
    NS_sp = NProt * 0.8

    molarP=31.0; # Molar mass of P


    # calculate the nutrient limitation (Frost and Franzen)
    no3lim_sp = no3/( kno3_sp + no3 + nh4*kno3_sp/knh4_sp );


    nh4lim_sp = nh4/(knh4_sp + nh4 + no3*knh4_sp/kno3_sp  );
    nlim_sp = no3lim_sp+nh4lim_sp;
    plim_stor_sp = po4/(k_p_stor_sp+po4)
    plim_sp = po4/(kpo4_sp+po4)
    PStor_sp = plim_stor_sp*PStorMax
     

    felim_sp = q_fe^2/(q_fe^2+q_fe_k_sp^2)

    nutlim_sp = min(nlim_sp,plim_sp,felim_sp)

    
    plim_stor_sp = po4/(k_p_stor_sp+po4)
    
    ENutLim_sp = nutlim_sp
    LNutLim_sp = (kE*ENutLim_sp*(1+zeta_sp)+bresp)/(DayLength*aPhoto/LFac)


    # Next estimate the alphaPhoto coefficient


    #LNutLim = -P_C_max/(kPhoto*IDay)*log(1 - (ENutLim_sp*(1+zeta_sp) + bresp/P_C_max)/D   )
    

    
    
    
          
                   
    # Check if the nutrient limited solution overinvests
    if (LNutLim_sp+ENutLim_sp>1-S_sp) ||  (LNutLim_sp<0)


        # If yes, find the max growth rate solution
        # We define an RHS function
	#EMaxGrowth = (1-S_sp)/(1 + (kE*(1+zeta_sp)+bresp)/(DayLength*aPhoto))
	
	EMaxGrowth = ((1-S_sp)*DayLength*aPhoto/LFac-bresp)/((kE*(1+zeta_sp))+(DayLength*aPhoto/LFac))
	LMaxGrowth =  (kE*EMaxGrowth*(1+zeta_sp)+bresp)/(DayLength*aPhoto/LFac)


	if EMaxGrowth<0
		EMaxGrowth = 0
		LMaxGrowth = 1-S_sp
	end


       		
	EOpt_sp = EMaxGrowth
	LOpt_sp = LMaxGrowth
        nutlimFlag_sp = 0
		
	else
	# If not, lock in the nutrient limited solution
        	EOpt_sp = ENutLim_sp
	        LOpt_sp = LNutLim_sp
        	nutlimFlag_sp = argmin([nlim_sp,plim_sp,felim_sp])
    end
    
       PStor_sp = plim_stor_sp*PStorMax # Level of P Storage

    	theta_sp = LOpt_sp * kPhoto/(alpha_sp*(S_sp+EOpt_sp+LOpt_sp)) # Chlorophyl to Carbon

    # Calculate NP_sp and structural NP. Note the abscence of storage term in the latter. Probably can disregard NP_struct_sp for Cobalt version unless you want to 
    # explicitly track storage

    NP_sp = ((S_sp*NS_sp+NProt*(EOpt_sp+LOpt_sp))/(PE*EOpt_sp+PExtra_sp+PStor_sp*(S_sp) ))*molarP/molarN

    NP_struct_sp = ((S_sp*NS_sp+NProt*(EOpt_sp+LOpt_sp))/(PE*EOpt_sp+PExtra_sp ))*molarP/molarN
    
    IrrLimInst = (PCM_s*LOpt_sp*(1-exp(-aph*phiM*c*4.57*IrrInst/PCM_s))-bresp)/(1+zeta_sp)
    
    if EOpt_sp<=0.0
	    mu_sp = PCM_s*LOpt_sp*((1-exp(-aph*phiM*c*4.57*IDay/PCM_s)))*DayLength*sperd-bresp*sperd
    else
	    mu_sp = EOpt_sp*kE*sperd
    end

    theta = thetaFac*LOpt_sp*c/(LOpt_sp+EOpt_sp+S_sp)
    #print([aPhoto,PCM_s,c,k1*(1-c),k2*c,ENutLim_sp,LNutLim_sp,no3lim_sp,nh4lim_sp,nlim_sp],"\n")
    spSol = @SVector [EOpt_sp,LOpt_sp,mu_sp,NP_sp,NP_struct_sp,nutlimFlag_sp,theta,IrrLimInst,LOpt_sp*c,nlim_sp];

end

