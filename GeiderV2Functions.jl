# Now for the Geider Version. 


function cOptPCM(I,tfac_Photo,liebig_lim)

        ICrit =  494.1251709153171*tfac_Photo*liebig_lim
	cCrit = 0.3034825870646766

	if I > ICrit
		return cCrit
	else
		ITrans = 2*I/ICrit-1

		a0 = 0.611510656625
		a1 = 0.22064845650928333
		a2 = -0.05382891709374325
		a3 =  0.01726321121940772
		a4 = -0.006761491369354537
		#print([0.5*a0, a1 * ITrans, a2 * (2*ITrans^2-1) , a3 * (4*ITrans^3-3*ITrans) , a4*(8*ITrans^4-8*ITrans^2+1)],"\n")


		return 1-sqrt(0.5*a0 + a1 * ITrans + a2 * (2*ITrans^2-1) + a3 * (4*ITrans^3-3*ITrans) + a4*(8*ITrans^4-8*ITrans^2+1))
	end
end







function cOpt(I,tfac_Photo)

        ICrit = 1482.3755127459513*tfac_Photo
	cCrit = 0.3034825870646766
	if I > ICrit
		return cCrit
	else
		ITrans = 2*I/ICrit-1

		a0 = 0.611510656625
		a1 = 0.22064845650928333
		a2 = -0.05382891709374325
		a3 =  0.01726321121940772
		a4 = -0.006761491369354537
		#print([0.5*a0, a1 * ITrans, a2 * (2*ITrans^2-1) , a3 * (4*ITrans^3-3*ITrans) , a4*(8*ITrans^4-8*ITrans^2+1)],"\n")

		return 1-sqrt(0.5*a0 + a1 * ITrans + a2 * (2*ITrans^2-1) + a3 * (4*ITrans^3-3*ITrans) + a4*(8*ITrans^4-8*ITrans^2+1))
	end
end


function cOptDiPCM(I,tfac_Photo,liebig)

        ICrit =  494.1251709153171*tfac_Photo*liebig

	cCrit = 0.3034825870646766



	if I > ICrit
		return cCrit
	else
		ITrans = 2*I/ICrit-1

		a0 = 0.611510656625
		a1 = 0.22064845650928333
		a2 = -0.05382891709374325
		a3 =  0.01726321121940772
		a4 = -0.006761491369354537
		return 1.0-sqrt(0.5*a0 + a1 * ITrans + a2 * (2*ITrans^2-1) + a3 * (4*ITrans^3-3*ITrans) + a4*(8*ITrans^4-8*ITrans^2+1))
	end
end



function cOptDi(I,tfac_Photo)

        ICrit =  494.1251709153171*tfac_Photo

	cCrit = 0.3034825870646766
	if I > ICrit
		return cCrit
	else
		ITrans = 2*I/ICrit-1

		a0 = 0.611510656625
		a1 = 0.22064845650928333
		a2 = -0.05382891709374325
		a3 =  0.01726321121940772
		a4 = -0.006761491369354537

		return 1.0-sqrt(0.5*a0 + a1 * ITrans + a2 * (2*ITrans^2-1) + a3 * (4*ITrans^3-3*ITrans) + a4*(8*ITrans^4-8*ITrans^2+1))
	end
end






# First we define a parameterization function. This will enable full simulation of the Cobalt food web model

function ZeroDInitGeiderV2()

#
#   COBALT parameter values (Default if global formulation of Stock et al.,
#   2014.  Note that all rates are given for 0 deg. C
#

sperd = 86400;                      # sec day-1
params= Dict([])
# sinking

params["mld"] = 50.0 
params["kpar"] = 0.04
params["par"] = 100.0
params["par24"] = 100.0
params["PStorMax_sp"] = 0.01

params["PStorMax_lp"] = 0.025
params["PStorMax_di"] = 0.01
params["LFac"]=1.2    
# These are units of gC/gC/hour

params["QF"]=0.17   # Nitrogen/Carbon ratio
params["k1"] = 6.1e-4*params["QF"]*3600*2.0
params["k2"] = 1.4e-3*params["QF"]*3600*2.0

params["aph_sm"] = 11.6*params["QF"] # m^2 /gC

params["aph_di"] = 11.6*params["QF"]/3.0 # m^2 /gC
params["aph_lp"] = 11.6*params["QF"]/3.0 # m^2 /gC
params["phiM"] = 1.0e-6  # gC/mu mol photons


    
    
params["wsnk"] = -100/sperd;           # m sec-1 (negative is down)
# refuge concentration for basal metabolic rates
params["ref_conc"] = 1e-5;             # mmoles m-3
# Phytoplankton
params["kappa_phyto"] = 0.063;        # deg C-1
params["kappa_photo"] = 0.063;
params["alpha_sp"] = 2.0e-5*2.77e18/6.022e17; # g C g Chl-1 m2 W-1 s-1..
params["alpha_lp"] = 1.0e-5*2.77e18/6.022e17; # g C g Chl-1 m2 W-2 s-1
params["alpha_di"] = 2.0e-5*2.77e18/6.022e17;
params["P_C_max_sp"] = 1.25/sperd;   # sec-1
params["P_C_max_lp"] = 1.25/sperd;    # sec-1
params["P_C_max_di"] = .6/sperd
params["thetamax_sp"] = 0.03;         # g Chl g C-1
params["thetamax_lp"] = 0.05;         # g Chl g C-1
params["thetamax_di"] = 0.05;
params["zeta_sp"] = 0.05;             # dimensionless
params["zeta_lp"] = 0.05;              # dimensionless
params["zeta_di"] = 0.05;


params["bresp_sp"] = 0.03/sperd;    # sec-1
params["bresp_lp"] = 0.05/sperd;   # sec-1
params["bresp_di"] = 0.05/sperd;
params["thetamin"] = 0.002;           # g Chl g C-1
params["PExtra_sp"] = 0.0035;
params["PExtra_lp"] = 0.005;

params["PExtra_di"] = 0.003;

# Here we will calculate the uptake parameters of the small and large
# phytoplankton. These will be based on the radius of the respective
# plankton types

#params["r_sp"] = 0.5 ;                  # radius of small phytoplankton in m-6
params["r_lp"] = 5.0  ;                 # radius of thelarge phytoplankton in m-6

params["r_di"] = 1.0  ;                 # radius of thelarge phytoplankton in m-6
params["r_sp"] = 1.0  ;                 # radius of thelarge phytoplankton in m-6

#params["r_di"] = 1.5
# We will calculate the half-saturation constant based on allometry and the
# radius
params["kno3_sp"] = 0.25;              # mmoles NO3 m-3
params["kno3_lp"] = 2.5;              # mmoles NO3 m-3
params["knh4_sp"] = 0.05;              # mmoles NH4 m-3
params["knh4_lp"] = 0.5;              # mmoles NH4 m-3
params["knh4_di"] = 0.25
params["kno3_di"] = 2.5


params["volume_sp"] = 4.0/3.0*pi*params["r_sp"]^3;

params["volume_lp"] = 4.0/3.0*pi*params["r_lp"]^3;

# Now we define the phosphorus uptake parameters, in a similar manner:

params["kpo4_sp"] = 0.01 ;
params["kpo4_lp"] = 0.1;
params["kpo4_di"] = 0.1;


params["k_p_stor_sp"] = .025
params["k_p_stor_lp"] = 1.0

params["k_p_stor_di"] = 1.0 

params["DayLength"] = 0.5
params["MeanIrradiance"] = 1000.0 
params["aFe"] = 10/100/(sperd*365.0)*1e-1

#


params["q_fe_k_sp"] = 3.0e-6*(106/16.0)
params["q_fe_k_lp"] = 6.0e-6*(106/16.0)
params["q_fe_k_di"] = 25.0e-6*(106/16.0)

params["k_fed_sp"] = 1e-4
params["k_fed_lp"] = 5e-4
params["k_fed_di"] = 5e-4

params["q_fe_sp_max"] = 50e-6*(106/16.0)
params["q_fe_lp_max"] = 50e-6*(106/16.0)
params["q_fe_di_max"] = 50e-6*(106/16.0)

params["FeN_sp_upt_fac"] = 30.0e-6
params["FeN_lp_upt_fac"] = 30.0e-6
params["FeN_di_upt_fac"] = 30.0e-6

params["thetaFactorSm"] = params["aph_sm"]*params["phiM"]/(6e-5*2.77e18/6.022e17)/params["QF"]/2.0
params["thetaFactorLg"] = params["aph_lp"]*params["phiM"]/(1.0e-5*2.77e18/6.022e17)/params["QF"]/2.0
params["thetaFactorDi"] = params["aph_di"]*params["phiM"]/(1.0e-5*2.77e18/6.022e17)/params["QF"]/2.0



# We make sure some other params["are defined

params["E_sp"] = .1;     # Biosynthesis allocation for small phytoplankton (placeholder)
params["E_lp"] =.1;     # Biosynthesis allocation for large phytoplankton (placeholder)
params["NPModelFlag"] = 1;

params["m_agg_sp"] = 0.1/sperd;       # sec-1 mmole N-1 m3
params["m_agg_lp"] = 0.3/sperd;       # sec-1 mmole N-1 m3
params["m_agg_di"] = 0.1/sperd;
params["m_vir_sp"] = 0.025/sperd;     # sec-1 mmole N-1 m3
params["m_vir_di"] = 0.05/sperd;
params["m_vir_lp"] = 0/sperd;         # sec-1 mmole N-1 m3
params["exu"] = 0.13;                 # dimensionless
params["mld_thresh"] = 0.01;          # kg m-3

# Zooplankton
params["kappa_zoo"] = 0.063;          # deg C-1
params["i_max_sz"] = 1.42/sperd;      # sec-1
params["i_max_mz"] = 0.57/sperd;      # sec-1
params["i_max_lz"] = 0.23/sperd;      # sec-1
params["ki_sz"] = 1.25;               # mmoles N m-3
params["ki_mz"] = 1.25;               # mmoles N m-3
params["ki_lz"] = 1.25;               # mmoles N m-3
params["gge_max_sz"] = 0.4;           # dimensionless
params["gge_max_mz"] = 0.4;           # dimensionless
params["gge_max_lz"] = 0.4;           # dimensionless
params["bresp_sz"] = 0.020/sperd;     # sec-1
params["bresp_mz"] = 0.008/sperd;     # sec-1
params["bresp_lz"] = 0.0032/sperd;    # sec-1

# detritus/dissolved organic matter partitioning, 30# goes to detritus,
# where it is partitioned between sinking and dissolved according to size.
# For the dissolved phase, partitioning between labile, semi-labile and
# semi-refractory is based on a tuning to global patterns.
params["phi_det_sz"] = 0.0;            # dimensionless
params["phi_det_mz"] = 0.20;            # dimensionless
params["phi_det_lz"] = 0.30;            # dimensionless
params["phi_ldon_sz"] = 0.30*0.57;      # dimensionless
params["phi_ldon_mz"] = 0.10*0.57;      # dimensionless
params["phi_ldon_lz"] = 0.0;            # dimensionless
params["phi_sldon_sz"] = 0.30*0.40;     # dimensionless
params["phi_sldon_mz"] = 0.10*0.40;     # dimensionless
params["phi_sldon_lz"] = 0.0;           # dimensionless
params["phi_srdon_sz"] = 0.30*0.03;     # dimensionless
params["phi_srdon_mz"] = 0.10*0.03;     # dimensionless
params["phi_srdon_lz"] = 0.0;           # dimensionless



params["phi_detp_sz"] = 0.0;
params["phi_detp_mz"] = 0.20;
params["phi_detp_lz"] = 0.30;
params["phi_ldop_sz"] = .8*.3;
params["phi_ldop_mz"] = .8*.1;
params["phi_ldop_lz"] = 0 ;
params["phi_sldop_sz"]= .2*.3;
params["phi_sldop_mz"]= .2*.1;
params["phi_sldop_lz"]= 0;
params["phi_srdop_sz"]= .0*.3;
params["phi_srdop_mz"]= .0*.1;
params["phi_srdop_lz"]= 0;


params["phi_detfe_sz"] = .3;
params["phi_detfe_mz"] = .3;
params["phi_detfe_lz"] = .3;

params["phi_fed_sz"] = .7;
params["phi_fed_mz"] = .7;
params["phi_fed_lz"] = .7;





# switching shape parameters, see: Stock, Powell and Levin (2008),
# Bottom-up and top-down forcing in a simple size structured plankton
# dynamics model.  Journal of Marine Systems 74(1) 134-152 for rationale
# and references.
params["ms_sz"] = 2;
params["ns_sz"] = 2;
params["ms_mz"] = 2;
params["ns_mz"] = 2;
params["ms_lz"] = 2;
params["ns_lz"] = 2;
# innate prey availability for sp, lp, di, sz, mz, lz, ndet; dimension nz x 7
ipa_matrix_sz = @SVector [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.25, 0.0];
ipa_matrix_mz = @SVector [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0];
ipa_matrix_lz = @SVector [0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0];

    

params["ipa_matrix_sz"] = ipa_matrix_sz;
params["ipa_matrix_mz"] = ipa_matrix_mz;
params["ipa_matrix_lz"] = ipa_matrix_lz;


# Higher predation closure
params["kappa_hp"] = 0.063;           # deg C-1
params["i_max_hp"] = 0.09/sperd;      # sec-1
params["ki_hp"] = 1.25;               # mmoles N m-3
params["ms_hp"] = 2;           # dimensionless
params["ns_hp"] = 2;           # dimensionless
ipa_matrix_hp = @SVector [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0];   # dimensionless

params["ipa_matrix_hp"] = ipa_matrix_hp;
params["coef_hp"] = 2;                # dimensionless
# partitioning of higher predator losses to different pools
params["phi_det_hp"] = 0.35;          # dimensionless
params["phi_ldon_hp"] = 0.0;          # dimensionless
params["phi_sldon_hp"] = 0.0;         # dimensionless
params["phi_srdon_hp"] = 0.0;         # dimensionless
params["phi_nh4_hp"] = 0.65;          # dimensionless


params["phi_detp_hp"] =0.35;
params["phi_ldop_hp"] =0;
params["phi_sldop_hp"] =0;
params["phi_srdop_hp"] =0;
params["phi_po4_hp"] =0.65;


params["phi_detfe_hp"] = 0.35;
params["phi_fed_hp"] = 0.65;



# Bacteria
params["kappa_bact"] = 0.063;         # deg C-1
params["mu_max_b"] = 1.0/sperd;       # sec-1             
params["kldon"] = 0.5;                # mmoles ldon m-3
params["gge_max_b"] = 0.4;            # dimensionless
params["bresp_b"] = 0.0075/sperd;     # sec-1
params["m_vir_b"] = 0.033/sperd;      # sec-1 mmole N-1 m3
# partitioning of virus losses amongest dissolved organics
params["phi_ldon_vir"] = 0.55;        # dimensionless
params["phi_sldon_vir"] = 0.40;       # dimensionless
params["phi_srdon_vir"] = 0.05;       # dimensionless


params["phi_ldop_vir"] =0.8;
params["phi_sldop_vir"] =0.2;
params["phi_srdop_vir"] =0.0;


# Detritus and dissolved organic material
params["gamma_det"] = 1.0/sperd;        # sec-1
params["gamma_srdon"] = 1.0/(18*365*sperd);    # sec-1
params["gamma_sldon"] = 1.0/(90*sperd);        # sec-1
params["gamma_nitrif"] = 1.0/(30*sperd);       # sec-1


params["gamma_detp"] = params["gamma_det"];
params["gamma_srdop"] =1.0/(18*365*sperd);
params["gamma_sldop"] =1.0/(90*sperd);



# Mixing params for 0D chemostat style model
params["N0"] = 50.0;
params["no30"] = 50.0;#mmole m3
params["nh40"] = 0.01; #mmole m3
params["ldon0"] =0.01;     #mmole m3
params["sldon0"] =0.05;    #mmole m3
params["srdon0"] =0.1 ;    #mmole m3
params["detSink"] =1.0/(1.0*sperd);     #1/s



params["po40"]=3.0;

params["ldop0"] =0.001;
params["sldop0"] =0.05;
params["srdop0"] =0.1;
params["QEm"] = 1.4;


params["fed0"] = 5e-4;
params["det_fe0"] = 0.0;

params["tau"] = 1.0/(sperd*365.0);
params["downwell"] = .00/(sperd*365);

params["np_b"] = 16.0;
params["kfe_eq_ll"] = 1.0e12;
params["kfe_eq_hl"] = 1.0e8;
params["io_fescav"] = 10.0;
params["felig_bkg"] = 1.0e-9;
params["alpha_fescav"] = 15.0/(sperd*365.0);

params["remin_eff_fedet"] = 0.5;

return params
end


