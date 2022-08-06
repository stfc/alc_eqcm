EQCM testing platform 
----------------------
----------------------

Author: i.scivetti  Oct 2020

ALC_EQCM uses the following tests, which provide reference data to be compared against. These tests are stored in ".tar" format. Via the do_test macro, these .tar files are unpacked automatically upon the execution of any to the xxx-test.sh inside the folder "tools" 
and run via Cmake in the corresponding directory. 

Description of tests
--------------------
####################
# EQCM data analysis
####################
test  1: Calibration of EQCM device for electrodeposition of Sb on Pt (printing raw EQCM data)
test  2: Calibration of EQCM device for electrodeposition of Sb on Pt (calculation of spectra in reciprocal space)
test  3: Calibration of EQCM device for electrodeposition of Sb on Pt (noise filtering and printing)
test  4: Electrodeposition of Ag on Pt (noise filtering and printing)
test  5: Intercalation/De-intercalation of Li+ ions into/from a (NiOH)2 host (characterization of the EQCM response)
test  6: Electrodeposition of Ag on Pt (mass calibration)
test  7: Intercalation/De-intercalation of Li+ ions into/from a (NiOH)2 host (characterization of the accumulated charge only)
test  8: Intercalation/De-intercalation of Li+ ions into/from a (NiOH)2 host (Print raw EQCM values between a given voltage range)
test 25: Compute massogram for intercalation of K+ ions into/from a (NiOH)2 host

########################################
# Solution of the stoichiometric problem
########################################
test  9: Intercalation/De-intercalation of Li+ ions into/from a (NiOH)2 host (Stoichiometry analysis: only Li+ and H+ participate in the reaction as dependent variables)
test 10: Intercalation/De-intercalation of Na+ ions into/from a (NiOH)2 host (Stoichiometry analysis: only Na+ and H+ participate in the reaction as dependent variables)
test 11: Intercalation/De-intercalation of Li+ ions into/from a (NiOH)2 host (Stoichiometry analysis: Li+ and H2O (dependent) and H+ (independent) participate in the reaction)
test 12: Intercalation/De-intercalation of Li+ ions into/from a (NiOH)2 host (Stoichiometry analysis: Li+ and H2O (dependent) and H+ (independent) participate in the reaction. No applied offset in the EQCM current)
test 13: Intercalation/De-intercalation of Na+ and Li+ ions into/from a (NiOH)2 host (Stoichiometric analysis for a fictitious example. Na+, Li+, H2O and H+ participate in the reaction. H2O and Li+ are dependent variables. Na+ and H+ are independent)
test 14: Electrodeposition of Sb on Pt
test 15: Intercalation/De-intercalation of K+ ions into/from a (NiOH)2 host  (Stoichiometry analysis: only K+ and H+ participate in the reaction as dependent variables)
test 16: Stoichiometric analysis with constraints: Intercalation/De-intercalation of Li+ ions into/from a (NiOH)2 host. Test for target_value 
test 17: Stoichiometric analysis with constraints: Intercalation/De-intercalation of Li+ ions into/from a (NiOH)2 host. Test for target_value and keep_ratio
test 18: Stoichiometric analysis with constraints: Intercalation/De-intercalation of K ions into/from a (NiOH)2 host. Test for target_min and keep_ratio
test 19: Stoichiometric analysis with constraints: Intercalation/De-intercalation of K ions into/from a (NiOH)2 host. Test for target_max and keep_ratio
test 20: Stoichiometric analysis with constraints: Intercalation/De-intercalation of K ions into/from a (NiOH)2 host. Test for target_range
test 21: Stoichiometric analysis with constraints: Intercalation/De-intercalation of K+ and Na+ ions into/from a (NiOH)2 host. Test for target_min and ratio_fixed
test 22: Stoichiometric analysis with constraints: Intercalation/De-intercalation of K+ and Na+ ions into/from a (NiOH)2 host. Test for target_min and ratio_fixed (settings differ from test 21)
test 23: Stoichiometric analysis with constraints: Intercalation/De-intercalation of K+, Na+ and Li+ ions into/from a (NiOH)2 host. Test for target_min, ratio_fixed, Target_value and keep_ratio
test 24: Stoichiometric analysis with constraints: Intercalation/De-intercalation of Li+ ions into/from a (NiOH)2 host. Test for constraint(s): target_value and ratio_fixed

###########
# Massogram
###########
test 25: Computation of massogram: Intercalation/De-intercalation of K+ ions into/from a (NiOH)2 host

########################
# Building atomic models 
########################
test 26: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of the pristine structure with target stoichiometry as set in the &Block_Species. Removing K+ and inserting water from an input model. 
test 27: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of the pristine structure with target stoichiometry [(NiO2)H_0.03 (H2O)_0.0023 K_0.654] 
test 28: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of the pristine structure with target stoichiometry [(NiO2)H_1.965 (H2O)_0.865 K_0.0885]
test 29: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of the pristine structure with target stoichiometry [(NiO2) (H2O)_0.971 Li_0.6362]
test 30: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of the pristine structure with target stoichiometry as set in the &Block_Species. Species involved are water, benzene, caffeine, Li+ and H+ (vasp format)
test 31: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of the pristine structure with target stoichiometry as set in the &Block_Species. Species involved are water, benzene, caffeine, Li+ and H+ (xyz format)
test 32: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of the pristine structure with target stoichiometry [(NiO2)H_0.20 (H2O)_0.33 K_0.354] 
test 33: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of a cycled sample with stoichiometric solutions computed via constraints. Li+, H2O and H+ participate in the reaction  
test 34: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of a cycled sample with stoichiometric solutions computed via constraints. K+, H2O and H+ participate in the reaction. Format for the generated models: cp2k
test 35: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of a cycled sample with stoichiometric solutions only for the first oxidation. 7 models are generated. Li+, H2O and H+ participate in the reaction. 7 models are generated.
test 36: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of a cycled sample with stoichiometric solutions only for the first oxidation. Na+, Li+, H2O and H+ participate in the (fictitious) reaction. 17 models are generated.
test 37: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of a cycled sample with stoichiometric solutions only for the first cycle. Na+, K+, H2O and H+ participate in the (fictitious) reaction. Only a single solution for oxidation. 17 models are generated for the subsequent reduction.
test 38: Electrodeposition of Sb on Pt - Building an atomistic model of the pristine structure with target stoichiometry as set in the &Block_Species.
test 39: Electrodeposition of Sb on Pt - Building an atomistic model of a cycled surface with stoichiometric solutions for the first 5 cycles. 
test 40: Electrodeposition of Al2O3 on Pt  (coating) - Building an atomistic model of the pristine structure with target stoichiometry as set in the &Block_Species.
test 41: Electrodeposition of Al2O3 on LMO (coating) - Building an atomistic model of the pristine structure with target stoichiometry as set in the &Block_Species. LMO is modelled with a reduction of Li, such that Li_{0.91667}Mn_{2}O_{4} 
test 42: Intercalation (bulk- nickel oxy-hydroxide) - Building an atomistic model of a cycled sample with stoichiometric solutions computed via constraints. K+, H2O and H+ participate in the reaction. The structure of the host has vacancies: Ni(OH)_{1.91667}.
test 43: The same as example 34, but defining a different atomistic composition for NiO2.
test 44: Electrodeposition of Li on amorphous Al2O3 - Building an atomistic model of the pristine structure with target stoichiometry as set in the &Block_Species.
test 45: Electrodeposition of Sb on Al2O3 - Building an atomistic model of a cycled surface with stoichiometric solutions for the first 2 cycles.
test 46: The same as test35 but both input and output structures are in onetep format 
test 47: The same as test35 but both input and output structures are in castep format
 
##########################################
# Building input files for DFT simulations 
##########################################
test 48: Modification of test 41. DFT files are in CP2K format 
test 49: Test 33, keeping the size of the input structure. Files are generated in CP2K format, only for the first cycle. Script files for the execution of HPC simulations in SCARF are also generated.
test 50: Test 35, only 3 models generated. Files are generated in CP2K format. 
test 51: Test 48 adapted to generate simulations files for CASTEP execution. 
test 52: Test 49 adapted to generate simulations files for CASTEP execution. 
test 53: Test 50 adapted to generate simulations files for CASTEP execution. 

####################################################
# Building an atomistic model of a disordered system 
####################################################
test 54: a single Li atom in a box of 32 water molecules

