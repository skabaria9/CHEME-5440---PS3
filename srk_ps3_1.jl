#How to add the package
#add https://github.com/varnerlab/CoreEcoliModelKit.git

#KEGG Diagram
#https://www.genome.jp/kegg-bin/show_pathway?map00220+C00624

#Converting from mM to mmol/gDW
vol_cell = 3.7e-12;                                       #L/cell (HeLa) BIND: 105879
fraction_cell_water = 0.798;                               # RBC, BIND: 101723
mass_of_single_cell = 2.3e-9;                              # g/cell BIND:103720
dryweight_cell = (1-fraction_cell_water)*mass_of_single_cell;
convert_factor(uM) = uM*(vol_cell)/dryweight_cell       #umol/gDW

#the metabolites 1 through 8 are described in attached Latex file

#Stoichiometic Matrix
stoichiometric_matrix = [
-1 0 0 0 0 0 0 1 0 0;
1 -1 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 -1 0;
0 1 -1 0 1 -1 0 0 0 0;
0 0 1 0 0 0 0 0 0 -1;
0 0 1 -1 0 0 0 0 0 0;
0 0 0 -1 0 0 0 0 0 0;
-1 0 0 1 -1 1 0 0 0 0;
];

#---------------------------------------------------------------------------
#Calculate upper and lower bounds of FLuxes
#-------------------------------------------------------------------------

#Lower Bound
flux_lower_bound = [0;0;0;0;0;0;0;0;0;0];

#upper bound
k_cat = [203;34.5;249;88.1;13.7;13.7;0;0;0;0]; #kcat for v1 to v5 (0 as placeholder for b's),s^-1
E = 0.01; #umol/gDW, enzyme concentration ratio, given in PS3
theta = 1; #regulatory control, assume is 1

metabolite_concentration = [ #uM converted to ummol/gDW
convert_factor(1.49E-2 * 10^6); #Aspartate
convert_factor(0);#Arginosuccitate
convert_factor(2.88E-4 * 10^6); #Fumarate, for E. coli (don't need it)
convert_factor(2.55e-4 * 10^6); #Arginine
convert_factor(0);#Urea, (don't need it)
convert_factor(0);#Ornithine
convert_factor(0);#Carbarmoyl Phosphate
convert_factor(0);#Citruline
];

Km = [ #uM converted to umol/gDW
convert_factor(0.18e-3 * 10^6); #v1, 6.3.4.5, Aspartate (BRENDA Mus musculus)
convert_factor(0.056e-3 * 10^6); #v1, 6.3.4.5, Citruline, (BRENDA, H. sapien, wild type)
convert_factor(0.1e-3 * 10^6); #v2, 4.3.2.1, Arginosuccitate, (BRENDA H. sapien)
convert_factor(2.05e-3 * 10^6); #v3, 3.5.3.1, Arginine, (BRENDA, Rattus norvegicus)
convert_factor(0.8e-3 * 10^6); #v4, 2.1.3.3, Ornithine, (BRENDA, Pseudomonas aeruginosa)
convert_factor(0.13e-3 * 10^6); #v4, 2.1.3.3, Carbamoyl Phosphate, (BRENDA, Homo sapiens, wild type)
convert_factor(0);#v5,1, 1.14.14.39, Citruline
convert_factor(3.50E-6 * 10^6); #v5,-1, 1.14.14.39, Arginine
];
];

saturation_term = [
(metabolite_concentration(1)/(Km(1)+metabolite_concentration(1)))*(metabolite_concentration(6)/(Km(2)+metabolite_concentration(6))); #v1, 6.3.4.5
(metabolite_concentration(2)/(Km(3)+metabolite_concentration(2))); #v2, 4.3.2.1
(metabolite_concentration(4)/(Km(4)+metabolite_concentration(4))); #v3, 4.3.2.1
#(metabolite_concentration(2)/(Km(5)+metabolite_concentration(2)))*(metabolite_concentration(2)/(Km(5)+metabolite_concentration(2))); #v4, 2.1.3.3
#(metabolite_concentration(2)/(Km(3)+metabolite_concentration(2))); #v5,1 1.14.14.39
#(metabolite_concentration(2)/(Km(3)+metabolite_concentration(2))); #v5,-1 1.14.14.39
]

v1_upper_bound = k_cat[1]*E;
v2_upper_bound = k_cat[2]*E;

flux_upper_bound =0;

#default_bounds_array = [flux_lower_bound flux_upper_bound];

#Species Bounds Array
species_bounds_array = [0;0;0;0;0;0;0;0;0;0];

#Objective Coefficient Array
species_bounds_array = [0;0;0;0;0;0;0;0;0;1];
