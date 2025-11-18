Validation dataset experiment


I. FORMAT

1. 48 compounds at 5 different concentrations
2. Quadruplicates for each interaction tested
3. Promoters: 3 sets of promoters that will tell us different things:
3.1. The same all libscreening promoters tested (26 promoters + EVC control)
3.2. New promoters, not tested in libscreening (19 promoters)
3.3. MarA and RamA single or double KOs backgrounds with 5 promoters each (EVC, MicF, ramA, marRAB and acrAB - 15 strains). 20 strains in total if considered the same promoters in the WT background among the session 3.1


II. EXPERIMENTAL DESIGN

Experimental design provided in the pdf "validation_overview.pdf" in this folder

The experiment was done using 4 biological replicates for each strain/promoter in a single 384 well plate, with 96 different conditions.

Each 384 well plate (found in Plates_OD or Plates_Lux - 183 total) provides information on 5 concentration for 16 compounds. Since we have 48 compounds, 3 plates were required to cover all compounds.

Therefore, for each promoter, metadata will be found with the following code:

e.g. PacrAB_vp1_OD, PacrAB_vp2_OD, PacrAB_vp3_OD

PacrAB - promoter tested
vp1/vp2/vp3 - each set of compounds tested. Name was chosen to differentiate from our srn_code "lp1", "lp2" etc. Every plate has unique compounds/concentrations
_OD or _Lux - OD and Lux curves for each interaction

IMPORTANT - all experiments were done with the same robot method used in the libscreening!
That means we have the same number of timepoints and incubation times (40 min) and output (OD and lux curves) using the same reader calibration settings


III. LIBRARY MAPS

The conditions for each well of each vp can be found at the excel tables/.txt files in this folder:

vp1 - Map_plate1
vp2 - Map_plate2
vp3 - Map_plate3

These tables come with 10 common columns:

$well - what is the well/position (e.g. A1, A2 ... P24 etc)

$drug - what is the name of the drug tested (e.g. Isoconazole, Trimethoprim etc)

$conc - concentration tested (in uM). The concentrations were defined accordingly to the growth effect in the libscreening and are specific to each compound. 3 variations are possible:
0 to 150 - 0, 30, 60, 90, 120 and 150 (no signifficant growth effect libscreening) - n=39
0 to 100 - 0, 20, 40, 60, 80 and 100 ("average" midgrower) - n=2
0 to 50 - 0, 10, 20, 30 and 50 - ("harsh" midgrower) - n=7

$concMock - padronization of the concentration in a 1 to 5 format, regardless of the range tested

$quad - quadrant of the 384 wp in the 96wp format (1 to 4). Matches the replicate

$rep - replicate for each interaction (r1 to r4). Matches the quadrant

$abreviation - 4 letters abreviation for each compound name.

$wellCode - Information of the compound, concentration and replicate for each well
e.g. "Isoc_150_r1"

where
"Isoc" - abreviation for "Isoconazole" - compound
"150" - concentration tested (150 uM)
"r1" - replicate tested for this condition. That means one can find 4 wells with the "Isoc_150" treatment in this plate

$lib_well - well code in the vp format. e.g. "vp1_A1", "vp1_A2" ... "vp3_P24"

$srn_code - compound well in our libscreening data! In case one wants do make parallels between datasets. e.g. "lp6_I16" (Isoconazole (nitrate)), "lp3_B1" (Zinc Pyrithione) etc



