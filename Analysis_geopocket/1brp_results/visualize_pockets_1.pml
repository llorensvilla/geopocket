bg_color black
load /home/llorenc/PYT_SBI/Analisis/1brp_results/1brp.pdb, protein
hide everything, protein

show cartoon, protein
color gray70, protein
set cartoon_transparency, 0.2, protein

show mesh, protein
color brown50, protein
set mesh_as_cylinders, 1
set mesh_width, 0.6
set surface_quality, 2
set two_sided_lighting, on

load /home/llorenc/PYT_SBI/Analisis/1brp_results/predicted_pockets_1.pdb, pockets
hide everything, pockets
show spheres, pockets
set sphere_scale, 0.4, pockets

color yellow, pockets and id 1
color red, pockets and id 2
color green, pockets and id 3
color blue, pockets and id 4
color magenta, pockets and id 5
color cyan, pockets and id 6
color orange, pockets and id 7

load /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_1_residues_1.pdb, residue_1
hide everything, residue_1
show surface, residue_1
color yellow, residue_1
set transparency, 0.3, residue_1

load /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_2_residues_1.pdb, residue_2
hide everything, residue_2
show surface, residue_2
color red, residue_2
set transparency, 0.3, residue_2

load /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_3_residues_1.pdb, residue_3
hide everything, residue_3
show surface, residue_3
color green, residue_3
set transparency, 0.3, residue_3

load /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_4_residues_1.pdb, residue_4
hide everything, residue_4
show surface, residue_4
color blue, residue_4
set transparency, 0.3, residue_4

load /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_5_residues_1.pdb, residue_5
hide everything, residue_5
show surface, residue_5
color magenta, residue_5
set transparency, 0.3, residue_5

load /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_6_residues_1.pdb, residue_6
hide everything, residue_6
show surface, residue_6
color cyan, residue_6
set transparency, 0.3, residue_6

load /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_7_residues_1.pdb, residue_7
hide everything, residue_7
show surface, residue_7
color orange, residue_7
set transparency, 0.3, residue_7

zoom all
