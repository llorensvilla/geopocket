background solid black
open /home/llorenc/PYT_SBI/Analisis/1brp_results/1brp.pdb
surface #0
surfrepr mesh #0

open /home/llorenc/PYT_SBI/Analisis/1brp_results/predicted_pockets_1.pdb
select #1 & :POC
rep sphere sel
color yellow sel
~select
focus

open /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_1_residues_1.pdb
surface #2
transparency 50 #2
color red #2

open /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_2_residues_1.pdb
surface #3
transparency 50 #3
color green #3

open /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_3_residues_1.pdb
surface #4
transparency 50 #4
color blue #4

open /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_4_residues_1.pdb
surface #5
transparency 50 #5
color magenta #5

open /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_5_residues_1.pdb
surface #6
transparency 50 #6
color cyan #6

open /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_6_residues_1.pdb
surface #7
transparency 50 #7
color orange #7

open /home/llorenc/PYT_SBI/Analisis/1brp_results/pocket_7_residues_1.pdb
surface #8
transparency 50 #8
color purple #8

focus
