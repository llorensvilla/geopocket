# Welcome to Geopocket!

GeoPocket is a geometry-Python based binding site predictor that combines surface and cavity detection with residue interaction scoring. It relies on classical geometric detection algorithms to identify concave regions on protein surfaces (Simões et al., 2017) and complements them with interaction propensities derived from surface triplet analysis (Mehio et al., 2010). 

With a single command line invocation, Geopocket will:

- Automatically adapt grid resolution to your protein size.
- Highlight and rank potential pockets by combining a Difference of Gaussians geometric score with empirically derived interaction scores.

- Output PDB files for pocket centroids and  residues, plus ready to run PyMOL and Chimera scripts for instant visualization.

For more details see the Geopocket_Tutorial which contains the installation theory and analysis of the Geopocket program. 

## Installation

1. Clone our repository:
All the material is available in our GitHub repository. For GeoPocket to work properly, it is recommended to clone the whole repository to have:
  - The main prediction file.
  - The setup.py file and the src/ folder with the code,
  - The environment.yml file contains the necessary dependencies.
PD: For the visualization it is necessary that you have previously installed in your operating system of preference, the programs Chimera and Pymol.

To clone the repository:
``git clone https://github.com/llorensvilla/geopocket.git``
``cd geopocket ``


3. Create the working environment with conda:
   ``conda env create -f environment.yml ``
  `` conda activate geopocket ``

4. Install geopocket as a local package (using setup.py):
``pip install . ``



After these steps, a command called geopocket will be generated and you can run it from any folder (as long as the geopocket environment is active).

   
## References 
- Mehio, W., Kemp, G. J., Taylor, P., & Walkinshaw, M. D. (2010). Identification of protein binding surfaces using surface triplet propensities. Proteins: Structure, Function, and Bioinformatics, 78(7), 1701–1712. https://doi.org/10.1002/prot.22623
- Simões, T., Lopes, D., Dias, S., Fernandes, F., Pereira, J. M. B., Jorge, J., Bajaj, C., & Gomes, A. J. (2017). Geometric detection algorithms for cavities on protein surfaces in molecular graphics: A survey. Computer Graphics Forum, 36(8), 643–683. https://doi.org/10.1111/cgf.13158
