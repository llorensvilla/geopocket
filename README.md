# Welcome to Geopocket

GeoPocket is a geometry based binding site predictor that combines surface and cavity detection with residue interaction scoring. It relies on classical geometric detection algorithms to identify concave regions on protein surfaces (Simões et al., 2017) and complements them with interaction propensities derived from surface triplet analysis (Mehio et al., 2010). 

## Installation
1. Create the working environment with conda:

``conda env create -f environment.yml ``
  `` conda activate geopocket ``

2. Install geopocket as a local package (using setup.py):
``pip install . ``

After these steps, a command called geopocket will be generated and you can run it from any folder (as long as the geopocket environment is active).

   
## References 
- Mehio, W., Kemp, G. J., Taylor, P., & Walkinshaw, M. D. (2010). Identification of protein binding surfaces using surface triplet propensities. Proteins: Structure, Function, and Bioinformatics, 78(7), 1701–1712. https://doi.org/10.1002/prot.22623
- Simões, T., Lopes, D., Dias, S., Fernandes, F., Pereira, J. M. B., Jorge, J., Bajaj, C., & Gomes, A. J. (2017). Geometric detection algorithms for cavities on protein surfaces in molecular graphics: A survey. Computer Graphics Forum, 36(8), 643–683. https://doi.org/10.1111/cgf.13158
