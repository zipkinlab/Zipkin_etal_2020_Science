## Community occurrence model

### Code
This folder contains the code to run the community occurrence model. The wrapper includes code to calcuate species richness, individual species occurrence rates, and detection rates for both observed and unobserved species.

1. WRAPPER_COMMUNITY_MODEL.R      
Description: This file contains the code to run the community occurrence model (calling the file JAGS_MODEL.R). The code reads in the data (species observations and transect effort) to create the objects for the Bayesian community occurrence model. The code to process the model results is also included (e.g., for Figure 1a, Figure 2).

1. JAGS_MODEL.R      
Description: This file contains the JAGS code for the community occurrence model.
