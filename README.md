# CRISTINA_script

Matlab script using CRISTINA_getMQC and CRISTINA_getFleysherMQC to compute the Multi quantum images, notably the Single Quantum (SQ) and Triple Quantum (TQ) images. Phantom measurement data are provided taken at 7T (SIEMENS Magnetom).

Multi-quantum sodium imaging offers additional insights compared to standard single-quantum (SQ) images with information beyond tissue sodium concentration. The characteristic temporal evolution of the SQ and triple quantum (TQ) sodium is of interest. A multi-echo 2D Cartesian Imaging sequence for SQ and TQ 23Na (CRISTINA) Phantom raw data, acquired at 7T (Siemens Magnetom , Erlangen, Germany) can be downloaded along with the reconstruction. 

The code corresponds to the results presented in:  # MRM-19-20481.R1.

## The following additional tools are neccessary:
- ff2c.m function of the SparseMRI V0.2 toolbox from:
  https://people.eecs.berkeley.edu/~mlustig/Software.html
  
- the pocs reconstruction from:
  https://fr.mathworks.com/matlabcentral/fileexchange/39350-mri-partial-fourier-reconstruction-with-pocs
  (version 1.2.0.0) by Michaela Völker
  
- the unwrap3.m function from:
  https://github.com/marcsous/unwrap/commits/master
  commit:50dac46f9295f108e41a5e86c879c5932ba60f47
  
- for visualization the following toolbox is nice and the function as.m was used:
  https://github.com/tsumpf/arrShow 
  


