 reho_csf takes a time series NIFTI volume set and runs a fourier transform 
 over the voxels highlighted in the mask volume

Building and installing reho_csf is done using the install-sh script 

Requires a Linux or Mac operating system with C and C++ compilers that allow for opemp

If your  using a Mac and your having trouble finding a C/C++ compiler that has openmp 
it's recommended that you use homebrew to download gcc and g++ versions that use openmp

To specify a specific C or C++ compiler and see all install-sh options use command:
install-sh --help
