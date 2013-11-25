###########################################################################
# Eventually, all the global variables should be initialized in this file #
# This file should be sourced at the beginning of toad.tcl                #
###########################################################################

set PI 3.1415926535897932384626433832795

# For FLICAM CCD, intensity goes up to 64000, which requires int* to hold data and 
# manipulate the data in NFIT.
# Setting FLICAM to 0 will make the pointer to the input data short*.  
set FLICAM 1

# paramNames: TCL list that contains the names of fitting parameters
# Useful for accessing elements in paramArray and freeOrFixed arrays
set paramNames {Kc B Lr Mz D mosaic edisp bFWHM s bc2b wavelength \
	              pixelSize qxzero nindex T}
	              
# paramArray: TCL array for the values of the fitting parameters
# The following lines simply initialize the array and the entries in the 
# fitting panel
set paramArray(Kc) 8e-13
set paramArray(B) 2e13
set paramArray(Lr) 5000
set paramArray(Mz) 10
set paramArray(D) 62.8
set paramArray(mosaic) 0
set paramArray(edisp) 0.0134
set paramArray(bFWHM) 2.3
set paramArray(s) 359.7
set paramArray(bc2b) 0
set paramArray(wavelength) 1.175
set paramArray(pixelSize) 0.07113
set paramArray(qxzero) 512
set paramArray(nindex) 0.999998
set paramArray(T) 37

# freeOrFixed: TCL array, if an element is 1, the correspoinding variable specified 
# by paramNames becomes a fitting parameter. If 0, the variable value will be fixed to 
# whatever value specified in its corresponding entry field
# Index key for this array is a parameter name. Use paramNames list to access it 
# conveniently through all the elements
foreach i $paramNames {
	set freeOrFixed($i) 0
}
set freeOrFixed(Kc) 1
set freeOrFixed(B) 1

set nthreads 2
set noscale 1
set dupe 1
set nocz 1
set iter 1
