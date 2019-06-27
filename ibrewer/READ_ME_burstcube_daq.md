# Documentation for the Energy_Resolution_Scope code

## Energy_Resolution_Scope is a Jupyter Notebook that pulls traces off the scope and stores these traces in arrays. These arrays are saved in h5py files. h5py files are convenient because you can splice into them, which is nice for plotting and data analysis. Each row in a data set corresponds to one scope trace, so the data sets returned have a number of rows equal to the number of traces desired for each run (column number is equal to the number of data points the scope records for each trace, usually 10000).

## Note on communicating to the scope:
### Make sure the scope and the lab computer are both plugged into the router (but not the Internet port). This sets the IP address of the scope.  In the code, enter the IP address of the scope (as assigned by the router) as a string into the hander "scope = sdaq.Scope(address="10.10.10.2")." To check the IP address of the scope, you can go to the Utility menu and check the LAN settings. Sometimes a LAN reset is required.

## Naming scheme for files and data sets:
### The naming scheme for files is ***BurstCube_name of test_source name_date***, ex "BurstCube_PostVibe_Cs137_061419." You can create multiple data sets within a file (you can take multiple runs and save them to one file).
### The naming scheme for data sets is ***filename_run#***, ex. "BurstCube_PostVibe_Cs137_061419_run1."

## How this code determines energy resolution:
### Once an array has been created, the code looks at each row of the array (or each pulse) and extracts a peak value. It does this by finding the maximum  value of the pulse and then subracting a noise value. The peak height is important because it is (approximately) a measure of energy.
### The peak values of each pulse are stored in an array that is as large as the number of traces taken with the scope. These peak values are then histogrammed to determine the energy resolution (how finely the detector can sense, in this case, a radioactive source). The full width-half max (FWHM) gives an indication of the energy resolution-- that is, we want the width of the peak of the histogram to be narrow in comparison to the energy (mean peak value) of the source.
### This code uses the scipy.optimize function to fit a Gaussian to the histogram and return the amplitude, mean peak value, and sigma. The energy resolution is calculated using the equation ***energy resolution  = (2.355xSigmax100)/Mu.***
### After the energy resolution is calculated, it will also calculate the error in the energy resolution.

## Note on h5py files:
### Make sure, once you are done with a file, you close the file! h5py gets very upset if you try and open a new file or kill the kernel if an h5py file is still open. Use the ***f.close()*** function to close files.
### Documentation for h5py can be found at http://docs.h5py.org/en/stable/.