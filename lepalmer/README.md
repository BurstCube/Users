This README document explains the outline of the larger documented README link for all information pertaining to the localization code making process!

The larger documented README can be found at this link: [ReadMe--Localization Code Information](https://docs.google.com/document/d/1x1kDDz5EkpokAVYvPcEUN7cDeyFaSucVGyu95ZPbbqQ/edit?usp=sharing)

Beginning Remark: **Information regarding what the directories in this directory are used for can be found in the larger README. In addition, if it occurs to the user that there are paths written to folders/files that are not correct, that is because when uploading the folders/files to gitHub, they will be coming directly from a local computer**.

The first section of the larger README, entitled Preliminary Process,
contains information on the preliminary process of making the localization code--specifically, figuring out how to create simulation files as well as making plots to understand how the effective area of one BurstCube detector compares to either a certain energy, a certain azimuth, or a certain zenith.

The second section, entitled Pre-Localization Codes, contains all of the preliminary codes created after determining how the effective area compared against the energy. The preliminary codes are placed in estimated order of completion.

The third section, entitled Extra comments about Preliminary Code, contains information about preliminary codes, but how they are related in the grand-scheme of things is unknown at the moment. This section could really be avoided. 

The fourth section, entitled Final Process (localDetFinder in SmallSteps/finalTest), contains information about A) How the code works and its documentation, B) Where the stored/simulated data and random burst data is stored, C) Background info pertaining to the process and understanding of the code, D) Comments on old data related to the code, E) Documentation on a script to easily update the configuration files relative to a certain problem, F) Extra code that used to be important relative to the final code, G) A listing of all the stored/simulated data with the random bursts employed, and H) An explanation for the two analyzing methods used on the stored/simulated data with the random bursts.

Final Code Summary Cell-by-cell:
The first portion of the code imports the sim files from the directory containing the configuration file, and the BurstCube module is used to extract elements from the files: energy, azimuth, and zenith. Astropy module is also imported to allow for creation of tables/updating tables.

The second portion of the code contains the detector counting variables, which are initialized per each sim file. The list hitsData is also initialized to store all data (either from one file to find the unknown location of a random burst or from multiple files to store simulated data in a table). Detector counts are found by extracting lines out of simulation files containing hit data (each hit should be located on the center of a detector because the source files set the DiscretizeHits variable to TRUE): each detector is located in each quadrant of the cartesian plane, the azimuth angle is measured counterclockwise from the x-axisand, the absolute value of each center x-y pair is: (5.52500, 5.52500)

The third portion/cell helps to analyze the hitsData array for future use.

The fourth portion/cell has multiple functions:
First, the user is asked whether an uncertainty factor will be applied to the data (if so, the entire cell is bypassed). If not, a series of statements allows the user to choose if they want to either determine the random burst location (inputDet will have the same value as hitsData), store simulated data loaded in from hitsData into a dictionary, or bypass all options. In the case that a dictionary is made, all elements partaining to the burst data are extracted: the energy, azimuth, zenith, and the numbers of hits per each detector (all per each simulation file).

The fifth portion/cell implements the uncertainty factor if the user has requested it, and an array variable trialArrayOH will contain all rows in the loaded table that have an energy of 100 keV. The table range to be extracted will be requested (using "first" and "second"), and then the number of hits per detector will be loaded from that data into the hitsArray variable. A list will then be generated corresponding to each set of detector data, and each element in each list will have the uncertainty factor applied, which is the square root of the number of triggers over the number of triggers (0.0316227766 or about 3 percent). The step size pertaining to the range of uncertainty, from a low of 3 percent lower than the current value to 3 percent higher than the current value, will also be inputted from the user. Uncertainty data will then be stored in a table.

The sixth portion/cell of the code has that if the random burst location is to be determined, the energy of the burst is requested by the user, and all rows containing that select data are generated. The user is then asked the zenith of the random burst. The user is then asked which elements are to be extracted from the table using their indexes ("begin" and "end") to compare against the random burst. The least squares calculation is then applied to the entire set of selected data against the random burst: the position where the least squares results in the lowest value is placed in the variable, "result". "result" will then be used to store burst location, categorizing it by specfied energy ("energyRB") and zenith ("zenithRB").

The last portion/cell contains a series of options that exist for the user in which case, if they are trying to determine the random burst location, they can either make a table for the first time to store the data, or they can add a new row, with their determined GRB info, to their pre-existing stored data.

Some Remarks:
On preliminary codes used, things that are in bold and say “comments”, are not actually comments like commented code--it is old code that is commented out.
If there is commented code not explained, or variables in the code that are not used/not mentioned, ignore that information or suggest further documentation.
In the user directory, there may be data in directories that are outliers: files or data not used. Please ignore!
Some data is unknown on where or how it was used because it was not documented at the time it was used--information regarding these inconsistencies will be found in the larger document.




