# How to match experimental data directly with L_spp ?

1. Select a paper, extract the SEM pictures using "Capture" tool in Linux. 

2. In ImageJ: set the scale using "Analyse > Set scale."

3. Installing the plugin OrientationJ (http://bigwww.epfl.ch/demo/orientation/).

4. Using OrientationJ > Orientation-Distribution tool, using the Fourier Gradient, capture the thickness of orientation distribution. 

5.  

6. Get metadata: magnification, distance / area of measurement, laser wavelength, material

7. Add the detailed data in ExperimentalData.csv files. 

TODO: put everything in one single file. Script will separate them per material for plotting (current GNUplot script). Change to Python. 

6. Scripted: Plot the new delta(angle) vs calculated L_spp. Scripts is here: ./SimpleSPP/Results/Air/Iaroslav/1030nm/plot.py.

# ORganization of the data

1. Dielectric permittivities 

2. SEM pictures

3. Analyzed and extracted information

4. SUmmary into the automatically generated paper picture. 

Done. 
