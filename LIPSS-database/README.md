# How to match experimental data directly with L_spp ?

1.1. Select a paper, extract the SEM pictures
1.2. Set the scale (ImageJ manual ? Python auto? --> Scripted Python-ImageJ interfaced)
1.3. Capture the thickness of orientation distribution. 
1.4. Get metadata: magnification, distance / area of measurement, laser wavelength, material
1.5. Add the in ExperimentalData.csv files. 
TODO: put everything in one single file. Script will separate them per material for plotting (current GNUplot script). Change to Python. 
1.6. Plot the new delta(angle) vs calculated L_spp. 

# ORganization of the data

1. Dielectric permittivities 
2. SEM pictures
3. Analyzed and extracted information
4. SUmmary into the automatically generated paper picture. 
