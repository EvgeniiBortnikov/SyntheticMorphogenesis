Here, we provide a MATLAB code used for numerical simulations within three-component reaction-diffusion model.

************
Description
************
A 2D numerical model is based on three primary components – substrate (S), autocatalyst (A), and fast inhibitor (I) – including their diffusion, supply, and reactions. The decrease in the diffusivity of the activator species was approximated as a fast interaction of the activator with the gel matrix (G), resulting in the reversible formation of immobilized gel-bonded activatory species (GA). The reaction medium is represented as a 300×300 pixel grid, with each pixel measuring 0.05 mm × 0.05 mm. At each simulation step, the concentration change of each reactant within a pixel was calculated as the cumulative effect of diffusion, chemical reactions, and exchange with a continuously stirred tank reactor (CSTR). Diffusion-driven concentration changes were computed using the diffusion coefficients and the local concentration gradients of each reactant. The exchange between the reaction medium and the CSTR was modeled based on the supply rate constant (k0 in the model, relating to the experimentally derived permeability coefficient) and the concentration difference between the medium and the CSTR. The chemical reaction component was described by chemical equations with their associated reaction rate constants. 


*******************
System requirements
*******************
The MATLAB script was successfully tested using Windows 10 and Windows 11 operating system and MATLAB software version 2021a and newer.


************
Instructions
************
- Install the latest version of MATLAB software
- Copy "Three-component-model.m" and "Template.jpg" into the same folder
- Open "Three-component-model.m" with MATLAB software
- Adjust parameters of the simulation (simulation time, grid size, concentration of components, kinetic rate constants, etc.)
- Run the script. The figures, containing the updating concentration of free and bound activator and/or concentrations of substrate, fast inhibitor, free activatior and free and bound activator will be shown.
- After reaching the end of simulation, the recorded video of the change in the concentration profiles of the components will be created in the folder with the MATLAB file
  

************
Demo scripts
************
We added three demo scripts with pre-set parameters (concentrations, rate constant, simulations times, and etc.) which were used in order to achieve the formation of patterns depicted in Figures 4e-g of original work. 
Running those scripts will generate video recordings that you can also find in the demo scripts folder. The completion of the script on "average" computer will take around 5-15 minutes.
