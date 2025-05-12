# 2D KMC

## Introduction
Implementation of a kinetic monte-carlo and stochastic simulation process are present in two separate branches of this repository. The "main" branch contains the stochastic simulation, 
which is more robust for larger energy ranges, especially considering dehalogenation effects. The kinetic monte-carlo simulation is present in the KMC_Complete branch, and an identical procedure is used for that case.

## Operation
The simulations, in both the KMC and stochastic cases, are conducted by running the main.py python script. In order to configure the settings of the simulations, open the main.py file and 
navigate to the final function (titled "main"). The beginning section of this function is shown below:

```
 # Initialize lattice and monomers
    diffusion_energies = np.linspace(0,1.5,6)
    coupling_energies = np.linspace(0,1.5,6)
    dehalogen_energies = np.linspace(0, 2.3, 6)
    rotation_energies = [0]

    width = 60 # only even numbers
    num_simulations_per_triplet = 3
```

The range of energies use for the simulations are controlled by the first four variables, where one can edit the numpy linspace to the required values. In this case, rotation energies were 
restricted to 0 at all times, but one simply needs to replace the variable with a linspace exactly the same as the previous three variables in order to create a range of rotation energies.
Next, the width of the simulation grid can be set (the width and height are taken to be equivalent values). The variable "num_simulations_per_triplet" indicates the number of redundant simulations
to perform for each set of energy values. The results from each of these extra simulations will be averaged together when exporting the data. At the bottom of the "main" function is a call to
"save_results_to_csv", with a file path passed into the function, which can be edited.

It is important to note that the KMC simulation struggles with even small energy ranges with the addition of dehalogenation, which is an algorithmic problem that needs to be addressed in the future.

## Analysis
In order to compare the outputs of the simulations to experimental results, a set of self-contained analysis files is required. For image connectivity analysis of experimental STM images, the 
STM-Island-Analysis.ipynb file may be used. It is found in the "data" folder on the "main" github branch, and only requires a path to the image (a .bmp file) in order to function. Before this analysis,
it is important to crop the STM image down to a single island and to sharply contrast the background. The existing STM-Island-Analysis.ipynb file has an example of the proper format for the image 
loaded into it. The easiest way to perform this preliminary image refinement is to use Gwyddion, a free software designed for analysis of AFM and STM images (https://gwyddion.net/presentations/tutorial-Gwyddion-basic-data-correction-Francois-Riguet-2012.pdf).

After running STM-Island-Analysis.ipynb on the experimental polymer island, a .csv file will be created containing all of the node connectivity information for that particular image. This .csv file can be
used alongside the .csv output from the simulation to perform Wasserstein distance analysis on the two sets of connectivity information in the "comparing-simulated-experimental-histograms.ipynb" file
also within the "data" folder. At the bottom of this file, one must submit the path to the .csv containing the simulation results as well as the directory containing both the (cleaned) experimental 
STM image and the results from the STM-Island-Analysis.ipynb file. Simply run the file to perform Wasserstein analysis. 

Once the Wasserstein analysis is done, one can find the minimum Wasserstein distance to identify the best-fit energy parameters that produce a polymer island with node connectivity most similar to that
of the experimental data. One can also visualize these results using matplotlib.pyplot.contourf, which generates interpolated phase diagrams indicating Wasserstein distance minima. An example of the use of
this plotting function is in "contour.py" within the data folder. This file will have to be altered to fit the needs of each simulation-experiment comparison procedure, but it is primarily a visualization
tool.

### For further references, the undergraduate thesis outlining the general theory and code architecture is provided in pdf form.


