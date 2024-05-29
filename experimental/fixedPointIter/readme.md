implementation of water flow in the plant and soil with X solute dissolved in water [0,2]
and (optionally) up to 4 microbial groups, solute adsorption and CO2 release by microbial activity (based on the TraiRhizo model).
We get therefore up 9 components. The full model is called "10c" as we count the static soil as an element.


# folder structure
- scripts:
	- XcGrowth: main script, to run simulation with growing plant for X components
	- cyl3plant: defines step of the fixed point iteration
	- scenario_setup: defines soil and platn characteristics
	- PhloemPhotosynthesis: defines plant water and carbon flow
	- printData: all functions to print to files
	- helpfull: other usefull functions
	- notused: funcitons kept just in case
- modules: 
	- python wrappers around dumux classes (solverbase, richars, richars_no_mpi)
	- rhizo_modelsPlant: class containing all the 1d models !!!ATT this scripts was not cleaned at all!!!
	- air_modelsPlant: dummy class for 1d models of shoot segments and root segments aboveground (no perirhizal zones)
	- vtk_plot_adapted: slightly adapted vtk_plot file that I still need to merge with the master branch

# launch a simulation
A simulation can be launch by going in the scripts folder and running:
`python3 XcGrowth.py startSim endSim`

# remarks 
- 'cyl' = 'cylinders' = '1d models' = perirhizal zone around the root
- bulk soil = 3d model
- scv: sub-control volume (volume of a cell of the grid)
- by setting a 0 concentration of microbes, we can go back to the 'richards2c' simulations
- we have a dynamic (=/= instantaneous) adsorption. It is necessary for the 1d-3d coupling, otherwise we get different adsorption values between the bulk
soil and the 1d models
- we set an empirical minimum number of loops (4). During our simulation, we usually saw a significant
decrease of the errors during the first 4 loops. Might need to adapt that for another setup
- currently plant-soil solute flow defined outside of iteration loop. will be different for N or P simulation

- it works with the CPB master branch
- currently, cpp printing when creating a 1d model is redirected to the file 'stdcout_cpp.txt', if there are a lot of 1d model being created, it might be better to just catch the sdt::cout call without printing it afterwards to no end up with a very long file.



# TODOs
- clean the file rhizomodels_Plant
- check/adapt the file to make it work for growing root systems
- implement solute uptake (N or P)
- add unit information
- separate phloem and plant water flow in different files, adapt transpiration
- make messages print to consol more readable
- merge fpit_clean with master branch, and merge file in fixedPointIter/modules with the default ones of the master branch

