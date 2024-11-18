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
- for 1ds and 3ds computation: manually adapt sink before hand to never send or take more than they can deal with
- for splitsoilval and set sink3d: look at "potentially available water": water beginning of time step + [a] pwu for splitsoilval, [b] max(Q3d3d, 0.) for set sink 3d
- resolution: higher dx makes it harder to solve but lower dx does not really help either. dxmin overwise not enought water.
- ATT: where the old dumux threw a non-convergence error, this one can give back bad (eg. negative values) better do more check on the simulation outputs.


# TODOs
- merge exudation_analysis_plant3 with master
- finish adding unit information + funciton description
- merge fpit_clean with master branch, and merge file in fixedPointIter/modules with the default ones of the master branch
- create test units for setup
- retry implementing the 1d1d flow
- implement analytic dumux solution for richardsnc problems
- setup the other solvers for dumux10c
- do MPI with shared memory? could make data exchange quicker
- [.. to be completed by users..]

# Other
## implementation TraiRhizo
- upwindweight 0.5 or 1?
- reset the update kr and kx for cplant box (use the adapted cplant box.

## implementation tillage:
- use the mappedplant class. check that the obtained plant is the same as when using the rootsystem class.
- set eps for vangenuchten functinos to 1e-10. see if creates issues.
- because of the soil van genuchten parameters, a small change in water content leads to big change in soil watre potential.
 might level of precision in python might become an issue? 
 level precision 1e-14 seems enough > max level of precision python (15–17)
 pressure_head(theta_wilting_point+1e-13, vg_soil) >  pressure_head(theta_wilting_point , vg_soil) (more 1cm diff)
 pressure_head(theta_wilting_point+1e-14, vg_soil) ~  pressure_head(theta_wilting_point , vg_soil)
 pressure_head(theta_wilting_point-1e-14, vg_soil) ~  pressure_head(theta_wilting_point , vg_soil) (less 1cm diff)
 also the higher rounding threshold in the setup (between 1e-11..1e-20) might need to be adapted for water
- for the same reason, we might need to set a lower max diff between 1d and 3d and/or finally set up the manual correction of water content after each time step to avoid divergence. How to make threshold of 1d-3d difference more adapted/dependen on VG params?
(psi + theta threshold?)
- had to put lower maxrelshift (1e-12 instead of 1e-8)
- controle the splitsoil val as well for water content (not done already?)
- check that splitSoilval works even if a thetaCyl is < 0 (i think already ok) and when pwu > cell water content
issue convergence.
possible solutions:
3) use the mean of the iterations => seems to work. we do get input rsi  ~ output rsi although the most important is that [a] PWU ~ PWUlim ~ 0 , [b] convergence. then ok if rsi_input != rsi_real (maybe check that rsi_input in [rsi_init, rsi_real]. __Need much more iteration__
7) maybe try and get the integrated kr directly from dumux to see what is going on <= almost like doing a smaller inner loop betweer
plant and 1ds


==> test at 25d, for lower wat pot, in parallel. adapt for higher time steps, use the sub-stepping for plant-1ds flow
so there is a small imprecision regarding the soil theta with the vg_soil [0.08, 0.43, 0.04, 1.6, 50]
error negligible for normal soil pset. todo: fix 