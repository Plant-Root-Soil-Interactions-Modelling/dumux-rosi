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



# TODOs
- implement solute uptake (N or P)
- finish adding unit information + funciton description
- separate phloem and plant water flow in different files, adapt transpiration
- make messages print to consol more readable
- merge fpit_clean with master branch, and merge file in fixedPointIter/modules with the default ones of the master branch
- create test units for setup
- limit in dumux the sink/source term (like for the bcs)
- [.. to be completed by users..]


# Other
## implementation tillage:
- because of the soil van genuchten parameters, a small change in water content leads to big change in soil watre potential.
 might level of precision in python might become an issue? 
 level precision 1e-14 seems enough > max level of precision python (15–17)
 pressure_head(theta_wilting_point+1e-13, vg_soil) >  pressure_head(theta_wilting_point , vg_soil) (more 1cm diff)
 pressure_head(theta_wilting_point+1e-14, vg_soil) ~  pressure_head(theta_wilting_point , vg_soil)
 pressure_head(theta_wilting_point-1e-14, vg_soil) ~  pressure_head(theta_wilting_point , vg_soil) (less 1cm diff)
 also the higher rounding threshold in the setup (between 1e-11..1e-20) might need to be adapted for water
- for the same reason, we might need to set a lower max diff between 1d and 3d and/or finally set up the manual correction of water content after each time step to avoid divergence. How to make threshold of 1d-3d difference more adapted/dependen on VG params?
(psi + theta threshold?)
- retry implementing the 1d1d flow
- even for high time step and transliration rates, it manages to converge BUT still cannot respect respect the PWU.
Do we switch to dirichlet BC when needed? y do psi remain so high when pwu should take all it has
- always convergence stoipes after ~ 19s. but issue starts earlier
- att python3 XcGrowth.py 0.5 10 0 none could go to the next time step even though PWU_error was very high
- python3 XcGrowth.py 0.5 10 0 none with dt = 20mn still fails when numcells ==1
- python3 XcGrowth.py 9.5 10 0 none with dt = 20mn still fails when numcells ==1 (no 3d-3d flow or 1d3d flow)
- high trans is ok for big root system it seems
- had to put lower maxrelshift (1e-12 instreaad of 1e-8)
- higher dx makes it harder to solve but lower dx does not really help either
- dxmin overwise not enought water
- controle the splitsoil val as well for water content (not done already?)
- for splitsoilval, need to do reset and then get also as input the futur plant water uptake
- check that splitSoilval works even if a thetaCyl is < 0 (i think already ok)

