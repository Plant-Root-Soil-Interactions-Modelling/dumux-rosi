#!/usr/bin/env python
import os
import subprocess
import multiprocessing as mp

if (not os.path.exists("resultsN4_R")):
	os.mkdir("resultsN4_R");

# the scenarios
scenarios = [
["-ParameterFile", "rootSystemRiceFluctuationV3_.input", "-Grid.File", "./grids/Rice_NERICA4_PlusP_FC_Mai_2019_V51.dgf", "-Problem.Name", "./resultsN4_R/N4_PlusP_FC_FlucV2_R", "-Soil.SpatialParams.SoilDataFileName", "BC_soilColumnSubP_V2.csv" ,"-Soil.BoundaryConditions.WeatherDataFileName", "BC3_N4_SubP_Drying.csv"],
["-ParameterFile", "rootSystemRiceFluctuationV3_.input", "-Grid.File", "./grids/Rice_NERICA4_PlusP_Drying_Mai_2019_V41.dgf", "-Problem.Name", "./resultsN4_R/N4_PlusP_Drying_FlucV2_R", "-Soil.SpatialParams.SoilDataFileName", "BC_soilColumnSubP_V2.csv", "-Soil.BoundaryConditions.WeatherDataFileName", "BC3_N4_SubP_Drying.csv"],
["-ParameterFile", "rootSystemRiceFluctuationV3_.input", "-Grid.File", "./grids/Rice_NERICA4_SubP_FC_Mai_2019V41.dgf", "-Problem.Name", "./resultsN4_R/N4_SubP_FC_FlucV2_R", "-Soil.SpatialParams.SoilDataFileName", "BC_soilColumnSubP_V2.csv", "-Soil.BoundaryConditions.WeatherDataFileName", "BC3_N4_SubP_Drying.csv"],
["-ParameterFile", "rootSystemRiceFluctuationV3_.input", "-Grid.File", "./grids/Rice_NERICA4_SubP_Drying_Mai_2019V51.dgf", "-Problem.Name", "./resultsN4_R/N4_SubP_Drying_FlucV2_R", "-Soil.SpatialParams.SoilDataFileName", "BC_soilColumnSubP_V2.csv", "-Soil.BoundaryConditions.WeatherDataFileName", "BC3_N4_SubP_Drying.csv"],
["-ParameterFile", "rootSystemRiceFluctuationV3_.input", "-Grid.File", "./grids/Rice_NERICA4_NoP_FC_Mai_2019V51.dgf", "-Problem.Name", "./resultsN4_R/N4_NoP_FC_FlucV2_R", "-Soil.SpatialParams.SoilDataFileName", "BC_soilColumnSubP_V2.csv", "-Soil.BoundaryConditions.WeatherDataFileName", "BC3_N4_SubP_Drying.csv"],
["-ParameterFile", "rootSystemRiceFluctuationV3_.input", "-Grid.File", "./grids/Rice_NERICA4_NoP_Drying_Mai_2019V51.dgf", "-Problem.Name", "./resultsN4_R/N4_NoP_Drying_FlucV2_R", "-Soil.SpatialParams.SoilDataFileName", "BC_soilColumnSubP_V2.csv", "-Soil.BoundaryConditions.WeatherDataFileName", "BC3_N4_SubP_Drying.csv"],
["-ParameterFile", "rootSystemRiceFluctuationV3_.input", "-Grid.File", "./grids/Rice_NERICA4_PlusP_FC_Mai_2019_V51.dgf", "-Problem.Name", "./resultsN4_R/N4_PlusP_FC_FlucV2", "-Soil.SpatialParams.SoilDataFileName", "BC_soilColumnPlusP_V2.csv" ,"-Soil.BoundaryConditions.WeatherDataFileName", "BC3_N4_PlusP_FC.csv"],
]

def run_scenario(args):
	subprocess.call(["./rosiHybridUG2c"] + args)

# run the scenarios using multiple processes
#print ("Computing data on {} processes.".format(mp.cpu_count()))
with mp.Pool() as pool:
        pool.map(run_scenario, scenarios)



