import os
import sys

sys.path.append('')
import mainExudate

scenarioData = {'soil_type': 'loam', 'res' : '4', 'simMax' : '60'}
results_dir = mainExudate.XcGrowth(scenarioData)


