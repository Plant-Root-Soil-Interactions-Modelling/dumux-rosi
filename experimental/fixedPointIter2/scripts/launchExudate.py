import os
import sys
os.chdir('/home/rbtl2404/dumuxMagda2/dumux/dumux-rosi/experimental/fixedPointIter2/scripts')
sys.path.append('')
import mainExudate

# spellData = {'scenario': 'none', 'spellStart': 30, 'spellEnd': 35,'soil_type': 'loam'}
scenarioData = {'soil_type': 'loam'}
results_dir = mainExudate.XcGrowth(100, scenarioData )


