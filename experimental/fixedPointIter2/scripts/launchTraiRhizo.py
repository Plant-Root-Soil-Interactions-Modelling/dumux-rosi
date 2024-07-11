import os
import sys
os.chdir('/home/rbtlm640/dumux10c38/dumux/dumux-rosi/experimental/fixedPointIter2/scripts')
sys.path.append('')
import mainTraiRhizo

spellData = {'scenario': 'baseline', 'spellStart': 11, 'spellEnd': 18,'condition': 'wet'}
results_dir = mainTraiRhizo.XcGrowth(10., 25.,44,spellData )
# (initsim, simMax,paramIndx_,spellData)

