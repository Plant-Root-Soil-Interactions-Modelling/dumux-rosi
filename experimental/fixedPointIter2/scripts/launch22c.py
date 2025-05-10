import os
import sys
#os.chdir('/home/rbtlm640/dumux10c38/dumux/dumux-rosi/experimental/fixedPointIter2/scripts')
sys.path.append('')
import main22c

spellData = {'scenario': 'none', 'spellStart': 30, 'spellEnd': 35,'condition': 'wet'}
results_dir = main22c.XcGrowth(5., 25.,44,spellData )
# (initsim, simMax,paramIndx_,spellData)

