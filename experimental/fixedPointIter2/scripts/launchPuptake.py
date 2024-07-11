import os
import sys
#os.chdir('/home/rbtlm640/dumux10c38/dumux/dumux-rosi/experimental/fixedPointIter2/scripts')
sys.path.append('')
import mainPuptake

spellData = {'scenario': 'none', 'spellStart': 92, 'spellEnd': 93,'condition': 'wet'}
results_dir = mainPuptake.XcGrowth(8, 30.,0,spellData )


