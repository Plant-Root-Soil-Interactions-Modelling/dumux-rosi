import os
import sys
os.chdir('/home/rbtlm640/dumux10c38/dumux/dumux-rosi/experimental/fixedPointIter2/scripts')
sys.path.append('')
import mainTillage

spellData = {'scenario': 'none', 'spellStart': 12, 'spellEnd': 13,'condition': 'wet'}
results_dir = mainTillage.XcGrowth(9.5, 10.,0,spellData )


