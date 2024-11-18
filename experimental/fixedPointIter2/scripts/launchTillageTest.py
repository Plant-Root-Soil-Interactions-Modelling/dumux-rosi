import os
import sys
#os.chdir('/home/rbtlm640/dumux10c38/dumux/dumux-rosi/experimental/fixedPointIter2/scripts')
sys.path.append('')
import mainTillage
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

doProfile = True
if doProfile:
    import cProfile
    import pstats, io
    pr = cProfile.Profile()
    pr.enable()
spellData = {'scenario': 'none', 'spellStart': 92, 'spellEnd': 93,'condition': 'wet'}
results_dir = mainTillage.XcGrowth(25.5, 25.51,0,spellData ,doProfile=doProfile)
if doProfile:
    pr.disable()
    filename = results_dir+'profile'+str(rank)+'.prof' 
    print(filename)
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
    ps.dump_stats(filename)


