#!/home/h.mai/share/tools/anaconda2/bin/python2
from HF5Data import HF5Data
import numpy as np
import scipy.stats as stats
import os
import sys
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plotSoil(inputName, variableName, unit):
	import numpy
	from matplotlib import pyplot as plt
	from scipy import interpolate
	from scipy import integrate
	from scipy import interp
	import pandas as pd
	import seaborn as sns
	import csv

	###
	#SETTING COLOR http://www.discoveryplayground.com/computer-programming-for-kids/rgb-colors/
	###
	cblue='#4169e1'
	cgreen='#228b22'
	cred='#ff0000'
	corange='#ff8c00'
	cpurple='#9370db'
	cpink='#ff69b4'
	cbrown='#a0522d'
	cgray='#696969'
	cyellow='#ffd700'


	######### LOADIND RESULTS
	NumberOfLayers = 2
	with open(inputName, 'r') as csvfile:
	    reader = csv.reader(csvfile)
	    temp_ = [[float(e) for e in r] for r in reader]
	time = []
	totalsoilValues = []
	topSoilsoilValues = []
	subSoilTypesoilValues = []
	for r in temp_:
		time.append(r[0])
		totalsoilValues.append(r[int(NumberOfLayers+1)])
		topSoilsoilValues.append(r[1])
		subSoilTypesoilValues.append(r[2])

	timeRef=numpy.asarray(time)
	f_soilValues = numpy.asarray(totalsoilValues)#*3600*24
	f_soilValuesTopSoil = numpy.asarray(topSoilsoilValues)#*3600*24
	f_soilValuesSubSoil = numpy.asarray(subSoilTypesoilValues)#*3600*24

	###################
	####
	###################
	#plot
	sns.set()
	#http://seaborn.pydata.org/tutorial/aesthetics.html
	ax = plt.subplot(111, xlabel='days', ylabel=unit)
	ax.set_title(variableName+'\n'+inputName,weight='bold')
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
	             ax.get_xticklabels() + ax.get_yticklabels()):
	    item.set_fontsize(12)

	#plt.plot(timeRef/24/3600, f_soilValues,  color=cred,ls="-",linewidth=2)
	plt.plot(timeRef/24/3600, f_soilValuesTopSoil, color= corange,ls="--",linewidth=2)
	plt.plot(timeRef/24/3600, f_soilValuesSubSoil, color= cblue,ls="--",linewidth=2,marker="o", markevery=0.3)

	plt.legend(["top", "bottom"], loc='lower left')
	#x1,x2,y1,y2 = plt.axis()
	#plt.axis((x1,20,y1,y2))
	plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	plt.savefig(inputName+'_'+variableName+'TopBottom.png', bbox_inches='tight', dpi=300)
	#sns.despine()
	#plt.show()

def filterHF5(H5_filename):
	SELF=HF5Data(H5_filename,"a")
	SELF.scanFileForInfosAndPositions()

#select variables
	selectedVarNums = [0,1]
	print ("")
	print ("YOU'VE SELECTED THESE VARIABLES:")
	for i in range(len(selectedVarNums)):
		print ("                                    ", SELF.varnames[selectedVarNums[i]])
	print ("")

	selectedTimeSteps = []
	for i in range(len(SELF.timeSteps)):
		selectedTimeSteps.append(i)

	print (selectedVarNums)

	NumberOfLayers = 2
	soilXYZ = SELF.readVarAtTime(1,3)

	path = H5_filename[0:H5_filename.index(".")]
	if not os.path.exists(path):
		os.mkdir( path, 0o755 );

	for j in selectedVarNums:
		soilValues=[]
		for itimeSel in (selectedTimeSteps):
			sys.stdout.write ("\ranalysing at time step: "+str(itimeSel)+'/'+str(len(selectedTimeSteps))+' .....')
			if SELF.vartype[j] == 0:
				N=SELF.nPoin
			else:
				N=SELF.elem
			#print N
			varPos=j
			value=SELF.readVarAtTime(itimeSel,varPos)
			soilValues_=np.zeros(NumberOfLayers+2)
			soilValues_[0]=SELF.timeSteps[itimeSel] # time
			volTopSoil = 0
			volSubSoil = 0
			for k in xrange (N*SELF.varcomp[j]):
				#print (soilXYZ[k*3+2], voxelVol[k])
				if (soilXYZ[k*3+2]>=-0.00):
					soilValues_[1] += value[k] # in topsoil
					volTopSoil += 1
				else:
					if (soilXYZ[k*3+2]<=-0.55):
						soilValues_[2] += value[k] # in Subsoil
						volSubSoil += 1
				#soilValues_[int(NumberOfLayers+1)] += value[k]*voxelVol[k] #total

			soilValues_[1] /= volTopSoil
			soilValues_[2] /= volSubSoil
			#soilValues_[3] /= (volTopSoil + volSubSoil)
			#print volTopSoil, volSubSoil
			soilValues.append(soilValues_)

		print ('saving analysis results and plotting..')
		# Write to file
		varName = ''
		unit = ''
		if (j==0):
			varName = "waterPressure"
			unit = "Pa"
		else:
			varName = "waterContent"
			unit = "[-]"
		outputName = ("./"+path+"/"+H5_filename+'_'+varName+'_byLayers.csv')
		#outputName = ('test.csv')
		with open(outputName, 'w') as csvfile:
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in soilValues]
		plotSoil(outputName, varName, unit)

		if (j==0):
			varName = "waterPotential"
			unit = "MPa"
			outputName = ("./"+path+"/"+H5_filename+'_'+varName+'_byLayers.csv')
			#outputName = ('test.csv')
			with open(outputName, 'w') as csvfile:
			    writer = csv.writer(csvfile)
			    for r in soilValues:
					r[1] -=1e5
					r[1] /=1e6
					r[2] -=1e5
					r[2] /=1e6
					r[3] -=1e5
					r[3] /=1e6
					writer.writerow(r)
			plotSoil(outputName, varName, unit)

if __name__=='__main__':
    sys.exit(filterHF5(sys.argv[1]))
