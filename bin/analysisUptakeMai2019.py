#!/home/h.mai/share/tools/anaconda2/bin/python2
from HF5Data import HF5Data
import numpy as np
#import scipy.stats as stats
import os
import sys
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate

import psutil
import gc

def plotRootCollarWaterPotential(inputName):
	import numpy
	from matplotlib import pyplot as plt
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

	rootPressure = np.genfromtxt(inputName, delimiter=',')
	timeRef = rootPressure[:,0]#/24/3600
	#hourlyTimeStep = numpy.asarray(range(int(timeRef[0]), int(timeRef[-1]), 3600))
	#f_rootPressure = interp(timeRef, rootPressure[:,0], (rootPressure[:,1]-1e5)/(1000*9.81))
	f_rootPressure = (rootPressure[:,1]-1e5)/1e6
	ax = plt.subplot(111, xlabel='days', ylabel='MPa ', title='Root water potential')
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
	             ax.get_xticklabels() + ax.get_yticklabels()):
		item.set_fontsize(12)
	plt.plot(rootPressure[:,0]/86400., f_rootPressure, ls="--",linewidth=2)
	plt.savefig(inputName+'V3.png', bbox_inches='tight', dpi=300)

def plotRootGrowth(inputName, variableName, unit, factor):
	import numpy
	from matplotlib import pyplot as plt
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

	NumberOfRoots = 3
	NumberOfLayers = 2

	values = np.genfromtxt(inputName, delimiter=',')
	timeRef = values[:,0]/24/3600
	#print (timeRef)
	total = values[:,NumberOfRoots*NumberOfLayers+NumberOfRoots+NumberOfLayers+1]*factor
	topSoilTotal = values[:,NumberOfRoots*NumberOfLayers+NumberOfRoots+1]*factor
	subSoilTotal = values[:,NumberOfRoots*NumberOfLayers+NumberOfRoots+2]*factor
	nodalTotal = values[:,NumberOfRoots*NumberOfLayers+1]*factor
	STypeTotal = values[:,NumberOfRoots*NumberOfLayers+2]*factor
	LTypeTotal = values[:,NumberOfRoots*NumberOfLayers+3]*factor
	nodalTopSoil = values[:,1]*factor
	STypeTopSoil = values[:,2]*factor
	LTypeTopSoil = values[:,3]*factor
	nodalSubSoil = values[:,NumberOfRoots+1]*factor
	STypeSubSoil = values[:,NumberOfRoots+2]*factor
	LTypeSubSoil = values[:,NumberOfRoots+3]*factor
	###################
	####
	###################
	#plot RATE
	sns.set()
	#http://seaborn.pydata.org/tutorial/aesthetics.html
	ax = plt.subplot(111, xlabel='days', ylabel=unit)
	ax.set_title(variableName,weight='bold')
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
	             ax.get_xticklabels() + ax.get_yticklabels()):
	    item.set_fontsize(12)

	plt.plot(timeRef, total,  color=cred,ls="-",linewidth=2, alpha=0.5)
	plt.plot(timeRef, nodalTotal, color= cgreen,ls="--",linewidth=2) #blue
	plt.plot(timeRef, STypeTotal, color= corange,ls="--",linewidth=2)
	plt.plot(timeRef, LTypeTotal, color= cblue,ls="--",linewidth=2)

	plt.legend(["total","nodals", "S Type", "L type"], loc='upper left')
	#x1,x2,y1,y2 = plt.axis()
	#plt.axis((x1,x2,y1,fixedY))
	plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	plt.show()
	plt.savefig(inputName+'_'+variableName+'RootTypeV3.png', bbox_inches='tight', dpi=300)
	#sns.despine()
	plt.clf()

	#plot RATE for layers
	sns.set()
	#http://seaborn.pydata.org/tutorial/aesthetics.html
	ax = plt.subplot(111, xlabel='days', ylabel=unit)
	ax.set_title(variableName,weight='bold')
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
	             ax.get_xticklabels() + ax.get_yticklabels()):
	    item.set_fontsize(12)

	plt.plot(timeRef, total,  color=cred,ls="-",linewidth=2, alpha=0.5)
	plt.plot(timeRef, topSoilTotal, color= corange,ls="--",linewidth=2)
	plt.plot(timeRef, subSoilTotal, color= cblue,ls="--",linewidth=2)

	plt.legend(["total","topSoil", "subSoil"], loc='upper left')
	#x1,x2,y1,y2 = plt.axis()
	#plt.axis((x1,x2,y1,fixedY))
	plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	plt.savefig(inputName+'_'+variableName+'SoilLayerV3.png', bbox_inches='tight', dpi=300)
	#sns.despine()
	plt.show()
	plt.clf()

#	#plot  for each type in each layers
#	sns.set()
#	#http://seaborn.pydata.org/tutorial/aesthetics.html
#	ax = plt.subplot(111, xlabel='days', ylabel=unit)
#	ax.set_title(variableName,weight='bold')
#	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#	             ax.get_xticklabels() + ax.get_yticklabels()):
#	    item.set_fontsize(12)
#
#	plt.plot(timeRef, total,  color=cred,ls="-",linewidth=2, alpha=0.5)
#	plt.plot(timeRef, nodalTopSoil, color= cgreen,ls="--",linewidth=2,marker="v", markevery=0.3) #blue
#	plt.plot(timeRef, STypeTopSoil, color= corange,ls="--",marker="v",linewidth=2, markevery=0.3)
#	plt.plot(timeRef, LTypeTopSoil, color= cblue,ls="--",linewidth=2,marker="v", markevery=0.3)
#
#	plt.plot(timeRef, nodalSubSoil, color= cgreen,ls="-",linewidth=2, alpha=0.5) #blue
#	plt.plot(timeRef, STypeSubSoil, color= corange,ls="-",linewidth=2, alpha=0.5)
#	plt.plot(timeRef, LTypeSubSoil, color= cblue,ls="-",linewidth=2, alpha=0.5)
#
#	plt.legend(["total","nodals TopSoil", "S Type TopSoil", "L type TopSoil","nodals SubSoil", "S Type SubSoil", "L type SubSoil"], loc='upper left')
#	#x1,x2,y1,y2 = plt.axis()
#	#plt.axis((x1,x2,y1,fixedY))
#	plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
#	plt.savefig(inputName+'_'+variableName+'RootTypeSoilLayerV3.png', bbox_inches='tight', dpi=300)
#	#sns.despine()
#	plt.show()
#	plt.clf()

def Analysing(H5_filename):

	process = psutil.Process(os.getpid())
	print(process.memory_info().rss)

	SELF=HF5Data(H5_filename,"a")
	SELF.scanFileForInfosAndPositions()

#select variables

	selectedVarNums = []
	print ("")
	print (u"Seclect Variables:")
	print ("  No     Variable  ")
	for i in range(len(SELF.varnames)):
		print (" ", i, "    ", SELF.varnames[i], "   ")
	print ("")

	waterPressureIdx = 0
	radiusIdx = 1
	lengthIdx = 2
	oderIdx = 3
	branchIdx = 4
	volumeFlowIdx = 6
	uptakeRateWaterIdx = 7
	uptakeRatePIdx = 8
	ageIdx = 9
	distanceFromOriginIdx = 10
	coordIdx = 11
	concentrationIdx = 5
	selectedVarNums = [0]#radius 0, length 1, Oder 2, Age 5, Coord 6

	selectedTimeSteps = []
	for i in range(len(SELF.timeSteps)):
		selectedTimeSteps.append(i)
	print (selectedVarNums)

	### Reading static geometrical properties in each root segments
	rootOrder = SELF.readVarAtTime(1, oderIdx)
	rootRadius = SELF.readVarAtTime(1, radiusIdx)
	rootLength = SELF.readVarAtTime(1, lengthIdx)
	rootBranch = SELF.readVarAtTime(1, branchIdx)
	rootXYZ=SELF.readVarAtTime(1, coordIdx)
	distanceFromOrigin = SELF.readVarAtTime(1, distanceFromOriginIdx)

	### Define the max branchIdx
	maxBranchIdx = 0
	for k in xrange (SELF.elem):
		if (rootBranch[k] > maxBranchIdx):
			maxBranchIdx = int(rootBranch[k])

	NumberOfRoots = 3
	NumberOfLayers = 2
	NumberOfValue = NumberOfRoots*NumberOfLayers+NumberOfRoots+NumberOfLayers
	for j in selectedVarNums:
		rootSurface=[]
		rootVolume=[]
		rootLengthOutput=[]
		rootTipNumber=[]
		uptakeRateWater=[]
		uptakeRateP=[]
		waterCollarPressure = []
		waterMassFlow = []
		concentrationFlow = []
		branchLength = np.zeros(maxBranchIdx+1)
		rootTipElemIdx = np.zeros(maxBranchIdx+1)
		for itimeSel in (selectedTimeSteps):
			sys.stdout.write ("\ranalysing at time step: "+str(itimeSel)+'/'+str(len(selectedTimeSteps))+' .....')
			if SELF.vartype[j] == 0:
				N=SELF.nPoin
			else:
				N=SELF.elem
			varPos=j

			rootAge=SELF.readVarAtTime(itimeSel,ageIdx)
			for k in xrange (SELF.elem):
				if (rootAge[k] > 0):
					#branchLengthIdx = branchIdxSet.index(rootBranch[k])
					branchIdx = int(rootBranch[k])
					if (distanceFromOrigin[k] > branchLength[branchIdx]):
						branchLength[branchIdx] = distanceFromOrigin[k]
						rootTipElemIdx[branchIdx] = k

			rootSurface_=np.zeros(NumberOfValue+2)
			rootSurface_[0]=SELF.timeSteps[itimeSel] # time
			rootVolume_=np.zeros(NumberOfValue+2)
			rootVolume_[0]=SELF.timeSteps[itimeSel] # time
			rootLength_=np.zeros(NumberOfValue+2)
			rootLength_[0]=SELF.timeSteps[itimeSel] # time
			rootTipNumber_=np.zeros(NumberOfValue+2)
			rootTipNumber_[0]=SELF.timeSteps[itimeSel] # time

			rootWaterUptakeRate=SELF.readVarAtTime(itimeSel,uptakeRateWaterIdx)
			uptakeRateWater_=np.zeros(NumberOfValue+2)
			uptakeRateWater_[0]=SELF.timeSteps[itimeSel] # time

			rootPUptakeRate=SELF.readVarAtTime(itimeSel,uptakeRatePIdx)
			uptakeRateP_=np.zeros(NumberOfValue+2)
			uptakeRateP_[0]=SELF.timeSteps[itimeSel] # time

			waterPressure=SELF.readVarAtTime(itimeSel,waterPressureIdx)
			waterCollarPressure_=np.zeros(2)
			waterCollarPressure_[0]=SELF.timeSteps[itimeSel]
			waterCollarPressure_[1]=waterPressure[0]

			waterVolumeFlowAtTime=SELF.readVarAtTime(itimeSel,volumeFlowIdx)
			waterMassFlowAtTime = waterVolumeFlowAtTime*1000
			waterMassFlow_=np.zeros(NumberOfValue+2)
			waterMassFlow_[0]=SELF.timeSteps[itimeSel]

			concentrationAtTime=SELF.readVarAtTime(itimeSel,concentrationIdx)
			concentrationFlowAtTime = concentrationAtTime * waterVolumeFlowAtTime
			concentrationFlow_=np.zeros(NumberOfValue+2)
			concentrationFlow_[0]=SELF.timeSteps[itimeSel]

			for k in xrange (N*SELF.varcomp[j]):
				branchIdx = int(rootBranch[k])
				rootTypeIdx = 1
				if (rootOrder[k] == 5):
					rootTypeIdx = 2
				if (rootOrder[k] == 2 or rootOrder[k] == 3):
					rootTypeIdx = 3
				if (rootAge[k]>0):
					if (rootXYZ[k*3+2]>=-0.21):
						rootSurface_[rootTypeIdx] += 2*3.14*rootRadius[k]*rootLength[k] #root surface per root types
						rootVolume_[rootTypeIdx] += 3.14*rootRadius[k]*rootRadius[k]*rootLength[k] #root volume per root types
						rootLength_[rootTypeIdx] += rootLength[k] #root length per root types

						rootSurface_[NumberOfRoots*NumberOfLayers+NumberOfRoots+1] += 2*3.14*rootRadius[k]*rootLength[k] #total in topsoil
						rootVolume_[NumberOfRoots*NumberOfLayers+NumberOfRoots+1] += 3.14*rootRadius[k]*rootRadius[k]*rootLength[k] #total in topsoil
						rootLength_[NumberOfRoots*NumberOfLayers+NumberOfRoots+1] += rootLength[k] #total in topsoil

						uptakeRateWater_[rootTypeIdx] += rootWaterUptakeRate[k] #uptake per root types
						uptakeRateP_[rootTypeIdx] += rootPUptakeRate[k]

						uptakeRateWater_[NumberOfRoots*NumberOfLayers+NumberOfRoots+1] += rootWaterUptakeRate[k]
						uptakeRateP_[NumberOfRoots*NumberOfLayers+NumberOfRoots+1] += rootPUptakeRate[k]
						waterMassFlow_[rootTypeIdx] += waterMassFlowAtTime[k]
						waterMassFlow_[NumberOfRoots*NumberOfLayers+NumberOfRoots+1] += waterMassFlowAtTime[k]

						concentrationFlow_[rootTypeIdx] += concentrationFlowAtTime[k]
						concentrationFlow_[NumberOfRoots*NumberOfLayers+NumberOfRoots+1] += concentrationFlowAtTime[k]

						#if (rootAge[k]<86400):
						if (rootTipElemIdx[branchIdx] == k):
							rootTipNumber_[rootTypeIdx] += 1 #root tips number per root types
							rootTipNumber_[NumberOfRoots*NumberOfLayers+NumberOfRoots+1] += 1 #total in topsoil
					else:
						rootSurface_[int(NumberOfRoots+rootTypeIdx)] += 2*3.14*rootRadius[k]*rootLength[k] #root surface per root types
						rootVolume_[int(NumberOfRoots+rootTypeIdx)] += 3.14*rootRadius[k]*rootRadius[k]*rootLength[k] #root volume per root types
						rootLength_[int(NumberOfRoots+rootTypeIdx)] += rootLength[k] #root length per root types

						rootSurface_[NumberOfRoots*NumberOfLayers+NumberOfRoots+2] += 2*3.14*rootRadius[k]*rootLength[k] #total in Subsoil
						rootVolume_[NumberOfRoots*NumberOfLayers+NumberOfRoots+2] += 3.14*rootRadius[k]*rootRadius[k]*rootLength[k] #total in Subsoil
						rootLength_[NumberOfRoots*NumberOfLayers+NumberOfRoots+2] += rootLength[k] #total in Subsoil

						uptakeRateWater_[int(NumberOfRoots+rootTypeIdx)] += rootWaterUptakeRate[k] #uptake per root types
						uptakeRateP_[int(NumberOfRoots+rootTypeIdx)] += rootPUptakeRate[k]
						uptakeRateP_[int(NumberOfRoots+rootTypeIdx)] += rootPUptakeRate[k]

						uptakeRateWater_[NumberOfRoots*NumberOfLayers+NumberOfRoots+2] += rootWaterUptakeRate[k]
						uptakeRateP_[NumberOfRoots*NumberOfLayers+NumberOfRoots+2] += rootPUptakeRate[k]

						waterMassFlow_[int(NumberOfRoots+rootTypeIdx)] += waterMassFlowAtTime[k]
						waterMassFlow_[NumberOfRoots*NumberOfLayers+NumberOfRoots+2] += waterMassFlowAtTime[k]

						concentrationFlow_[int(NumberOfRoots+rootTypeIdx)] += concentrationFlowAtTime[k]
						concentrationFlow_[NumberOfRoots*NumberOfLayers+NumberOfRoots+2] += concentrationFlowAtTime[k]

						#if (rootAge[k]<86400):
						if (rootTipElemIdx[branchIdx] == k):
							rootTipNumber_[int(NumberOfRoots+rootTypeIdx)] += 1 #root tip number per root types
							rootTipNumber_[NumberOfRoots*NumberOfLayers+NumberOfRoots+2] += 1 #total in Subsoil

					rootSurface_[int(NumberOfRoots*NumberOfLayers+rootTypeIdx)] += 2*3.14*rootRadius[k]*rootLength[k] #root surface per root types
					rootVolume_[int(NumberOfRoots*NumberOfLayers+rootTypeIdx)] += 3.14*rootRadius[k]*rootRadius[k]*rootLength[k] #root volume per root types
					rootLength_[int(NumberOfRoots*NumberOfLayers+rootTypeIdx)] += rootLength[k] #root length per root types

					rootSurface_[int(NumberOfValue+1)] += 2*3.14*rootRadius[k]*rootLength[k] #total root surface
					rootVolume_[int(NumberOfValue+1)] += 3.14*rootRadius[k]*rootRadius[k]*rootLength[k] #total root volume
					rootLength_[int(NumberOfValue+1)] += rootLength[k] #total root length

					uptakeRateWater_[int(NumberOfRoots*NumberOfLayers+rootTypeIdx)] += rootWaterUptakeRate[k]
					uptakeRateP_[int(NumberOfRoots*NumberOfLayers+rootTypeIdx)] += rootPUptakeRate[k]

					uptakeRateWater_[int(NumberOfValue+1)] += rootWaterUptakeRate[k]
					uptakeRateP_[int(NumberOfValue+1)] += rootPUptakeRate[k]

					waterMassFlow_[int(NumberOfValue+1)] += waterMassFlowAtTime[k]
					waterMassFlow_[int(NumberOfRoots*NumberOfLayers+rootTypeIdx)] += waterMassFlowAtTime[k]

					concentrationFlow_[int(NumberOfValue+1)] += concentrationFlowAtTime[k]
					concentrationFlow_[int(NumberOfRoots*NumberOfLayers+rootTypeIdx)] += concentrationFlowAtTime[k]

					#if (rootAge[k]<86400):
					if (rootTipElemIdx[branchIdx] == k):
						rootTipNumber_[int(NumberOfRoots*NumberOfLayers+rootTypeIdx)] += 1 #root length per root types
						rootTipNumber_[int(NumberOfValue+1)] += 1 #total root tip number

			rootSurface.append(rootSurface_) #total
			rootVolume.append(rootVolume_) #total
			rootLengthOutput.append(rootLength_) #total
			rootTipNumber.append(rootTipNumber_) #total

			uptakeRateWater.append(uptakeRateWater_)
			uptakeRateP.append(uptakeRateP_)
			waterMassFlow.append(waterMassFlow_)
			concentrationFlow.append(concentrationFlow_)

			waterCollarPressure.append(waterCollarPressure_)

		# Calculated cumulative uptake and its ratios
		rootSurface_ = np.asarray(rootSurface)
		rootVolume_ = np.asarray(rootVolume)
		rootLength_ = np.asarray(rootLengthOutput)
		rootTipNumber_ = np.asarray(rootTipNumber)
		uptakeRateWater_ = np.asarray(uptakeRateWater)
		uptakeRateP_ = np.asarray(uptakeRateP)

		waterMassFlow_ = np.asarray(waterMassFlow)
		concentrationFlow_ = np.asarray(concentrationFlow)

		time = uptakeRateWater_[:,0]
		#print(time.shape, uptakeRate.shape)
		c_uptakeRateWater = integrate.cumtrapz(uptakeRateWater_,time, initial=0, axis = 0)
		c_uptakeRateWater[:,0] = time

		c_uptakeRateP = integrate.cumtrapz(uptakeRateP_,time, initial=0, axis = 0)
		c_uptakeRateP[:,0] = time

		c_waterMassFlow = integrate.cumtrapz(waterMassFlow_,time, initial=0, axis = 0)
		c_waterMassFlow[:,0] = time

		c_concentrationFlow = integrate.cumtrapz(concentrationFlow_,time, initial=0, axis = 0)
		c_concentrationFlow[:,0] = time

		rootsurfacePerVolume = np.nan_to_num(rootSurface_/rootVolume_)
		rootsurfacePerVolume[:,0] = time

		POverWaterCum = np.nan_to_num(c_uptakeRateP/c_uptakeRateWater)
		POverWaterCum[POverWaterCum < 0] = 0
		POverWaterCum[:,0] = time

		PCumPerNumberRootTips = np.nan_to_num(c_uptakeRateP/rootTipNumber_)
		WaterCumPerNumberRootTips = np.nan_to_num(c_uptakeRateWater/rootTipNumber_)
		PCumPerNumberRootTips[:,0] = time
		WaterCumPerNumberRootTips[:,0] = time

		PRatePerVolume = np.nan_to_num(uptakeRateP_/rootVolume_)
		PRatePerSurface = np.nan_to_num(uptakeRateP_/rootSurface_)
		PRatePerLength = np.nan_to_num(uptakeRateP_/rootLength_)

		PCumPerVolume = np.nan_to_num(c_uptakeRateP/rootVolume_)
		PCumPerSurface = np.nan_to_num(c_uptakeRateP/rootSurface_)
		PCumPerLength = np.nan_to_num(c_uptakeRateP/rootLength_)

		WaterRatePerVolume = np.nan_to_num(uptakeRateWater_/rootVolume_)
		WaterRatePerSurface = np.nan_to_num(uptakeRateWater_/rootSurface_)
		WaterRatePerLength = np.nan_to_num(uptakeRateWater_/rootLength_)

		waterMassFlowPerVolume = np.nan_to_num(waterMassFlow_/rootVolume_)
		waterMassFlowPerSurface = np.nan_to_num(waterMassFlow_/rootSurface_)
		waterMassFlowPerLength = np.nan_to_num(waterMassFlow_/rootLength_)

		ConcentrationFlowPerVolume = np.nan_to_num(waterMassFlow_/rootVolume_)
		ConcentrationFlowPerSurface = np.nan_to_num(waterMassFlow_/rootSurface_)
		ConcentrationFlowPerLength = np.nan_to_num(waterMassFlow_/rootLength_)

		WaterCumPerVolume = np.nan_to_num(c_uptakeRateWater/rootVolume_)
		WaterCumPerSurface = np.nan_to_num(c_uptakeRateWater/rootSurface_)
		WaterCumPerLength = np.nan_to_num(c_uptakeRateWater/rootLength_)

		WaterTransportCumPerVolume = np.nan_to_num(c_waterMassFlow/rootVolume_)
		WaterTransportCumPerSurface = np.nan_to_num(c_waterMassFlow/rootSurface_)
		WaterTransportCumPerLength = np.nan_to_num(c_waterMassFlow/rootLength_)

		NutrientTransportCumPerVolume = np.nan_to_num(c_concentrationFlow/rootVolume_)
		NutrientTransportCumPerSurface = np.nan_to_num(c_concentrationFlow/rootSurface_)
		NutrientTransportCumPerLength = np.nan_to_num(c_concentrationFlow/rootLength_)

		ratioWaterUptakePerTransport = np.nan_to_num(uptakeRateWater_/waterMassFlow_)
		ratioNutrientUptakePerTransport = np.nan_to_num(uptakeRateP_/concentrationFlow_)

		PRatePerVolume[:,0] = time
		PRatePerSurface[:,0] = time
		PRatePerLength[:,0] = time
		PCumPerVolume[:,0] = time
		PCumPerSurface[:,0] = time
		PCumPerLength[:,0] = time
		WaterRatePerVolume[:,0] = time
		WaterRatePerSurface[:,0] = time
		WaterRatePerLength[:,0] = time
		WaterCumPerVolume[:,0] = time
		WaterCumPerSurface[:,0] = time
		WaterCumPerLength[:,0] = time

		waterMassFlowPerVolume[:,0] = time
		waterMassFlowPerSurface[:,0] = time
		waterMassFlowPerLength[:,0] = time
		WaterTransportCumPerVolume[:,0] = time
		WaterTransportCumPerSurface[:,0] = time
		WaterTransportCumPerLength[:,0] = time

		ConcentrationFlowPerVolume[:,0] = time
		ConcentrationFlowPerSurface[:,0] = time
		ConcentrationFlowPerLength[:,0] = time
		NutrientTransportCumPerVolume[:,0] = time
		NutrientTransportCumPerSurface[:,0] = time
		NutrientTransportCumPerLength[:,0] = time

		ratioWaterUptakePerTransport[:,0] = time
		ratioNutrientUptakePerTransport[:,0] = time

		# Write to file
		path = H5_filename[0:H5_filename.index(".")]
		if not os.path.exists(path):
			os.mkdir( path, 0o755 );

		print ('')
		print(process.memory_info().rss)
		print ('saving analysis results..')
		gc.collect()
		print(process.memory_info().rss)

		np.savetxt("./"+path+"/"+H5_filename+'_rootCollarPressureV3.csv', waterCollarPressure, delimiter=",")

		np.savetxt("./"+path+"/"+H5_filename+'_rootSurfaceV3.csv', rootSurface_, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_rootLengthV3.csv', rootLengthOutput, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_rootVolumeV3.csv', rootVolume_, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_rootsurfacePerVolumeV3.csv', rootsurfacePerVolume, delimiter=",")

		np.savetxt("./"+path+"/"+H5_filename+'_rootTipNumberV3.csv', rootTipNumber_, delimiter=",")

		np.savetxt("./"+path+"/"+H5_filename+'_uptakeRateWaterV3.csv', uptakeRateWater_, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_uptakeRatePV3.csv', uptakeRateP_, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_uptakeCumWaterV3.csv', c_uptakeRateWater, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_uptakeCumPV3.csv', c_uptakeRateP, delimiter=",")

		np.savetxt("./"+path+"/"+H5_filename+'_PCumPerNumberRootTipsV3.csv', PCumPerNumberRootTips, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_WaterCumPerNumberRootTipsV3.csv', WaterCumPerNumberRootTips, delimiter=",")

		np.savetxt("./"+path+"/"+H5_filename+'_POverWaterCumV3.csv', POverWaterCum, delimiter=",")

		np.savetxt("./"+path+"/"+H5_filename+'_PRatePerVolumeV3.csv', PRatePerVolume, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_PRatePerSurfaceV3.csv', PRatePerSurface, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_PRatePerLengthV3.csv', PRatePerLength, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_PCumPerVolumeV3.csv', PCumPerVolume,  delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_PCumPerSurfaceV3.csv', PCumPerSurface, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_PCumPerLengthV3.csv', PCumPerLength, delimiter=",")

		np.savetxt("./"+path+"/"+H5_filename+'_WaterRatePerVolumeV3.csv', WaterRatePerVolume, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_WaterRatePerSurfaceV3.csv', WaterRatePerSurface, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_WaterRatePerLengthV3.csv', WaterRatePerLength, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_WaterCumPerVolumeV3.csv', WaterCumPerVolume, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_WaterCumPerSurfaceV3.csv', WaterCumPerSurface, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_WaterCumPerLengthV3.csv', WaterCumPerLength, delimiter=",")

		np.savetxt("./"+path+"/"+H5_filename+'_waterMassFlowPerVolumeV3.csv', waterMassFlowPerVolume, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_waterMassFlowPerSurfaceV3.csv', waterMassFlowPerSurface, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_waterMassFlowPerLengthV3.csv', waterMassFlowPerLength, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_WaterTransportCumPerVolumeV3.csv', WaterTransportCumPerVolume, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_WaterTransportCumPerSurfaceV3.csv', WaterTransportCumPerSurface, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_WaterTransportCumPerLengthV3.csv', WaterTransportCumPerLength, delimiter=",")

		np.savetxt("./"+path+"/"+H5_filename+'_ConcentrationFlowPerVolumeV3.csv', ConcentrationFlowPerVolume, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_ConcentrationFlowPerSurfaceV3.csv',ConcentrationFlowPerSurface, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_ConcentrationFlowPerLengthV3.csv', ConcentrationFlowPerLength, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_NutrientTransportCumPerVolumeV3.csv', NutrientTransportCumPerVolume, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_NutrientTransportCumPerSurfaceV3.csv',NutrientTransportCumPerSurface, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_NutrientTransportCumPerLengthV3.csv', NutrientTransportCumPerLength, delimiter=",")

		np.savetxt("./"+path+"/"+H5_filename+'_ratioWaterUptakePerTransport.csv', ratioWaterUptakePerTransport, delimiter=",")
		np.savetxt("./"+path+"/"+H5_filename+'_ratioNutrientUptakePerTransport.csv', ratioNutrientUptakePerTransport, delimiter=",")

		#Plotting(H5_filename)

def Plotting(H5_filename):
		process = psutil.Process(os.getpid())
		print(process.memory_info().rss)
		print ('ploting analysis results..')
		path = H5_filename[0:H5_filename.index(".")]
		plotRootCollarWaterPotential("./"+path+"/"+H5_filename+'_rootCollarPressureV3.csv')
		plotRootGrowth("./"+path+"/"+H5_filename+'_rootSurfaceV3.csv',  "rootSurface",'m2',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_rootLengthV3.csv', "rootLength",'m',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_rootVolumeV3.csv', "rootVolume",'m3',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_rootTipNumberV3.csv', "rootTipsNumber",'-',1)

		print(process.memory_info().rss)
		plotRootGrowth("./"+path+"/"+H5_filename+'_rootsurfacePerVolumeV3.csv', "rootSurface Per rootVolume",'m2/m3',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_rootVolumeV3.csv', "rootMass",'kg', 0.15*1e-3*1e6)
		plotRootGrowth("./"+path+"/"+H5_filename+'_rootTipVolumeV3.csv', "rootTipMass",'kg', 0.15*1e-3*1e6)

		print(process.memory_info().rss)
		plotRootGrowth("./"+path+"/"+H5_filename+'_uptakeRateWaterV3.csv', "Water uptake rate",'kg/s',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_uptakeRatePV3.csv', "P uptake rate",'kg/s',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_uptakeCumWaterV3.csv', "Cumulative water uptake",'kg',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_uptakeCumPV3.csv', "Cumulative P uptake",'kg',1)

		plotRootGrowth("./"+path+"/"+H5_filename+'_uptakeRateWaterTipsV3.csv', "Root tips - Water uptake rate",'kg/s',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_uptakeRatePTipsV3.csv', "Root tips - P uptake rate",'kg/s',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_uptakeCumWaterTipsV3.csv', "Root tips - Cumulative water uptake",'kg',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_uptakeCumPTipsV3.csv', "Root tips - Cumulative P uptake",'kg',1)

		plotRootGrowth("./"+path+"/"+H5_filename+'_POverWaterCumV3.csv', "Ratio cumulative uptake of P over water",'kg/kg',1)

		plotRootGrowth("./"+path+"/"+H5_filename+'_PRatePerVolumeV3.csv', "P uptake rate per root volume",'kg/s/m3',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_PRatePerSurfaceV3.csv', "P uptake rate per root surface",'kg/s/m2',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_PRatePerLengthV3.csv', "P uptake rate per root length",'kg/s/m',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_PRatePerVolumeV3.csv', "P uptake rate per root mass",'kg/s/kg',1/(0.15*1e-3*1e6))

		print(process.memory_info().rss)
		plotRootGrowth("./"+path+"/"+H5_filename+'_PCumPerVolumeV3.csv', " Cumulative P uptake per root volume",'kg/m3',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_PCumPerSurfaceV3.csv', "Cumulative P uptake per root surface",'kg/m2',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_PCumPerLengthV3.csv', "Cumulative P uptake per root length",'kg/m',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_PCumPerVolumeV3.csv', "Cumulative P uptake per root mass",'kg/kg',1/(0.15*1e-3*1e6))

		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterRatePerVolumeV3.csv', "water uptake rate per root volume",'kg/s/m3',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterRatePerSurfaceV3.csv', "water uptake rate per root surface",'kg/s/m2',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterRatePerLengthV3.csv', "water uptake rate per root length",'kg/s/m',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterRatePerVolumeV3.csv', "water uptake rate per root mass",'kg/s/kg',1/(0.15*1e-3*1e6))

		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterCumPerVolumeV3.csv', " Cumulative water uptake per root volume",'kg/m3',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterCumPerSurfaceV3.csv', "Cumulative water uptake per root surface",'kg/m2',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterCumPerLengthV3.csv', "Cumulative water uptake per root length",'kg/m',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterCumPerVolumeV3.csv', "Cumulative water uptake per root mass",'kg/kg',1/(0.15*1e-3*1e6))
###
		plotRootGrowth("./"+path+"/"+H5_filename+'_waterMassFlowPerVolumeV3.csv',"water flow rate per root volume",'kg/s/m3',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_waterMassFlowPerSurfaceV3.csv', "water flow rate per root surface",'kg/s/m2',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_waterMassFlowPerLengthV3.csv',"water flow rate per root length",'kg/s/m',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_waterMassFlowPerVolumeV3.csv', "water flow rate per root mass",'kg/s/kg',1/(0.15*1e-3*1e6))

		print(process.memory_info().rss)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterTransportCumPerVolumeV3.csv',  " Total water transport per root volume",'kg/m3',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterTransportCumPerSurfaceV3.csv', "Total water transport per root surface",'kg/m2',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterTransportCumPerLengthV3.csv', "Total water transport per root length",'kg/m',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterTransportCumPerVolumeV3.csv', "Total water transport  per root mass",'kg/kg',1/(0.15*1e-3*1e6))
###
		plotRootGrowth("./"+path+"/"+H5_filename+'_ConcentrationFlowPerVolumeV3.csv',"Nutrient flow rate per root volume",'m3/s/m3',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_ConcentrationFlowPerSurfaceV3.csv', "Nutrient flow rate per root surface",'m3/s/m2',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_ConcentrationFlowPerLengthV3.csv',"Nutrient flow rate per root length",'m3/s/m',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_ConcentrationFlowPerVolumeV3.csv', "Nutrient flow rate per root mass",'m3/s/kg',1/(0.15*1e-3*1e6))

		print(process.memory_info().rss)
		plotRootGrowth("./"+path+"/"+H5_filename+'_NutrientTransportCumPerVolumeV3.csv', " Total nutrient transport per root volume",'kg/m3',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_NutrientTransportCumPerSurfaceV3.csv', "Total nutrient transport per root surface",'kg/m2',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_NutrientTransportCumPerLengthV3.csv', "Total nutrient transport per root length",'kg/m',1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_NutrientTransportCumPerVolumeV3.csv', "Total nutrient transport  per root mass",'kg/kg',1/(0.15*1e-3*1e6))

		plotRootGrowth("./"+path+"/"+H5_filename+'_PCumPerNumberRootTipsV3.csv', "PCumPerNumberRootTips", "kg/tip", 1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_WaterCumPerNumberRootTipsV3.csv', "WaterCumPerNumberRootTips", "kg/tip", 1)

		plotRootGrowth("./"+path+"/"+H5_filename+'_ratioWaterUptakePerTransport.csv', "ratioWaterUptakePerTransport", "-", 1)
		plotRootGrowth("./"+path+"/"+H5_filename+'_ratioNutrientUptakePerTransport.csv', "ratioNutrientUptakePerTransport", "-", 1)

def filterHF5(H5_filename):
	Analysing(H5_filename)
	Plotting(H5_filename)

if __name__=='__main__':
    sys.exit(filterHF5(sys.argv[1]))
