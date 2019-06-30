#!/home/h.mai/share/tools/anaconda2/bin/python2
import sys
import numpy
import matplotlib
matplotlib.use('Agg')
import math
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from scipy import interpolate
from scipy import integrate
from scipy import interp
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import os.path

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

def plotRhizosVariable(fileName, soilLayer, variableName, unit, color):
	values = numpy.genfromtxt(fileName, delimiter=',')
	timeSteps = values[:,0]
	#https://matplotlib.org/examples/pylab_examples/multicolored_line.html
	bound05 = values[:,1]
	bound25 = values[:,2]
	bound50 = values[:,3]
	bound75 = values[:,4]
	bound95 = values[:,5]

	#meanPoints = numpy.array([timeSteps/24/3600, bound50]).T.reshape(-1, 1, 2)
	#lowerPoints25 = numpy.array([timeSteps/24/3600, bound25]).T.reshape(-1, 1, 2)
	#upperPoints75 = numpy.array([timeSteps/24/3600, bound75]).T.reshape(-1, 1, 2)
	#lowerPoints05 = numpy.array([timeSteps/24/3600, bound05]).T.reshape(-1, 1, 2)
	#upperPoints95 = numpy.array([timeSteps/24/3600, bound95]).T.reshape(-1, 1, 2)
#
	#meanSegments = numpy.concatenate([meanPoints[:-1], meanPoints[1:]], axis=1)
	#lowerSegments25 = numpy.concatenate([lowerPoints25[:-1], lowerPoints25[1:]], axis=1)
	#upperSegments75 = numpy.concatenate([upperPoints75[:-1], upperPoints75[1:]], axis=1)
	#lowerSegments05 = numpy.concatenate([lowerPoints05[:-1], lowerPoints05[1:]], axis=1)
	#upperSegments95 = numpy.concatenate([upperPoints95[:-1], upperPoints95[1:]], axis=1)
#
	## SETTING CMAP
	##https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
	#Z = [[0,0],[0,0]]
	##soilLayer = [-0.00,-0.2]
	#levels = numpy.arange(soilLayer[1], soilLayer[0], (soilLayer[0]-soilLayer[1])/12)
	#cmap = plt.get_cmap('jet')
	#CS3 = plt.contourf(Z, levels, cmap=cmap)
#
	#norm = plt.Normalize(0, 1)
	#lc_mean = LineCollection(meanSegments, cmap=cmap, norm=norm)
	##lc_lower = LineCollection(lowerSegments, cmap=cmap, norm=norm)
	##lc_upper = LineCollection(upperSegments, cmap=cmap, norm=norm)
	## Set the values used for colormapping
	#meanDepth = values[:,-1]
	#lc_mean.set_array(numpy.asarray(meanDepth))
	#lc_mean.set_linewidth(2)
#
#
	## Ploting
	#sns.set()
	##ax = plt.subplot(111, xlabel='days', ylabel='concentrations at root surface', title='uptake rate - '+simulationName)
	#ax = plt.subplot(111, xlabel='days', ylabel=unit, title= variableName+'\n'+fileName)
	#for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
	#             ax.get_xticklabels() + ax.get_yticklabels()):
	#    item.set_fontsize(12)
	#plt.gca().add_collection(lc_mean)
#
	#for i in range(len(timeSteps)-1):
	#	plt.fill_between(timeSteps[i:i+2]/24/3600, [lowerSegments05[i][0][1],lowerSegments05[i][1][1]], [upperSegments95[i][0][1],upperSegments95[i][1][1]], color=cmap(meanDepth[i]), alpha=.3, facecolor="none", linewidth=0.0)
	#	plt.fill_between(timeSteps[i:i+2]/24/3600, [lowerSegments25[i][0][1],lowerSegments25[i][1][1]], [upperSegments75[i][0][1],upperSegments75[i][1][1]], color=cmap(meanDepth[i]), alpha=.5, facecolor="none", linewidth=0.0)
	#plt.axis([min(timeSteps/24/3600),max(timeSteps/24/3600),min(bound05)*0.75,max(bound95)*1.1])
	#clb = plt.colorbar(CS3)
	#clb.ax.set_title('depth (m)')
	#plt.show()

	# Plotting
	sns.set()
	title = variableName+'\n'+fileName[fileName.rfind('/')+1:len(fileName)]
	ax = plt.subplot(111, xlabel='days', ylabel=unit, title= title)
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
	             ax.get_xticklabels() + ax.get_yticklabels()):
	    item.set_fontsize(12)

	endPoint = len(timeSteps)-3
	plt.plot(timeSteps[0:endPoint]/24/3600, bound50[0:endPoint], ls="-", color= color, linewidth=2)
	plt.fill_between(timeSteps[0:endPoint]/24/3600, bound05[0:endPoint], bound95[0:endPoint], color= color, alpha = 0.3)
	plt.fill_between(timeSteps[0:endPoint]/24/3600, bound25[0:endPoint], bound75[0:endPoint], color= color, alpha = 0.5)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.savefig(fileName+'.png', bbox_inches='tight', dpi=300)
	##plt.colorbar(cmap)
	#clb = plt.colorbar(CS3)
	#clb.ax.set_title('depth (m)')
	#clb.set_label('depth (m)', rotation=270)
	#cbar.set_label('depth (m)', rotation=270)
	#plt.show()

def statisticalAnalysisRhizo(simulationName, variableIdx, soilLayer, rootTypes, outputFile):
	#Macroscale Data
	uptakeRateMacroScale = numpy.loadtxt(simulationName+"-soil_uptakeRate.log")
	timeRefMacroScale = numpy.unique(uptakeRateMacroScale[:,0])#/24/3600
	####

	#reading Rhizosphere information
	#rootEIdx	soilEIdx	Position	Length	Radius	Vector
	ReadRhizo = numpy.loadtxt(simulationName+"-root_pointsource.log",skiprows=1)
	rootEIdxCol = 0
	soilEIdxCol = 1
	PositionCol = [2,3,4]
	LengthCol = 5
	rootRadiusCol = 6
	VectorCol = [7,8,9]
	rootOrderCol = 10
	rootBranchCol = 11
	rootBornTimeCol = 12
	distanceFromOriginCol = 13

	# Filtering rhizo files
	RhizoInfo = ReadRhizo
	print (RhizoInfo.shape)
	rowIdx = 0
	while (rowIdx < RhizoInfo.shape[0]):
		if not (((RhizoInfo[rowIdx][4] < soilLayer[0]) and (RhizoInfo[rowIdx][4] > soilLayer[1])) and (RhizoInfo[rowIdx][rootOrderCol] in rootTypes)):#and(RhizoInfo[rowIdx][0]>=1):
			#print RhizoInfo[rowIdx][0]
			RhizoInfo = numpy.delete(RhizoInfo, (rowIdx), axis=0)
			rowIdx -=1
		rowIdx +=1
	print (RhizoInfo.shape)

	#RhizoInfo = ReadRhizo[ReadRhizo[:,PositionCol[2]].argsort()][::-1] #sorting array based on Z from high to low value

	# Collecting rhizo files
	RhizosphereFileName = []
	RhizosphereUptake =[]
	for i in range(0, RhizoInfo.shape[0]):
		filename_ = simulationName+"_rhizo/r"+str(int(RhizoInfo[i][rootEIdxCol]))+"_s"+str(int(RhizoInfo[i][soilEIdxCol]))+"/r_uptakeRate.log"
		if (os.path.isfile(filename_)):
			RhizosphereFileName.append(filename_)
			RhizosphereUptake.append(numpy.loadtxt(filename_))
	#values = [[]]*len(timeRefMacroScale)
	#https://stackoverflow.com/questions/13447882/python-2-7-creating-a-multidimensional-list
	for k in range(0,len(variableIdx)):
		#initialize lists for value
		values = [[[] for i in range(0)] for i in range(len(timeRefMacroScale))]
		depth = [[[] for i in range(0)] for i in range(len(timeRefMacroScale))]
		#c = numpy.arange(0, 1)
		#for i in range(0, RhizoInfo.shape[0]):

		# collecting values
		for i in range(1, len(RhizosphereUptake)):
			relativeDepth = 1-(RhizoInfo[i][4]-soilLayer[0])/(soilLayer[1]-soilLayer[0])
			#data = RhizosphereUptake[i]
			if (RhizosphereUptake[i].size > 15):
				for j in range(0,len(RhizosphereUptake[i][:,0])):
					timeIndex = numpy.where(timeRefMacroScale==RhizosphereUptake[i][j,0])[0]
					if (len(timeIndex)!= 0):
						#print data[j,2], timeIndex[0], values[timeIndex[0]]
						values[timeIndex[0]].extend([RhizosphereUptake[i][j,variableIdx[k]]])
						depth[timeIndex[0]].extend([relativeDepth])

		#statistical analysis
		meanValue = []
		bound05 = []
		bound25 = []
		bound75 = []
		bound95 = []
		meanDepth = []
		for i in range(0,len(values)):
			#print len(values[i])
			if (len(values[i])==0):
				meanValue.append(0)
				meanDepth.append(0)
				bound05.append(0)
				bound25.append(0)
				bound75.append(0)
				bound95.append(0)
			else:
			#if (len(values[i])!=0):
				#meanValue.append(numpy.percentile(values[i],50))
				meanValue.append(numpy.mean(values[i]))
				meanDepth.append(numpy.mean(depth[i]))
				bound05.append(numpy.percentile(values[i],5))
				bound25.append(numpy.percentile(values[i],25))
				bound75.append(numpy.percentile(values[i],75))
				bound95.append(numpy.percentile(values[i],95))

		outPut = numpy.array([timeRefMacroScale, bound05, bound25, meanValue, bound75, bound95, meanDepth]).T
		numpy.savetxt(outputFile[k], outPut, delimiter=",")

def plotRhizos(fileName):
	simulationName = fileName[1]
	if (not os.path.exists('extract')):
		os.mkdir('extract');

	if (not os.path.exists('extract/'+simulationName+'_rhizo')):
		os.mkdir('extract/'+simulationName+'_rhizo');

	#variableIdxs = [2,3,4,7,6,8,9,10,11]
	#variableNames = ['rhizoNutrientUptakeRates', 'rhizoConcentrations', 'rhizoDepletionRadii','rhizoDomain','relativeRhizoDepletionRadii','EffDiffCoeff',"waterPressureInterface","waterContentInterface","waterPressureMacro"]
	#units = ['kg/s','kg/m3','m','m','m']
	#colors = ['red', 'blue', 'brown','green','orange']

	#variableIdxs = [7,6]
	#variableNames = ['rhizoDomain','relativeRhizoDepletionRadii']
	#units = ['m','-']
	#colors = ['green','orange']

	variableIdxs = [2,3,4,7,6,8,10,12,13]
	variableNames = ['rhizoNutrientUptakeRates', 'rhizoConcentrations', 'rhizoDepletionRadii','rhizoDomain','relativeRhizoDepletionRadii','EffDiffCoeff','WaterContentAtRootSurface','WaterContentMacro','PelectNo']
	units = ['kg/s/m2','kg/m3','m','m','-','m2/s','-','-','-']
	colors = ['red', 'blue', 'brown','green','orange','black','pink','gray','yellow']

	topSoil = [0,-0.2]
	subSoil = [-0.2,-0.6]
	totalSoil = [0,-0.6]
	#soilLayers = [topSoil, subSoil,totalSoil]
	#soilLayerNames = ['TopSoil','SubSoil','TotalSoil']
	soilLayers = [totalSoil]
	soilLayerNames = ['TotalSoil']

	totalRoot = [0,1,2,3,4,5]
	nodalRoot = [0,1,4]
	STypeRoot = [5]
	LTypeRoot = [2,3]
	#rootTypes = [nodalRoot, STypeRoot, LTypeRoot, totalRoot]
	#rootTypeNames = ['NodalRoot', 'STypeRoot', 'LTypeRoot', 'TotalRoot']
	rootTypes = [totalRoot]
	rootTypeNames = ['TotalRoot']

	for i in range(0,len(soilLayers)):
		for j in range(0,len(rootTypes)):
			print ('statisticalAnalysisRhizos.. '+soilLayerNames[i]+' '+rootTypeNames[j])
			outputFileNames = []
			for variableName_ in variableNames:
				outputFileNames.append('extract/'+simulationName+'_rhizo/'+simulationName+'_'+variableName_+soilLayerNames[i]+rootTypeNames[j]+'.csv')
			#print outputFileNames
			statisticalAnalysisRhizo(simulationName, variableIdxs, soilLayers[i], rootTypes[j], outputFileNames)

	for i in range(0,len(soilLayers)):
		for j in range(0,len(rootTypes)):
			print ('plotting RhizosVariables.. '+soilLayerNames[i]+' '+rootTypeNames[j])
			outputFileNames = []
			for variableName_ in variableNames:
				outputFileNames.append('extract/'+simulationName+'_rhizo/'+simulationName+'_'+variableName_+soilLayerNames[i]+rootTypeNames[j]+'.csv')
			#print outputFileNames
			for k in range(0,len(variableNames)):
				if (os.path.exists(outputFileNames[k])):
					plotRhizosVariable(outputFileNames[k],soilLayers[i], variableNames[k]+soilLayerNames[i]+rootTypeNames[j], units[k], colors[k])
				else:
					print('	'+outputFileNames[k]+' does not exist !!')
if __name__=='__main__':
    sys.exit(plotRhizos(sys.argv))