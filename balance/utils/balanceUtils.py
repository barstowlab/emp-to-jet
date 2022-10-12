from numpy import matmul
import scipy.optimize
from scipy.optimize import fmin_slsqp
from utils.specutils12 import GenerateFileList


# ------------------------------------------------------------------------------------------------ #
def ImportStoichiometricMatrix(fileName):
	
	import numpy
	import gc
	from pdb import set_trace
	import csv
		
	fHandle = open(fileName, 'r')
		
	i = 0
	datareader = csv.reader(fHandle)
	matrix = []

	for row in datareader:
		matrix.append(row)

	fHandle.close()
	
	return matrix

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CheckConcentrationChangeVector(concentrationChangeVector, inputOutputStatus):
# Check the change of concentration vector

	i = 0
	indicesWhereIntermediateNotStable = []
	intermediateStableVector = []

	while i < len(concentrationChangeVector):
		if inputOutputStatus[i] == 'Intermediate':
			if concentrationChangeVector[i] == 0:
				intermediateStableVectorEntry = True
			else:
				intermediateStableVectorEntry = False
				indicesWhereIntermediateNotStable.append(i)
		else:
			intermediateStableVectorEntry = 'NA'
	
		intermediateStableVector.append(intermediateStableVectorEntry)
		i += 1
	
	return intermediateStableVector, indicesWhereIntermediateNotStable
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateStoichiometricMatrix(inputMatrix, startIndex, endIndex):
	
	from numpy import transpose, matrix, array
	import pdb
	
	matrixT = transpose(inputMatrix)

	
	reactions = matrixT[1][startIndex:]
	reactants = inputMatrix[0][startIndex:endIndex]
	ioStatus = inputMatrix[1][startIndex:endIndex]
	sMatrixT = []

	i = startIndex
	while i < len(inputMatrix):
		j = startIndex
		row = []
		while j < endIndex:
			element = inputMatrix[i][j]
			if element == '':
				elementFloat = 0
			else:
				try:
					elementFloat = float(element)
				except:
					pdb.set_trace()
			
			row.append(elementFloat)
			j += 1
	
		sMatrixT.append(row)
	
		i += 1

	sMatrix = transpose(sMatrixT)
	sMatrix = array(sMatrix)

	return sMatrix, reactions, reactants, ioStatus
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def w_rss1(fluxVector, sMatrix, ioStatus):
    
    from numpy import array
    import pdb
    import numpy
    
    predictedCDotVector = numpy.dot(sMatrix, fluxVector)
     
    errors = []
    
    i = 0
    while i < len(predictedCDotVector):
    	if ioStatus[i] == 'Intermediate':
    		error = 0 - predictedCDotVector[i]
    		#print('Intermediate')
    		errors.append(error)
    	#print(i)
    	i += 1
    
    errorArray = array(errors)
    
    rss = (errorArray**2).sum()
    
    #pdb.set_trace()
    
    return rss
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ieqconstraint(fluxVector, sMatrix, ioStatus):
    
    from numpy import array
    import pdb
    import numpy
    
    predictedCDotVector = numpy.dot(sMatrix, fluxVector)
     
    nTargetCompound = 0
        
    i = 0
    while i < len(predictedCDotVector):
    	if ioStatus[i] == 'Target':
    		nTargetCompound = predictedCDotVector[i]
    		#print(nTargetCompound)
    	i += 1
    
    
    return nTargetCompound
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def PrintStoichiometry(cDotVector, reactants, ioStatus, printIntermediates=True):
	
	i = 0
	while i < len(cDotVector):
		outputString = ''
		
		
		if printIntermediates == True:
			outputString += reactants[i] + '\t' + str(cDotVector[i]) + '\t' + ioStatus[i]
			print(outputString)
		else:
			if ioStatus[i] != 'Intermediate':
				outputString += reactants[i] + '\t' + str(cDotVector[i])
				print(outputString)
		
		i += 1
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
class FindTargetIndexFailure(Exception):
	pass
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def FindTargetIndex(ioStatus):
	
	import pdb
	
	targetIndices = []
	
	i = 0
	while i < len(ioStatus):
		
		if ioStatus[i] == 'Target':
			targetIndices.append(i)
				
		i += 1
	
	if len(targetIndices) == 0:
		print('No target compound')
		ex = FindTargetIndexFailure('Number of prPool entries wrong')
		raise ex 
	elif len(targetIndices) == 1:
		returnIndex = targetIndices[0]
	else:
		ex = FindTargetIndexFailure('Number of prPool entries wrong')
		raise ex 
		
	
	return returnIndex

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def SolveFluxBalanceEquation(sMatrix, reactions, reactants, ioStatus):
	
	from numpy import ones, dot
	import pdb
	
	try:
		targetIndex = FindTargetIndex(ioStatus)
	except:
		pdb.set_trace()
		
	# Calculate the change of concentration vector
	fVector0 = ones(len(reactions), float)
	cDotVector0 = matmul(sMatrix, fVector0)

	# Basically, I need to find the flux vector that results in the minimum value of concentration
	# change vector.

	rss = w_rss1(fVector0, sMatrix, ioStatus)

	result = fmin_slsqp(w_rss1, fVector0, args=(sMatrix, ioStatus), iprint=0, \
	ieqcons=[ieqconstraint])

	fVectorOpt = result / 1.0

	cDotVectorOpt = dot(sMatrix, fVectorOpt)
	
	normalizationFactor = cDotVectorOpt[targetIndex]
	
	cDotVectorOptNorm = cDotVectorOpt / normalizationFactor


	return [fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, result]
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def BalanceStoichiometricMatrix(matrix, startIndex, endIndex):

	sMatrix, reactions, reactants, ioStatus = \
	GenerateStoichiometricMatrix(matrix, startIndex, endIndex)
	
	[fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, result] = \
	SolveFluxBalanceEquation(sMatrix, reactions, reactants, ioStatus)

	return [sMatrix, reactions, fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, reactants, ioStatus, \
	result]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportBalanceAndReportStoichiometricMatrix(directory, fileName, startIndex, endIndexOffset, \
outputDict, cDotVectorDict, printIntermediates=False):

	import pdb

	dictKey = fileName.split('-')[1]
	dictKey = dictKey.split('.')[0]
	
	print(dictKey)


	matrix = ImportStoichiometricMatrix(directory + '/' + fileName)
	startIndex = 2
	endIndex = len(matrix[0]) - endIndexOffset

	[sMatrix, reactions, fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, reactants, ioStatus, \
	result] = \
	BalanceStoichiometricMatrix(matrix, startIndex, endIndex)
	
	outputDict[dictKey] = [sMatrix, reactions, fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, \
	reactants, ioStatus]
	cDotVectorDict[dictKey] = cDotVectorOpt
	
	# PrintStoichiometry(cDotVectorOptNorm, reactants, ioStatus, printIntermediates=True)
# 	print()
# 	
	PrintStoichiometry(cDotVectorOptNorm, reactants, ioStatus, printIntermediates=printIntermediates)
	print()
	
	atpIndex = reactants.index('ATP')
	nadhIndex = reactants.index('NADH')
	try:
		fdredIindex = reactants.index('Fdred')
	except:
		pdb.set_trace()
	
	nATP = cDotVectorOptNorm[atpIndex]
	nFdred = cDotVectorOptNorm[fdredIindex]
	nNADH = cDotVectorOptNorm[nadhIndex]
	
	
	return [dictKey, nATP, nFdred, nNADH]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportBalanceAndReportStoichiometricMatrix2(directory, fileName, startIndex, endIndexOffset, \
outputDict, cDotVectorDict, printIntermediates=False, reactantsToGet=None):

	import pdb

	dictKey = fileName.split('-')[1]
	dictKey = dictKey.split('.')[0]
	
	print(dictKey)


	matrix = ImportStoichiometricMatrix(directory + '/' + fileName)
	startIndex = 2
	endIndex = len(matrix[0]) - endIndexOffset

	[sMatrix, reactions, fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, reactants, ioStatus, \
	result] = \
	BalanceStoichiometricMatrix(matrix, startIndex, endIndex)
	
	outputDict[dictKey] = [sMatrix, reactions, fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, \
	reactants, ioStatus]
	cDotVectorDict[dictKey] = cDotVectorOpt
	
	PrintStoichiometry(cDotVectorOptNorm, reactants, ioStatus, \
	printIntermediates=printIntermediates)
	print()
	
	nReactantsDict = {}
	nTargetsDict = {}
	targetsToGet = []
	
	
	if reactantsToGet == None:
		reactantsToGet = []
		i = 0
		while i < len(reactants):
			ioStatusReactant = ioStatus[i]
			if ioStatusReactant != 'Intermediate':
				reactantsToGet.append(reactants[i])
			
			if ioStatusReactant == 'Target':
				targetsToGet.append(reactants[i])
			
			i += 1
		
	for reactant in reactantsToGet:	
		try:
			reactantIndex = reactants.index(reactant)
			nReactant = cDotVectorOptNorm[reactantIndex]
			
		except:
			nReactant = 0
			
		nReactantsDict[reactant] = nReactant
	
	for target in targetsToGet:
		targetIndex = reactants.index(target)
		nTarget = cDotVectorOptNorm[targetIndex]
		nTargetsDict[target] = nTarget 
	
	
	
	return dictKey, nReactantsDict, nTargetsDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportBalanceAndOutputStoichiometricMatrices(directory, regex=".*CBB\.csv", ignoreCase=False, \
printIntermediates=False):
	
	from specutils11 import GenerateFileList

	
	fileList = GenerateFileList(directory=directory, regex=regex, ignoreCase=ignoreCase)

	cDotVectorDict = {}
	outputDict = {}
	dictKeyArray = []
	nATPArray = []
	nFdredArray = []
	nNADHArray = []

	for fileName in fileList:
	
		[dictKey, nATP, nFdred, nNADH] = \
		ImportBalanceAndReportStoichiometricMatrix(directory, fileName, 2, 3, outputDict, \
		cDotVectorDict, printIntermediates=printIntermediates)
	
		dictKeyArray.append(dictKey)
		nATPArray.append(nATP)
		nFdredArray.append(nFdred)
		nNADHArray.append(nNADH)
	
	
	return [dictKeyArray, nATPArray, nNADHArray, nFdredArray]

# ------------------------------------------------------------------------------------------------ #


# ----------------------------------------------------------------------------------------------- #
def FindUniqueTerms(sideSplit):

	compoundArray = []
	for term in sideSplit:
		termSplit = term.split('*')
		if len(termSplit) == 1:
			compound = termSplit[0].strip()
			compoundArray.append(compound)
		elif len(termSplit) == 2:
			compound = termSplit[1].strip()
			compoundArray.append(compound)
		else:
			print("Something weird is going on, length of termSplit is longer than 2.")
	
	return compoundArray
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def GenerateUniqueCompoundsList(dataSplit):
	
	from numpy import unique
	
	compoundArray = []

	for line in dataSplit:
		lhs = line[0].split('+')
		rhs = line[1].split('+')
	
		lhsCompounds = FindUniqueTerms(lhs)
		rhsCompounds = FindUniqueTerms(rhs)
	
		compoundArray += lhsCompounds
		compoundArray += rhsCompounds

	uniqueCompounds = unique(compoundArray)

	return uniqueCompounds
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def UpdateSMatrixWithMultipliers(sideSplit, signMultiplier, sMatrixT, rowIndex):

	import pdb

	for term in sideSplit:
		termSplit = term.split('*')
		if len(termSplit) == 1:
			compound = termSplit[0].strip()
			multiplier = 1*signMultiplier
		elif len(termSplit) == 2:
			compound = termSplit[1].strip()
			multiplier = int(termSplit[0].strip())*signMultiplier
		else:
			print("Something weird is going on, length of termSplit is longer than 2.")
	
		try:
			sMatrixT[rowIndex][compound] += multiplier
		except:
			pdb.set_trace()
	
	return
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def ConvertIndexedSMatrix(sMatrixTKeyIndexed, uniqueCompounds):
	
	from numpy import int8, zeros
	import pdb
	
	nCols = len(uniqueCompounds)
	nRows = sMatrixTKeyIndexed.shape[0]
	sMatrixT = zeros((nRows, nCols), dtype=int8)
	
	
	i = 0
	while i < nRows:
		j = 0
		while j < nCols:
			sMatrixT[i][j] = sMatrixTKeyIndexed[i][uniqueCompounds[j]]
			j += 1
		i += 1
	

	return sMatrixT
# ----------------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------------- #
def GenerateIndexedSMatrixT(compoundsUnique, reactions, reactionArrow='→'):
	
	from numpy import int8, zeros


	# Start splitting the chemical equations into left and right hand sides
	dataSplit = []

	for line in reactions:
		lineSplit = line.split(reactionArrow)
		dataSplit.append(lineSplit)
	
	# Make a list of all of the compounds involved
	uniqueCompounds = GenerateUniqueCompoundsList(dataSplit)


	# Generate a dtypeArray for the stoichiometric matrix
	dtypeArray = []
	i = 0
	while i < len(uniqueCompounds):
		dtypeArray.append((uniqueCompounds[i], int8))
		i += 1

	sMatrixTKeyIndexed = zeros(len(reactions), dtype=dtypeArray)


	# Populate the stoichiometric matrix
	i = 0
	while i < len(dataSplit):
		line = dataSplit[i]
		lhs = line[0].split('+')
		rhs = line[1].split('+')
	
		UpdateSMatrixWithMultipliers(lhs, -1, sMatrixTKeyIndexed, i)
		UpdateSMatrixWithMultipliers(rhs, 1, sMatrixTKeyIndexed, i)
	
		i += 1

	
	return sMatrixTKeyIndexed
# ----------------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------------- #
def ImportReactionFile(filename, reactionArrow='→'):
	
	from numpy import int8, zeros
	
	# Open up the file
	fileHandle = open(filename, 'r')
	data = fileHandle.readlines()


	# Clean up the file a bit
	dataClean = []

	for line in data:
		dataClean.append(line.strip())


	# Start splitting the chemical equations into left and right hand sides
	dataSplit = []

	for line in dataClean:
		lineSplit = line.split(reactionArrow)
		dataSplit.append(lineSplit)
	
	# Make a list of all of the compounds involved
	uniqueCompounds = GenerateUniqueCompoundsList(dataSplit)


	# Generate a dtypeArray for the stoichiometric matrix
	dtypeArray = []
	i = 0
	while i < len(uniqueCompounds):
		dtypeArray.append((uniqueCompounds[i], int8))
		i += 1

	sMatrixTKeyIndexed = zeros(len(dataClean), dtype=dtypeArray)


	# Populate the stoichiometric matrix
	i = 0
	while i < len(dataSplit):
		line = dataSplit[i]
		lhs = line[0].split('+')
		rhs = line[1].split('+')
	
		UpdateSMatrixWithMultipliers(lhs, -1, sMatrixTKeyIndexed, i)
		UpdateSMatrixWithMultipliers(rhs, 1, sMatrixTKeyIndexed, i)
	
		i += 1
	
	reactions = dataClean
		
	return uniqueCompounds, reactions, sMatrixTKeyIndexed
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def ExportUniqueCompoundsWithIOStatus(fileName, uniqueCompoundsIOStatus):
	
	fileHandle = open(fileName, 'w')
	
	for line in uniqueCompoundsIOStatus:
		outputStr = '"' + line[0] + '"' + ',' + line [1] + '\n' 
		fileHandle.write(outputStr)
	
	fileHandle.close()
	
	return
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def GenerateIOStatusList(uniqueCompounds):
	
	uniqueCompoundsIOStatus = []
	
	for compound in uniqueCompounds:
		
		if compound == 'ATP':
			ioStatus = 'Input'
		elif compound == 'NADH':
			ioStatus = 'Input'
		elif compound == 'Fdred':
			ioStatus = 'Input'
		elif compound == 'CO2':
			ioStatus = 'Input/Output'
		elif compound == 'HCO3-':
			ioStatus = 'Input/Output'
		elif compound == 'HCOO-':
			ioStatus = 'Input/Output'
		elif compound == 'H2O':
			ioStatus = 'Input/Output'
		elif compound == 'N2':
			ioStatus == 'Input'
		elif compound == 'HCO2-':
			ioStatus == 'Input/Output'
		elif compound == 'H2':
			ioStatus == 'Input/Output'
		else:
			ioStatus = 'Intermediate'
		
		uniqueCompoundsIOStatus.append([compound, ioStatus])
	
	return uniqueCompoundsIOStatus
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def ImportIOStatus(fileName):
	
	import csv
	
	fHandle = open(fileName, 'r')
	
	poolColumnToHeaderIndexDict = {}	
	datareader = csv.reader(fHandle)

	uniqueCompoundsIOStatus = []
	
	for row in datareader:
		uniqueCompoundsIOStatus.append(row)

	return uniqueCompoundsIOStatus
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def GenerateMergedIOStatusList(uniqueCompounds, iostatusCO2, iostatusAA):

	ioStatusCO2Dict = {}
	ioStatusAADict = {}

	for line in iostatusCO2:
		ioStatusCO2Dict[line[0]] = line[1]
		
	for line in iostatusAA:
		ioStatusAADict[line[0]] = line[1]
		
	
	ioStatusCO2DictKeys = ioStatusCO2Dict.keys()
	ioStatusAADictKeys = ioStatusAADict.keys()
	
	mergedIOStatus = []
	
	for compound in uniqueCompounds:
		
		if compound in ioStatusCO2DictKeys and compound not in ioStatusAADictKeys:
			mergedIOStatus.append(ioStatusCO2Dict[compound])
		
		elif compound not in ioStatusCO2DictKeys and compound in ioStatusAADictKeys:
			mergedIOStatus.append(ioStatusAADict[compound])
		
		elif compound in ioStatusCO2DictKeys and compound in ioStatusAADictKeys:
			mergedIOStatus.append(ioStatusCO2Dict[compound])
			
		else:
			mergedIOStatus.append('Intermediate')
	
		
	return mergedIOStatus, ioStatusCO2Dict, ioStatusAADict
	
	
	
# ----------------------------------------------------------------------------------------------- #

