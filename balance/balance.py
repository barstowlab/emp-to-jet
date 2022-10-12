from utils.vectorOutput import generateOutputMatrixWithHeaders, writeOutputMatrix
from utils.specutils12 import GenerateFileList, ensure_dir
import pdb
from utils.balanceUtils import ImportStoichiometricMatrix, BalanceStoichiometricMatrix, \
PrintStoichiometry
from numpy import array


#directory = 'input/Amino Acid Synthesis Pathways/'
directory = 'input/glucose/'

outputFileName = 'output/glucose/GLUC_CBB.csv'
ensure_dir(outputFileName)

# Find all of the matrices using the Calvin cycle
regText = r'.*CBB\.csv'
reactantsToGet = ['ATP', 'NADH', 'Fdred', 'CO2', 'HCO3-', 'HCO2-']


fileList = GenerateFileList(directory=directory, regex=regText, ignoreCase=True)
# fileList = fileList[0]

startIndex = 2
endIndexOffset = 3
printIntermediates = False

# ------------------------------------------------------------------------------------------------ #
# Balance stoichiometric matrices
dictKeyArray = []
nReactantsDict = {}
nTargetsArray = []
targetsArray = []


# Initialize storage for number of reactants 
for reactant in reactantsToGet:
	nReactantsDict[reactant] = []


# Balance each stoichiometric matrix and record reactant numbers
for fileName in fileList:

	dictKey = fileName.split('.')[0]
	
	print(dictKey)
	dictKeyArray.append(dictKey)


	matrix = ImportStoichiometricMatrix(directory + '/' + fileName)
	endIndex = len(matrix[0]) - endIndexOffset
	
	[sMatrix, reactions, fVectorOpt, cDotVectorOpt, cDotVectorOptNorm, reactants, ioStatus, \
	result] = \
	BalanceStoichiometricMatrix(matrix, startIndex, endIndex)

	PrintStoichiometry(cDotVectorOptNorm, reactants, ioStatus, \
	printIntermediates=printIntermediates)
	print()
	
	
	for reactantToGet in reactantsToGet:
		# If the reactant to query isn't in the list of reactants for this matrix, record 0
		if reactantToGet not in reactants:
			nReactantsDict[reactantToGet].append(0)
			
		# On the other hand, if it is, add the number of if to its array
		elif reactantToGet in reactants:
			reactantIndex = reactants.index(reactantToGet)
			nReactant = cDotVectorOptNorm[reactantIndex]
			nReactantsDict[reactantToGet].append(nReactant)
	
	i = 0
	while i < len(reactants):
		ioStatusReactant = ioStatus[i]
				
		if ioStatusReactant == 'Target':
			targetsArray.append(reactants[i])
			nTarget = cDotVectorOptNorm[i]
			nTargetsArray.append(nTarget)
		i += 1
	
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Output results of balancing 

headers = ['Scenario'] + reactantsToGet + ['Target Molecule', 'Target']
vectorList = [dictKeyArray]

for reactantToGet in reactantsToGet:
	vectorList.append(array(nReactantsDict[reactantToGet], int)*-1)

vectorList.append(targetsArray)
vectorList.append(array(nTargetsArray, int))

oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
writeOutputMatrix(outputFileName, oMatrix)
# ------------------------------------------------------------------------------------------------ #

