from utils.balanceUtils import ImportReactionFile, ExportUniqueCompoundsWithIOStatus, \
GenerateIOStatusList
import os.path
from utils.specutils12 import ensure_dir



dirName = 'input/jetfuel_reactions/'
exportDirName = 'output/jetfuel_iostatus/'

fileNames = [\
'pinene_reactions.txt', \
'limonene_reactions.txt',
'farnesene_reactions.txt',
'Bisabolene_reactions'
]

for fileName in fileNames:
	
	print(fileName)
	
	filePath = os.path.join(dirName, fileName)
	
	uniqueCompounds, reactions, sMatrixTKeyIndexed = \
	ImportReactionFile(filePath, reactionArrow='→')

	uniqueCompoundsIOStatus = GenerateIOStatusList(uniqueCompounds)

	baseName = os.path.basename(fileName)

	exportFileName = os.path.splitext(baseName)[0].split('_')[0] + '_iostatus.csv'
	exportFilePath = os.path.join(exportDirName, exportFileName)
	
	ensure_dir(exportFilePath)

	ExportUniqueCompoundsWithIOStatus(exportFilePath, uniqueCompoundsIOStatus)
