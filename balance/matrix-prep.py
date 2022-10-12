from utils.balanceUtils import ImportReactionFile, ExportUniqueCompoundsWithIOStatus, \
GenerateIOStatusList
import os.path
from utils.specutils12 import ensure_dir



dirName = 'input/jetfuel_reactions/'
exportDirName = 'output/jetfuel_iostatus/'

fileNames = [\
#'pinene_reactions.txt', \
#'limonene_reactions.txt', \
#'farnesene_reactions.txt', \
#'bisabolene_reactions.txt', \
#'geraniol_reactions.txt', \
#'fattyacid_reactions.txt', \
'hexanoicacid_reactions.txt', \
'heptanoicacid_reactions.txt', \
'octanoicacid_reactions.txt', \
'decanoicacid_reactions.txt', \
'dodecanoicacid_reactions.txt', \
'tetradecanoicacid_reactions.txt', \
'hexadecanoicacid_reactions.txt', \
'LLeucine_reactions.txt', \
'LValine_reactions.txt', \
'LIsoleucine_reactions.txt'
#'1longchain_reactions.txt', \
#'2longchain_reactions.txt', \
#'3longchain_reactions.txt', \
#'4longchain_reactions.txt'
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
