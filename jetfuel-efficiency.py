#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# jetfuel-efficiency.py
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Bargraph, Generate_EfficienciesDict_Keys_Sorted_by_Efficiency, \
Export_Efficiency_Bargraph

from rewiredcarbon.utils import ensure_dir



scenarioTableFileName = 'input/jetfuel-efficiency - termination 2.csv'
outputFilenameEff = 'output/jetfuel-efficiency - termination 2.csv'
outputFilenameFuelMassEff = 'output/jetfuel-mass - termination 2.csv'
outputFilenameFuelMoleEff = 'output/jetfuel-mole - termination 2.csv'

ensure_dir(outputFilenameEff)
ensure_dir(outputFilenameFuelMassEff)
ensure_dir(outputFilenameFuelMoleEff)

scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict)



keysArray = list(efficienciesDict.keys())

Plot_Efficiency_Bargraph(efficienciesDict, 'effTotalElectricalToFuel', \
'effTotalElectricalToFuel_lowerError', 'effTotalElectricalToFuel_upperError', keysToPlot=keysArray)


Plot_Efficiency_Bargraph(efficienciesDict, 'effTotalElectricalFuelMassEfficiency', \
'effTotalElectricalFuelMassEfficiency_lowerError', \
'effTotalElectricalFuelMassEfficiency_upperError', keysToPlot=keysArray)


#Plot_Efficiency_Bargraph(efficienciesDict, 'effTotalElectricalFuelMoleEfficiency', \
#'effTotalElectricalFuelMoleEfficiency_lowerError', \
#'effTotalElectricalFuelMoleEfficiency_upperError', keysToPlot=keysArray)



Export_Efficiency_Bargraph(outputFilenameEff, efficienciesDict, scenarioDict, \
'effTotalElectricalToFuel', 'effTotalElectricalToFuel_lowerError', \
'effTotalElectricalToFuel_upperError', keysToPlot=keysArray)

Export_Efficiency_Bargraph(outputFilenameFuelMassEff, efficienciesDict, scenarioDict, \
'effTotalElectricalFuelMassEfficiency', 'effTotalElectricalFuelMassEfficiency_lowerError', \
'effTotalElectricalFuelMassEfficiency_upperError', keysToPlot=keysArray)

#Export_Efficiency_Bargraph(outputFilenameFuelMoleEff, efficienciesDict, scenarioDict, \
#'effTotalElectricalFuelMoleEfficiency', 'effTotalElectricalFuelMoleEfficiency_lowerError', \
#'effTotalElectricalFuelMoleEfficiency_upperError', keysToPlot=keysArray)
