import importlib
from Simulation_Code import simulator
from Simulation_Code import propensity_fcn
from Simulation_Code import reacRatesCalc
from Simulation_Code import dataSetLoader
from Simulation_Code import matchNcreate
importlib.reload(simulator)
importlib.reload(propensity_fcn)
importlib.reload(reacRatesCalc)
importlib.reload(dataSetLoader)
importlib.reload(matchNcreate)

from Simulation_Code.simulator import simulator

totFiles = int(input("Enter total files to train from: "))
simulationTime = int(input("Enter total simulation time in timesteps: "))
simFileNumber = int(input("Enter the file number that KMC is being used on: "))

simulator.startup(simulator(), totFiles, simulationTime, simFileNumber)