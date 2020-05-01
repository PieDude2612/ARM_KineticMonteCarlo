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
import numpy as np

totFiles = int(input("Enter total files to train from: "))
simulationTime = int(input("Enter total simulation time in timesteps: "))
simFileNumber = int(input("Enter the file number that KMC is being used on: "))

filesArr = np.arange(totFiles).astype(int)
filesArr[0:len(filesArr)] = filesArr[0:len(filesArr)] + 1
filesArr = np.delete(filesArr, np.where(filesArr==simFileNumber))

simulator.startup(simulator(), filesArr, simulationTime, simFileNumber)