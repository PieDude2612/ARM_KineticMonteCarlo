from Simulation_Code.simulator import simulator

totFiles = int(input("Enter total files to train from: "))
simulationTime = int(input("Enter total simulation time in timesteps: "))
simFileNumber = int(input("Enter the file number that KMC is being used on: "))

simulator.startup(simulator(), totFiles, simulationTime, simFileNumber)