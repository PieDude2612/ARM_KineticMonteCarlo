import numpy as np
import math
import time
from matplotlib import pyplot as plt
from Simulation_Code.propensity_fcn import propensity_fcn
from Simulation_Code.reacRatesCalc import reacRatesCalc
from Simulation_Code.dataSetLoader import dataSetLoader
from Simulation_Code.matchNcreate import matchNcreate

class simulator():

    def startup(self, filesArray, simTime, simFileNum):
        print(str(time.ctime(time.time())) + ": Starting Program Execution")
        t = np.array([0, simTime])
        thresholdReact = 5

        masterReactionArr = np.array([]).astype(str)
        masterMoleculeArr = np.array([]).astype(str)
        reactionRateConstants = np.array([])

        dictofdicts = dict()
        keys = np.array([])
        for i in range(len(filesArray)): keys = np.append(keys, i + 1)
        dictofdicts.fromkeys(keys)
#######################################################################################################################
        for filenum in range(len(filesArray)):
            reacdictN = np.array([])
            theReacsFile = open('D:\\ARM_KineticMonteCarlo\\Data Files\\reacdict_all'
                                    + str(filesArray[filenum]) + '.dat', 'r')
            for line in theReacsFile:
                reacdictN = np.append(reacdictN, line.replace('\n', ''))
            dictofdicts[filenum] = list(reacdictN)
            # find the file, which should be named a certain way

            for reac in range(len(reacdictN)):
                if(filenum == 0):
                    masterReactionArr = np.append(masterReactionArr, reacdictN[reac])
                else:
                    try:
                        x = list(masterReactionArr).index(reacdictN[reac])
                    except ValueError or AttributeError:
                        masterReactionArr = np.append(masterReactionArr, reacdictN[reac])
            # loading processes for reaction and molecule dictionaries
        finalTimesHapp = np.zeros((len(masterReactionArr))).astype(int)
        finalTimesPoss = np.zeros((len(masterReactionArr))).astype(int)
        print(str(time.ctime(time.time())) + ": All dictionaries loaded")

        for filenum in range(len(filesArray)):
            timesHapp, timesPoss = reacRatesCalc.calcrr(reacRatesCalc(), dataSetLoader.xi(dataSetLoader(), filesArray[filenum]),
                                                        dataSetLoader.rpfc(dataSetLoader(), filesArray[filenum]),
                                                        dataSetLoader.sms(dataSetLoader(), filesArray[filenum], 1),
                                                        dataSetLoader.sms(dataSetLoader(), filesArray[filenum], 2),
                                                        t[0], t[1], thresholdReact,
                                                        dataSetLoader.sms(dataSetLoader(), filesArray[filenum], 3),
                                                        dataSetLoader.sms(dataSetLoader(), filesArray[filenum], 4))
            # get the reaction rate constants for each MD file

            timesHappened = np.array([]).astype(int)
            timesPossible = np.array([]).astype(int)
            for ind in range(len(masterReactionArr)): # rearrange reactions for final alignment
                try:
                    timesHappened = np.append(timesHappened, timesHapp[dictofdicts[filenum].index(masterReactionArr[ind])])
                    timesPossible = np.append(timesPossible, timesPoss[dictofdicts[filenum].index(masterReactionArr[ind])])
                except ValueError:
                    timesHappened = np.append(timesHappened, 0)
                    timesPossible = np.append(timesPossible, 0)
                except AttributeError:
                    timesHappened = np.append(timesHappened, 0)
                    timesPossible = np.append(timesPossible, 0)

            finalTimesHapp = finalTimesHapp + timesHappened
            finalTimesPoss = finalTimesPoss + timesPossible
            print(str(time.ctime(time.time())) + ": Calculated and appended timesHapp and timesPoss " + str(filesArray[filenum]))
#######################################################################################################################
            theMolesFile = open('D:\\ARM_KineticMonteCarlo\\Data Files\\moleculedict_all'
                                         + str(filesArray[filenum]) + '.dat', 'r')
            moledictN = np.array([])
            for line in theMolesFile:
                moledictN = np.append(moledictN, line.replace('\n', ''))

            for reac in range(len(moledictN)):
                if (filenum == 0):
                    masterMoleculeArr = np.append(masterMoleculeArr, moledictN[reac])
                else:
                    try:
                        x = list(masterMoleculeArr).index(moledictN[reac])
                    except ValueError or AttributeError:
                        masterMoleculeArr = np.append(masterMoleculeArr, moledictN[reac])
            print(str(time.ctime(time.time())) + ": Appended molecule file " + str(filesArray[filenum]))
#######################################################################################################################
        for reaction in range(len(masterReactionArr)):
            reactionRateConstants = np.append(reactionRateConstants, ((finalTimesHapp[reaction] / finalTimesPoss[reaction]) / 0.012))
        print(str(time.ctime(time.time())) + ": Finished calculating reaction rate constants")

        masterStoichMat = matchNcreate.doStringMatch(matchNcreate(), list(masterReactionArr), list(masterMoleculeArr))
        print(str(time.ctime(time.time())) + ": Stoich matrix created")

        xtoCompare = dataSetLoader.xi(dataSetLoader(), simFileNum)
        xtoTake, plotInds, speciesNotTrained = dataSetLoader.createTestMD(dataSetLoader(), simFileNum, masterMoleculeArr, xtoCompare[0, :])
        xtoCompare = np.delete(xtoCompare, speciesNotTrained, axis=1)
        print(str(time.ctime(time.time())) + ": Calculated test MD")

        simulator.iterateNplot(simulator(), masterStoichMat, t, xtoTake, xtoCompare, masterMoleculeArr, plotInds,
                               reactionRateConstants)

    def iterateNplot(self, stoich_matrix, tspan, x0, xcomp, masterMolArr, pltInds, reaction_rates):
        print(str(time.ctime(time.time()))+ ": Started simulation")
        num_species = stoich_matrix.shape[1]
        T = np.zeros((tspan[1], 1))  # time step array
        X = np.zeros((tspan[1], num_species))  # molecules that exist over time
        MU = np.zeros((tspan[1], 1))
        T[0] = tspan[0]
        X[0, :] = x0

        rxnCount = 0
        print(str(time.ctime(time.time()))+ ": Simulating propensity function...")
        ##################################################################################################
        while (rxnCount < tspan[1] - 1):  # as long as the time step stays within allocated time
            a = propensity_fcn.calc_propensity(propensity_fcn(), stoich_matrix, X[rxnCount, :], reaction_rates)
            asum = np.sum(a)  # take the entire sum of the prop function

            r1 = np.random.uniform(0, 1)
            tau = (1 / asum) * math.log((1 / r1), math.e)
            mu = 0
            ai = a[0]  # for first loop

            while (ai < (r1 * asum)):
                mu = mu + 1
                ai = ai + a[mu]

                if ((mu + 1) == len(a)):
                    break

            T[rxnCount + 1] = T[rxnCount] + tau  # the next time step determined by PRV
            X[rxnCount + 1, :] = X[rxnCount, :] + stoich_matrix[mu, :]  # the next molecule count (PRV)
            MU[rxnCount + 1] = mu
            rxnCount = rxnCount + 1
        ###################################################################################################
        print(str(time.ctime(time.time())) + ": Finished simulation.")
        print("Type in any number in given range to see graphs. STOP to exit.")
        while 1 == 1:
            spectoSee = str(input("Enter species to analyse (0-" + str(len(pltInds)) + ") : "))
            if (spectoSee == "STOP"):
                break
            print(str(time.ctime(time.time())) + ": Generating picture...")

            t = T[0:rxnCount]
            yvalsKMC = X[0:rxnCount, pltInds[int(spectoSee)]]
            yvalsMD = xcomp[0:rxnCount, int(spectoSee)]

            plt.figure()
            plt.plot(t, yvalsKMC, label=str(masterMolArr[pltInds[int(spectoSee)]]) + " KMC") # KMC plot
            plt.plot(t, yvalsMD, label=str(masterMolArr[pltInds[int(spectoSee)]]) + " MD") # MD Plot
            plt.xlabel("Time (ps)")
            plt.ylabel("Molecules")
            plt.legend(loc='upper right')
            plt.title("Graph showing KMC and MD trajectories of " + masterMolArr[pltInds[int(spectoSee)]])
            plt.show()