import numpy as np
import math
from matplotlib import pyplot as plt
from Simulation_Code.propensity_fcn import propensity_fcn
from Simulation_Code.reacRatesCalc import reacRatesCalc
from Simulation_Code.dataSetLoader import dataSetLoader
from Simulation_Code.matchNcreate import matchNcreate

class simulator():

    def startup(self, totalFiles, simTime, simFileNum):
        t = np.array([0, simTime])
        thresholdReact = 5

        masterReactionArr = np.array([]).astype(str)
        masterMoleculeArr = np.array([]).astype(str)
        reactionRateConstants = np.array([])

        dictofdicts = dict()
        keys = np.array([])
        for i in range(totalFiles + 1): keys = np.append(keys, i + 1)
        dictofdicts.fromkeys(keys)
#######################################################################################################################
        for filenum in range(totalFiles):
            reacdictN = np.array([])
            theReacsFile = open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reacdict_all'
                                    + str(filenum + 1) + '.dat', 'r')
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

        for filenum in range(totalFiles):
            timesHapp, timesPoss = reacRatesCalc.calcrr(reacRatesCalc(), dataSetLoader.xi(dataSetLoader(), filenum + 1),
                                                        dataSetLoader.rpfc(dataSetLoader(), filenum + 1),
                                                        dataSetLoader.sms(dataSetLoader(), filenum + 1, 1),
                                                        dataSetLoader.sms(dataSetLoader(), filenum + 1, 2),
                                                        t[0], t[1], thresholdReact,
                                                        dataSetLoader.sms(dataSetLoader(), filenum + 1, 3),
                                                        dataSetLoader.sms(dataSetLoader(), filenum + 1, 4))
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
#######################################################################################################################
            theMolesFile = open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\moleculedict_all'
                                         + str(filenum + 1) + '.dat', 'r')
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
#######################################################################################################################
        for reaction in range(len(masterReactionArr)):
            reactionRateConstants = np.append(reactionRateConstants, ((finalTimesHapp[reaction] / finalTimesPoss[reaction]) / 0.012))

        masterStoichMat = matchNcreate.doStringMatch(matchNcreate(), list(masterReactionArr), list(masterMoleculeArr))

        xtoCompare = dataSetLoader.xi(dataSetLoader(), 3)
        xtoTake, plotInds = dataSetLoader.createTestMD(dataSetLoader(), simFileNum, masterMoleculeArr, xtoCompare[0, :])
        keyreacSpeciesnum = dataSetLoader.sms(dataSetLoader(), simFileNum, 1).shape[1]

        simulator.iterateNplot(simulator(), masterStoichMat, t, xtoTake, xtoCompare, keyreacSpeciesnum, plotInds,
                               reactionRateConstants, 3000)

    def iterateNplot(self, stoich_matrix, tspan, x0, xcomp, keySpecs, pltInds, reaction_rates, max_output_length):
        num_species = stoich_matrix.shape[1]
        T = np.zeros((max_output_length, 1))  # time step array
        X = np.zeros((max_output_length, num_species))  # molecules that exist over time
        MU = np.zeros((max_output_length, 1))
        T[0] = tspan[0]
        X[0, :] = x0

        rxnCount = 0

        ##################################################################################################
        while (T[rxnCount] < tspan[1]):  # as long as the time step stays within allocated time
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

            if ((rxnCount + 1) >= max_output_length):  # If time allocated is still not exceeded and loop
                t = T[0:rxnCount - 1]
                graphs = 1
                colorTally = 0
                species = pltInds
                specs = 0
                limit = 0

                fig = plt.figure(1)
                comparefig = plt.figure(2)
                # Since the graph legend only has 10 colors in store, they are repeated making it hard to read
                # graph. So I designed something that will create separate graphs everytime 10 colours are used
                # This way it is easier to keep track of stuff

                while (colorTally < keySpecs):
                    if (keySpecs - colorTally < 10):
                        limit = limit + (keySpecs - colorTally)
                    else:
                        limit = limit + 10

                    ax = fig.add_subplot(10, 10, graphs)
                    ax2 = comparefig.add_subplot(10, 10, graphs)
                    while (specs < limit):
                        ax.plot(t, X[0:(rxnCount - 1), species[specs]], label='x' + str(species[specs]))
                        ax2.plot(t, xcomp[0:rxnCount - 1, specs], label='xo' + str(specs))
                        specs = specs + 1

                    colorTally = colorTally + 10
                    graphs = graphs + 1
                    fig.xlabel("Time (ps)")
                    comparefig.xlabel("Time (ps)")
                    fig.ylabel("Molecules")
                    comparefig.ylabel("Molecules")
                    fig.legend(loc='upper right')
                    comparefig.legend(loc='upper right')
                    fig.show()
                    comparefig.show()

                raise Exception("Simulation terminated because max output length has been reached.")

        t = T[0:rxnCount - 1]
        graphs = 1
        colorTally = 0
        species = pltInds
        specs = 0
        limit = 0

        fig = plt.figure(1)
        comparefig = plt.figure(2)
        # Since the graph legend only has 10 colors in store, they are repeated making it hard to read
        # graph. So I designed something that will create separate graphs everytime 10 colours are used
        # This way it is easier to keep track of stuff

        while (colorTally < keySpecs):
            if (keySpecs - colorTally < 10):
                limit = limit + (keySpecs - colorTally)
            else:
                limit = limit + 10

            ax = fig.add_subplot(10, 10, graphs)
            ax2 = comparefig.add_subplot(10, 10, graphs)
            while (specs < limit):
                ax.plot(t, X[0:(rxnCount - 1), species[specs]], label='x' + str(species[specs]))
                ax2.plot(t, xcomp[0:rxnCount - 1, specs], label='xo' + str(specs))
                specs = specs + 1

            colorTally = colorTally + 10
            graphs = graphs + 1
            fig.xlabel("Time (ps)")
            comparefig.xlabel("Time (ps)")
            fig.ylabel("Molecules")
            comparefig.ylabel("Molecules")
            fig.legend(loc='upper right')
            comparefig.legend(loc='upper right')
            fig.show()
            comparefig.show()
