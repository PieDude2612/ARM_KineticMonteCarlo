import numpy as np
import math
from matplotlib import pyplot as plt
from Simulation_Code.propensity_fcn import propensity_fcn
from Simulation_Code.reacRatesCalc import reacRatesCalc
from Simulation_Code.dataSetLoader import dataSetLoader

# TODO: Call reaction calc on every MD separately
#  Master table of union reaction of [numerator, denominator]
#  Keep track of right indices of reactions in separate MD
#  Master dictionary plus master table of reaction links (Union to MD file)
#  Run reacRatesCalc for each file and fill entire column
#  Add up numerators and denominators, 0 if no repeat over any other file

class simulator():

    def startup(totalFiles):
        t = np.array([0, 20000])
        thresholdReact = 5
        count = 0
        listofKeys = np.array([])
        masterArr = np.array([])
        masterArrIndices = np.array([])
        listofKeys = np.append(listofKeys, 'Master')

        for filenum in range(totalFiles + 1):
            listofKeys = np.append(listofKeys, 'MD' + str(filenum + 1) + 'TH')
            listofKeys = np.append(listofKeys, 'MD' + str(filenum + 1) + 'TP')
        listofKeys = np.append(listofKeys, 'Final Reaction Rate Constants')

        master = dict.fromkeys(listofKeys, 0)

        for filenum in range(totalFiles + 1):
            reacdictN = str(np.loadtxt(open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reacdict_all'
                                    + str(filenum + 1) + '.dat', 'r')))
            for reac in range(len(reacdictN)):
                if((reacdictN[reac] != any(masterArr))):
                    masterArr = np.append(masterArr, reacdictN[reac])
                    masterArrIndices = np.append(masterArrIndices, count)
                    count = count + 1

            master['Master'] = masterArrIndices
            dataSetLoader.loadFiles(filenum + 1)
            timesHapp, timesPoss = reacRatesCalc.calcrr(dataSetLoader.xi, dataSetLoader.rfpc,
                                                    dataSetLoader.stoich_matrix, dataSetLoader.stoich_pos,
                                                    t[0], t[1], thresholdReact, dataSetLoader.mols_neg_id,
                                                    dataSetLoader.expConc)
            timesHappened = np.array([len(timesHapp)])
            timesPossible = np.array([len(timesPoss)])

            for ind in range (masterArr):
                timesHappened[ind] = timesHapp[reacdictN.index(masterArr[ind])]
                timesPossible[ind] = timesPoss[reacdictN.index(masterArr[ind])]


            master['MD' + str(filenum + 1) +'TH'] = timesHappened
            master['MD' + str(filenum + 1) + 'TP'] = timesPossible
        #TODO: Add up the TH and TP to fill in final constants key in dict. This will conclude the master table

        simulator.directMethod(dataSetLoader.stoich_matrix, t, dataSetLoader.xi[0, :], master['Final Reaction Rate Constants'], 3000)

    def directMethod(stoich_matrix, tspan, x0, reaction_rates, max_output_length):
        num_species = stoich_matrix.shape[1]
        T = np.zeros((max_output_length, 1))  # time step array
        X = np.zeros((max_output_length, num_species))  # molecules that exist over time
        MU = np.zeros((max_output_length, 1))
        T[0] = tspan[0]
        X[0, :] = x0

        rxnCount = 0

        ##################################################################################################
        while (T[rxnCount] < tspan[1]):  # as long as the time step stays within allocated time
            a = np.double([])

            a = propensity_fcn.calc_propensity(stoich_matrix, X[rxnCount, :], reaction_rates)
            asum = np.sum(a)  # take the entire sum of the prop function

            # a = np.append(a, reaction_rates[0] * X[rxnCount, 0] * X[rxnCount, 1])
            # a = np.append(a, reaction_rates[1] * X[rxnCount, 2])
            # a = np.append(a, reaction_rates[2] * X[rxnCount, 2])
            # asum = np.sum(a)

            r1 = np.random.uniform(0, 1)
            r2 = np.random.uniform(0, 1)

            tau = (1 / asum) * math.log((1 / r1), math.e)
            # tau = (math.log(1 / r1)) / asum
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
                t = T[1:rxnCount]
                graphs = 1
                colorTally = 0
                species = 0
                limit = 0

                fig = plt.figure(1)
                # Since the graph legend only has 10 colors in store, they are repeated making it hard to read
                # graph. So I designed something that will create separate graphs everytime 10 colours are used
                # This way it is easier to keep track of stuff

                while(colorTally < num_species):
                    if(num_species - colorTally < 10):
                        limit = limit + (num_species - colorTally)
                    else:
                        limit = limit + 10

                    while(species < limit):
                        ax = fig.add_subplot(1, 3, graphs)
                        ax.plot(t, X[0:(rxnCount - 1), species], label='x' + str(species))
                        species = species + 1

                    colorTally = colorTally + 10
                    graphs = graphs + 1
                    plt.xlabel("Time (s)")
                    plt.ylabel("Molecules")
                    plt.legend(loc='upper right')
                    plt.show()

                raise Exception("Simulation terminated because max output length has been reached.")
                break

        t = T[1:rxnCount]
        graphs = 1
        colorTally = 0
        species = 0
        limit = 0

        fig = plt.figure(1)
        # Since the graph legend only has 10 colors in store, they are repeated making it hard to read
        # graph. So I designed something that will create separate graphs everytime 10 colours are used
        # This way it is easier to keep track of stuff

        while (colorTally < num_species):
            if (num_species - colorTally < 10):
                limit = limit + (num_species - colorTally)
            else:
                limit = limit + 10

            while (species < limit):
                ax = fig.add_subplot(1, 3, graphs)
                ax.plot(t, X[0:(rxnCount - 1), species], label='x' + str(species))
                species = species + 1

            colorTally = colorTally + 10
            graphs = graphs + 1
            plt.xlabel("Time (s)")
            plt.ylabel("Molecules")
            plt.legend(loc='upper right')
            plt.show()
