import numpy as np
import math
from matplotlib import pyplot as plt
from Simulation_Code import simulator
from Simulation_Code import propensity_fcn
from Simulation_Code import reacRatesCalc


class simulator():
    def startup(num):
        # import all data from external files
        all_data = np.loadtxt(open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactbasis_all.dat', 'r'),
                              usecols=range(3))
        rows = (all_data[:, 0]).astype(int)
        cols = (all_data[:, 1]).astype(int)
        mols = (all_data[:, 2]).astype(int)

        stoich_matrix = np.zeros((np.amax(rows), np.amax(cols)))

        index = 0  # load in all data into the specific positions in complete matrices

        while (index < (len(rows) - 1)):  # load data into actual sm
            sparr = rows[index] - 1
            sparc = cols[index] - 1
            sparv = mols[index]
            stoich_matrix[sparr, sparc] = sparv
            index = index + 1

        all_data.clear()
        rows.clear()
        cols.clear()

        all_data = np.loadtxt(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\molhistperframe_1.dat', 'r'),
            usecols=range(3))
        rows = (all_data[:, 0]).astype(int)
        cols = (all_data[:, 1]).astype(int)
        mols = (all_data[:, 2]).astype(int)

        xi = np.zeros((np.amax(rows), np.amax(cols)))

        index = 0

        while (index < (len(rows) - 1)):  # load data into actual xi
            sparr = rows[index] - 1
            sparc = cols[index] - 1
            sparv = mols[index]
            xi[sparr, sparc] = sparv
            index = index + 1

        all_data.clear()
        rows.clear()
        cols.clear()

        all_data = np.loadtxt(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactperframe_1.dat', 'r'),
            usecols=range(3))
        rows = (all_data[:, 0]).astype(int)
        cols = (all_data[:, 1]).astype(int)
        mols = (all_data[:, 2]).astype(int)

        rfpc = np.zeros((np.amax(rows), np.amax(cols)))

        index = 0

        while (index < (len(rows) - 1)):  # load data into actual sm
            sparr = rows[index] - 1
            sparc = cols[index] - 1
            sparv = mols[index]
            rfpc[sparr, sparc] = sparv
            index = index + 1

        all_data.clear()
        rows.clear()
        cols.clear()

        # recycle the variables to reduce total memory allocated in the loading process by clearing all data from vars
        # every time the matrices have been loaded

        nA = 6.023e23
        volu = 1e-15
        t = np.array([0, 50])

        reaction_rates = reacRatesCalc.calcrr(xi, rfpc, stoich_matrix, t[0], t[1])

        simulator.directMethod(stoich_matrix, t, xi, reaction_rates, volu, 1000)

    def directMethod(stoich_matrix, tspan, x0, reaction_rates, volume, max_output_length):
        num_species = stoich_matrix.shape[1]
        T = np.zeros((max_output_length, 1))  # time step array
        X = np.zeros((max_output_length, num_species))  # molecules that exist over time
        MU = np.zeros((max_output_length, 1))
        T[0] = tspan[0]
        X[0, 0:len(x0)] = x0

        rxnCount = 0

        ##################################################################################################
        while (T[rxnCount] < tspan[1]):  # as long as the time step stays within allocated time
            a = np.double([])

            a = propensity_fcn.calc_propensity(stoich_matrix, X[rxnCount, :], reaction_rates, volume)
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

            if ((rxnCount + 1) >= max_output_length):  # If time allocated is still not exceeded and loop
                t = T[1:rxnCount]
                x = X[1:rxnCount, 0]
                x2 = X[1:rxnCount, 1]
                x3 = X[1:rxnCount, 2]
                x4 = X[1:rxnCount, 3]

                plt.figure(1)
                plt.xlabel("Time (s)")
                plt.ylabel("Molecules")
                plt.plot(t, x, t, x2, t, x3, t, x4)
                plt.legend(["x0", "x1", "x3", "x4"])
                plt.show()

                raise Exception("Simulation terminated because max output length has been reached.")
                break
        ###################################################################################################

        t = T[1:rxnCount]
        x = X[1:rxnCount, 0]
        x2 = X[1:rxnCount, 1]
        x3 = X[1:rxnCount, 2]
        x4 = X[1:rxnCount, 3]

        plt.figure(1)
        plt.xlabel("Time (s)")
        plt.ylabel("Molecules")
        plt.plot(t, x, t, x2, t, x3, t, x4)
        plt.legend(["x0", "x1", "x3", "x4"])
        plt.show()
