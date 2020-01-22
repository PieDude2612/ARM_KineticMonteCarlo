import numpy as np
import math
from matplotlib import pyplot as plt
from Simulation_Code.propensity_fcn import propensity_fcn
from Simulation_Code.reacRatesCalc import reacRatesCalc


class simulator():

    def startup(self):
        # import all data from external files
        all_data = np.loadtxt(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactbasis_all.dat', 'r'),
            usecols=range(3))
        rows = (all_data[:, 0]).astype(int)
        cols = (all_data[:, 1]).astype(int)
        mols = (all_data[:, 2]).astype(int)

        mols_pos = np.array([])
        mols_neg_id = np.array([])
        mols_neg_id_col = np.array([])
        expConc = np.array([])
        expConc_col = np.array([])
        for n in range(len(mols)):
            if (mols[n] < 0):
                mols_pos = np.append(mols_pos, np.absolute(mols[n]))
                continue
            else:
                mols_pos = np.append(mols_pos, 0)
                continue

        stoich_matrix = np.zeros((np.amax(rows), np.amax(cols)))
        stoich_pos = np.zeros((np.amax(rows), np.amax(cols)))

        index = 0  # load in all data into the specific positions in complete matrices

        while (index < (len(rows) - 1)):  # load data into actual sm
            sparr = rows[index] - 1
            sparc = cols[index] - 1
            sparv = mols[index]
            sparv_pos = mols_pos[index]

            stoich_matrix[sparr, sparc] = sparv
            stoich_pos[sparr, sparc] = sparv_pos
            index = index + 1

        for r in range(stoich_matrix.shape[0]):
            cVal = 0
            for c in range(stoich_matrix.shape[1]):
                if (stoich_matrix[r, c] < 0):
                    mols_neg_id_col = np.append(mols_neg_id_col, c)  # take all index values for reactants
                    expConc_col = np.append(expConc_col, stoich_matrix[r, c])  # take their value
                    cVal = cVal + 1
                    continue
            mols_neg_id = np.append(mols_neg_id, mols_neg_id_col)  # add the column value array to the big array
            expConc = np.append(expConc, expConc_col)  # creates a matrix sort of

        # all_data.clear()
        # rows.clear()
        # cols.clear()
        # mols.clear()
        # mols_neg_id_col.clear()
        # expConc_col.clear()

        all_data2 = np.loadtxt(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\molhistperframe_1.dat', 'r'),
            usecols=range(3))
        rows2 = (all_data2[:, 0]).astype(int)
        cols2 = (all_data2[:, 1]).astype(int)
        mols2 = (all_data2[:, 2]).astype(int)

        xi = np.zeros((np.amax(rows2), np.amax(cols2)))

        index = 0

        while (index < (len(rows) - 1)):  # load data into actual xi
            sparr = rows2[index] - 1
            sparc = cols2[index] - 1
            sparv = mols2[index]
            xi[sparr, sparc] = sparv
            index = index + 1

        # all_data.clear()
        # rows.clear()
        # cols.clear()
        # mols.clear()

        all_data3 = np.loadtxt(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactperframe_1.dat', 'r'),
            usecols=range(3))
        timestep = (all_data3[:, 0]).astype(int)
        reacNum = (all_data3[:, 1]).astype(int)
        reacCount = (all_data3[:, 2]).astype(int)

        rfpc = np.zeros((np.amax(timestep), np.amax(reacNum)))
        print(rfpc.shape[0])
        print(rfpc.shape[1])

        for tsp in range(len(timestep)):
            rfpc[timestep[tsp] - 1, reacNum[tsp] - 1] = reacCount[tsp]

        # all_data.clear()
        # rowNum.clear()
        # reacCount.clear()

        # recycle the variables to reduce total memory allocated in the loading process by clearing all data from vars
        # every time the matrices have been loaded

        nA = 6.023e23
        volu = 1e-15
        t = np.array([0, 50])
        thresholdReact = 5

        reaction_rates = reacRatesCalc.calcrr(xi, rfpc, stoich_matrix, stoich_pos, t[0], t[1], thresholdReact,
                                              mols_neg_id, expConc)

        simulator.directMethod(stoich_matrix, t, xi, reaction_rates, 1000)

    def directMethod(stoich_matrix, tspan, x0, reaction_rates, max_output_length):
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
                plt.figure(1)
                plt.xlabel("Time (s)")
                plt.ylabel("Molecules")

                for species in range(num_species):
                    plt.plot(t, X[0:(rxnCount - 1), species])

                plt.show()

                raise Exception("Simulation terminated because max output length has been reached.")
                break

        t = T[1:rxnCount]
        plt.xlabel("Time (s)")
        plt.ylabel("Molecules")

        for species in range((num_species / 13)):
            plt.plot(t, X[0:(rxnCount - 1), species])

        plt.show()
