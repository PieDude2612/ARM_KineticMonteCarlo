import numpy as np
import math
from matplotlib import pyplot as plt
from Simulation_Code.propensity_fcn import propensity_fcn
from Simulation_Code.reacRatesCalc import reacRatesCalc

# TODO: Call reaction calc on every MD separately
#  Master table of union reaction of [numerator, denominator]
#  Keep track of right indices of reactions in separate MD
#  Master dictionary plus master table of reaction links (Union to MD file)
#  Run reacRatesCalc for each file and fill entire column
#  Add up numerators and denominators, 0 if no repeat over any other file

class simulator():

    def startup(self):
        # import all data from external files
        all_data = np.loadtxt(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactbasis_all.dat', 'r'),
            usecols=range(3))
        rows = (all_data[:, 0]).astype(int)
        cols = (all_data[:, 1]).astype(int)
        mols = (all_data[:, 2]).astype(int)

        mols_pos = np.array([]).astype(int)
        rows_neg = np.array([]).astype(int)
        cols_neg = np.array([]).astype(int)
        for n in range(len(mols)):
            if (mols[n] < 0):
                mols_pos = np.append(mols_pos, np.absolute(mols[n]))
                rows_neg = np.append(rows_neg, rows[n])
                cols_neg = np.append(cols_neg, cols[n])

        stoich_matrix = np.zeros((np.amax(rows), np.amax(cols))).astype(int)
        stoich_pos = np.zeros((np.amax(rows), np.amax(cols))).astype(int)

        for index in range(len(rows)):  # load data into actual sm
            sparr = rows[index] - 1
            sparc = cols[index] - 1
            sparv = mols[index]

            stoich_matrix[sparr, sparc] = sparv

        for indexx in range(len(rows_neg)):
            sparr_pos = rows_neg[indexx] - 1
            sparc_pos = cols_neg[indexx] - 1
            sparv_pos = mols_pos[indexx]

            stoich_pos[sparr_pos, sparc_pos] = sparv_pos

        mols_neg_id = np.zeros((np.amax(rows), np.amax(cols))).astype(int)
        expConc = np.zeros((np.amax(rows), np.amax(cols))).astype(int)

        for r in range(stoich_matrix.shape[0]):
            ind = 0
            for c in range(stoich_matrix.shape[1]):
                if (stoich_matrix[r, c] < 0):
                    exponent = 0
                    for num in range(np.absolute(stoich_matrix[r, c])):
                        mols_neg_id[r, ind] = c + 1
                        #have to append ID as many times as exponent for math
                        expConc[r, ind] = exponent
                        #need to go from 0 to whatever exponent for math
                        exponent = exponent + 1
                        ind = ind + 1

        all_data = None
        rows = None
        cols = None
        mols = None
        mols_pos = None
        rows_neg = None
        cols_neg = None

        all_data = np.loadtxt(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\molhistperframe_1.dat', 'r'),
            usecols=range(3))
        rows = (all_data[:, 0]).astype(int)
        cols = (all_data[:, 1]).astype(int)
        mols = (all_data[:, 2]).astype(int)

        xi = np.zeros((np.amax(rows), np.amax(cols))).astype(int)

        index = 0

        while (index < (len(rows) - 1)):  # load data into actual xi
            sparr = rows[index] - 1
            sparc = cols[index] - 1
            sparv = mols[index]
            xi[sparr, sparc] = sparv
            index = index + 1

        all_data = None
        rows = None
        cols = None
        mols = None

        all_data = np.loadtxt(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactperframe_1.dat', 'r'),
            usecols=range(3))
        timestep = (all_data[:, 0]).astype(int)
        reacNum = (all_data[:, 1]).astype(int)
        reacCount = (all_data[:, 2]).astype(int)

        rfpc = np.zeros((np.amax(timestep), np.amax(reacNum))).astype(int)

        for tsp in range(len(timestep)):
            rfpc[timestep[tsp] - 1, reacNum[tsp] - 1] = reacCount[tsp]

        all_data = None
        timestep = None
        reacNum = None
        reacCount = None

        # recycle the variables to reduce total memory allocated in the loading process by clearing all data from vars
        # every time the matrices have been loaded

        nA = 6.023e23
        volu = 1e-15
        t = np.array([0, 20000])
        thresholdReact = 5

        reaction_rates = reacRatesCalc.calcrr(xi, rfpc, stoich_matrix, stoich_pos, t[0], t[1], thresholdReact,
                                              mols_neg_id, expConc)

        simulator.directMethod(stoich_matrix, t, xi[0, :], reaction_rates, 3000)

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
