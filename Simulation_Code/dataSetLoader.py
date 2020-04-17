import numpy as np

class dataSetLoader():

    def loadFiles(fileNumber):
        theData = open(
            'D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactbasis_all' + str(fileNumber) + '.dat', 'r')
        all_data = np.array(theData.readlines())
        # TODO: Finish correct formatting and loading of all_data
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
                        # have to append ID as many times as exponent for math
                        expConc[r, ind] = exponent
                        # need to go from 0 to whatever exponent for math
                        exponent = exponent + 1
                        ind = ind + 1

        all_data = None
        rows = None
        cols = None
        mols = None
        mols_pos = None
        rows_neg = None
        cols_neg = None

        all_data = np.load(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\molhistperframe_' + str(fileNumber) + '.dat',
                 'r'),
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

        all_data = np.load(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactperframe_' + str(fileNumber) + '.dat', 'r'),
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

    def loadMDFile(theFile):
        all_data = np.load(
            open(
                'D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\molhistperframe_' + str(theFile) + '.dat',
                'r'),
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

        return xi

    def loadStoichMat(theFile):
        all_data = np.load(
            open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactbasis_all' + str(theFile) + '.dat',
                 'r'),
            usecols=range(3))
        rows = (all_data[:, 0]).astype(int)
        cols = (all_data[:, 1]).astype(int)
        mols = (all_data[:, 2]).astype(int)

        stoich_matrix = np.zeros((np.amax(rows), np.amax(cols))).astype(int)

        for index in range(len(rows)):  # load data into actual sm
            sparr = rows[index] - 1
            sparc = cols[index] - 1
            sparv = mols[index]

            stoich_matrix[sparr, sparc] = sparv

        all_data = None
        rows = None
        cols = None
        mols = None

        return stoich_matrix