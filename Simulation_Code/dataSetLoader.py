import numpy as np

class dataSetLoader():

    def sms(self, fileNumber, type):
        theData = open(
            'D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactbasis_all' + str(fileNumber) + '.dat', 'r')
        all_data = theData.readlines()
        rows = np.array([]).astype(int)
        cols = np.array([]).astype(int)
        mols = np.array([]).astype(int)
        for line in all_data:
            theLine = np.array(line.split(' ')).astype(int)
            rows = np.append(rows, theLine[0])
            cols = np.append(cols, theLine[1])
            mols = np.append(mols, theLine[2])

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

        mols_pos = None
        rows_neg = None
        cols_neg = None
        if(type == 1):
            return stoich_matrix
        elif(type == 2):
            return stoich_pos
        elif(type == 3):
            return mols_neg_id
        elif(type == 4):
            return expConc
        return 1

    def xi(self, fileNumber):
        theData = open(
            'D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\molhistperframe_' + str(fileNumber) + '.dat',
            'r')
        all_data = theData.readlines()
        rows = np.array([]).astype(int)
        cols = np.array([]).astype(int)
        mols = np.array([]).astype(int)
        for line in all_data:
            theLine = np.array(line.split(' ')).astype(int)
            rows = np.append(rows, theLine[0])
            cols = np.append(cols, theLine[1])
            mols = np.append(mols, theLine[2])

        xi = np.zeros((np.amax(rows), np.amax(cols))).astype(int)
        index = 0

        while (index < (len(rows) - 1)):  # load data into actual xi
            sparr = rows[index] - 1
            sparc = cols[index] - 1
            sparv = mols[index]
            xi[sparr, sparc] = sparv
            index = index + 1
        return xi

    def rpfc(self, fileNumber):
        theData = open(
            'D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\reactperframe_' + str(fileNumber) + '.dat', 'r')
        all_data = theData.readlines()
        timestep = np.array([]).astype(int)
        reacNum = np.array([]).astype(int)
        reacCount = np.array([]).astype(int)
        for line in all_data:
            theLine = np.array(line.split(' ')).astype(int)
            timestep = np.append(timestep, theLine[0])
            reacNum = np.append(reacNum, theLine[1])
            reacCount = np.append(reacCount, theLine[2])

        rpfc = np.zeros((np.amax(timestep), np.amax(reacNum))).astype(int)

        for tsp in range(len(timestep)):
            rpfc[timestep[tsp] - 1, reacNum[tsp] - 1] = reacCount[tsp]

        timestep = None
        reacNum = None
        reacCount = None
        return rpfc

    def createTestMD(self, inputFilenum, masterMolArr, xin):
        inputFile = open('D:\\PythonProgramming\\ARM_KineticMonteCarlo\\Data Files\\moleculedict_all' +
                         str(inputFilenum) + '.dat', 'r')
        testMD = np.zeros((len(masterMolArr)))
        plotIndices = np.array([]).astype(int)
        count = 0
        for line in inputFile:
            try:
                indexMatch = list(masterMolArr).index(line.replace('\n', ''))
                plotIndices = np.append(plotIndices, indexMatch)
                testMD[indexMatch] = xin[count]
                count = count + 1
            except ValueError:
                count = count + 1

        return testMD, plotIndices