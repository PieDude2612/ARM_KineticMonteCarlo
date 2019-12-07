import numpy as np
import math

class reacRatesCalc():
    def calcrr(x0, reactperFrame, stoich_mat, stoich_mat_pos, startFrame, endFrame):
        global ireactmat
        global concexp

        numreacts = stoich_mat.shape[1]
        reactRates = np.zeros(numreacts, 1)
        tsteps = (endFrame - startFrame) + 1
        backt = 1  # what is this?

        sumr = np.zeros(len(reactRates))
        sumc = np.zeros(len(reactRates))
        stdkr = np.zeros(len(reactRates, 2))
        lenyesr = np.zeros(len(reactRates))  # this?

        dkdt = np.array([])

        for reac in range(numreacts):
            for i in range(reactperFrame.shape[0]):
                if (reactperFrame[i, 0] == reac):
                    dkdt = np.append(dkdt, reactperFrame[i, 1])
                    continue

                if (np.sum(dkdt) == 0):
                    break

                ireact = np.nonzero(stoich_mat_pos[i, :])  # this?
                concexpirereact = stoich_mat_pos[i, ireact]  # what is this line doing?
                # conexpirereactmat = np.ones(length(1:tsteps - backt), 1) *conexpirereact
