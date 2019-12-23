import numpy as np
import math

class reacRatesCalc():
    def calcrr(x0, reactperFrame, stoich_mat, stoich_mat_pos, startFrame, endFrame):
        global ireactmat
        global concexp

        numreacts = stoich_mat.shape[1]
        reactRates = np.zeros(numreacts)
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
                    dkdt = np.append(dkdt, reactperFrame[i, 1])  # array of 1 and 0
                    continue  # how many times it occurs found by adding up array

                if (np.sum(dkdt) == 0):  # if reaction given does not occur
                    break
                # TODO: Does it make more sense to keep everything below this outside second for loop?

                ireact = np.nonzero(stoich_mat_pos[reac, :])  # take nonzero value indices out
                reactionOrder = stoich_mat_pos[i, ireact]  # input their reac num here.

                reacOrderMatrix = np.ones(1, len(tsteps - backt)) * reactionOrder  # matrix operation

                # TODO: Understand the matrix operations in MATLAB code before proceeding. Might need loops.
