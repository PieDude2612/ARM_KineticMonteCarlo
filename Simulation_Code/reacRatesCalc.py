import numpy as np
import math


# rows are reactions. columns are species involved in the reaction
class reacRatesCalc():
    def calcrr(x0, reactperFrame, stoich_mat, stoich_mat_pos, startFrame, endFrame, rarereactLimit):
        global ireactmat
        global concexp

        numreacts = stoich_mat.shape[0]
        reactRates = np.zeros(numreacts)
        tsteps = (endFrame - startFrame) + 1
        backt = 1  # what is this?

        sumr = np.zeros(len(reactRates))
        sumc = np.zeros(len(reactRates))
        stdkr = np.zeros(len(reactRates, 2))  # what is this?
        lenyesr = np.zeros(len(reactRates))
        ##############################################################################################################################
        for reac in range(numreacts):
            dkdt = np.array([])

            for i in range(reactperFrame.shape[0]):
                if (reactperFrame[i, 0] == reac):
                    dkdt = np.append(dkdt, reactperFrame[i, 1])  # array of 1 and 0
                    continue  # how many times it occurs found by adding up array

                if (np.sum(dkdt) == 0):  # if reaction given does not occur
                    break

            reactionOrderIndex = np.nonzero(stoich_mat_pos[reac, :])  # take nonzero value indices out
            reactionOrder = stoich_mat_pos[
                reac, reactionOrderIndex[0:len(reactionOrderIndex)]]  # input their reac num here.

            reacOrderMatrix = np.multiply(np.ones(1, len(tsteps - backt)), reactionOrder)  # matrix operation

            reactionNotReady = x0[startFrame:endFrame, reactionOrderIndex[0:len(reactionOrderIndex)]] < \
                               reacOrderMatrix[0:reacOrderMatrix.shape[0], 0:reacOrderMatrix.shape[1]]

            reactionReady = np.where(reactionNotReady == 0)[0]
            yesreact = reactionReady + backt

            xr = np.ones(len(yesreact), np.sum(reactionOrder))
            xr = x0[reactionReady + startFrame - 1, ireactmat[reac, 0:np.sum(reactionOrder)]]
            xr = xr - (np.multiply(np.ones(len(yesreact), 1), concexp[reac, 0:np.sum[reactionOrder]]))
            concreact = np.prod(xr, axis=1)

            for ind in range(len(concreact)):
                if (concreact[ind] < 0):
                    concreact[ind] = 0

            sumr[reac] = sumr[reac] + np.sum(dkdt[yesreact])
            lenyesr[reac] = lenyesr[reac] + len(yesreact)
            sumc[reac] = sumc[reac] + np.sum(concreact)
        ##############################################################################################################################
        for j in range(len(reactRates)):  # k = (integral) dkdt
            reactRates[j] = (sumr[j] / sumc[j]) / tsteps

        for val in range(len(reactRates)):  # eliminate rare events through threshold
            if (sumr[val] < rarereactLimit):
                reactRates[val] = 0

        # Check for invalid reaction rates
        for k in range(len(reactRates)):
            if (reactRates[k] >= math.inf):
                print("Error. Contains NaN k")
                break
            if (sum(reactRates) > (1 / tsteps)):
                print("Error. Some k values are too large.")
                break
