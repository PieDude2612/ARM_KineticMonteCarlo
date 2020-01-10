import numpy as np
import math

# rows are reactions. columns are species involved in the reaction
class reacRatesCalc():
    def calcrr(x0, reactperFrame, stoich_mat, stoich_mat_pos, startFrame, endFrame, rarereactLimit, negativeID,
               concenExp):
        ireactmat = negativeID  # negative indices from stoich matrix
        concexp = concenExp  # values in the negative indices from stoich matrix (in n-1 power)

        numreacts = stoich_mat.shape[0]
        reactRates = np.zeros(numreacts)
        tsteps = (endFrame - startFrame) + 1
        backt = 1

        timesHappened = np.zeros(len(reactRates))
        timesPossible = np.zeros(len(reactRates))
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

            reactionReady = np.where(reactionNotReady[0:reactionNotReady.shape[0], 0:reactionNotReady.shape[1]] == 0)
            yesreact = reactionReady + backt

            xr = x0[reactionReady + startFrame - 1, ireactmat[reac, 0:np.sum(reactionOrder)]]
            xr = xr - (np.multiply(np.ones(len(yesreact), 1), concexp[reac, 0:np.sum[reactionOrder]]))
            concreact = np.prod(xr, axis=1)

            for ind in range(len(concreact)):
                if (concreact[ind] < 0):
                    concreact[ind] = 0

            timesHappened[reac] = timesHappened[reac] + np.sum(dkdt[yesreact])
            lenyesr[reac] = lenyesr[reac] + len(yesreact)
            timesPossible[reac] = timesPossible[reac] + np.sum(concreact)
        ##############################################################################################################################
        for j in range(len(reactRates)):  # k = (integral) dkdt
            reactRates[j] = (timesHappened[j] / timesPossible[j]) / tsteps

        for val in range(len(reactRates)):  # eliminate rare events through threshold
            if (timesHappened[val] < rarereactLimit):
                reactRates[val] = 0

        # Check for invalid reaction rates
        for k in range(len(reactRates)):
            if (reactRates[k] >= math.inf):
                print("Error. Contains NaN k")
                break
            if (sum(reactRates) > (1 / tsteps)):
                print("Error. Some k values are too large.")
                break
