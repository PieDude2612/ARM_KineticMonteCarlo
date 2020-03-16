import numpy as np
import math
import sys

# rows are reactions. columns are species involved in the reaction
class reacRatesCalc():
    def calcrr(x0, reactperFrame, stoich_mat, stoich_mat_pos, startFrame, endFrame, rarereactLimit, negativeID,
               concenExp):
        ireactmat = negativeID  # negative indices from stoich matrix
        concexp = concenExp  # values in the negative indices from stoich matrix (in n-1 power)

        numreacts = stoich_mat.shape[0]
        reactRates = np.zeros(numreacts)

        timesHappened = np.zeros(len(reactRates))
        timesPossible = np.zeros(len(reactRates))
        ##############################################################################################################################
        for reac in range(numreacts):
            dkdt = np.sum(reactperFrame[:, reac])
            reactionOrderIndex = np.nonzero(stoich_mat_pos[reac, :])
            reactionOrder = stoich_mat_pos[reac, reactionOrderIndex]

            reactionNotReady = np.array([]).astype(int)
            reactionReady = np.array([]).astype(int)

            for timestep in range(x0.shape[0]): #find out if a reaction is ready or not
                isIt = all(x0[timestep, reactionOrderIndex[0]] >= reactionOrder[0])
                if (isIt):
                    reactionReady = np.append(reactionReady, timestep)
                else:
                    reactionNotReady = np.append(reactionNotReady, timestep)

            fortheCalc = ireactmat[reac, np.nonzero(ireactmat[reac, :])]
            fortheCalc = fortheCalc - 1
            xr = x0[reactionReady, :][:, fortheCalc[0]].astype(int)

            reactantLengths = np.array(np.nonzero(ireactmat[reac, :]))
            # get all reactants of reactions that are ready
            theOrderMatrix = np.multiply(np.ones((len(reactionReady), 1)), concexp[reac, 0:reactantLengths.shape[1]])
            xr = np.subtract(xr, theOrderMatrix) # operands could not be broadcast together with shapes 0, 26
            concreact = np.array([])
            concreact = np.append(concreact, np.prod(xr, axis=1))
            #TODO: Somewhere in the math there is a 0 for concreact. Find it and fix it.

            for ind in range(len(concreact)):
                if (concreact[ind] < 0):
                    concreact[ind] = 0

            timesHappened[reac] = dkdt
            timesPossible[reac] = np.sum(concreact)
            # 5CH3 + 3H
        ##############################################################################################################################
        for j in range(len(reactRates)):  # k = (integral) dkdt
            reactRates[j] = (timesHappened[j] / timesPossible[j]) / 0.012

        for val in range(len(reactRates)):  # eliminate rare events through threshold
            if (timesHappened[val] < rarereactLimit):
                reactRates[val] = 0

        # Check for invalid reaction rates
        if (any(reactRates >= math.inf)):
            print("Error. Contains NaN k")
            sys.exit()
        if (any(reactRates > (1 / 0.012))):
            # if the time of reaction is greater than the time allocated for simulation
            print("Error. Some k values are too large.")
            sys.exit()

        return reactRates