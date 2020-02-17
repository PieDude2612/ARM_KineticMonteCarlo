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

        timesHappened = np.zeros(len(reactRates))
        timesPossible = np.zeros(len(reactRates))
        concreact = np.array([])
        ##############################################################################################################################
        for reac in range(numreacts):
            dkdt = np.sum(reactperFrame[:, reac])
            reactionOrderIndex = np.array([]).astype(int)
            reactionOrder = np.array([]).astype(int)

            for reactants in range(len(stoich_mat_pos[reac, :])):
                if(stoich_mat_pos[reac, reactants] > 0):
                    reactionOrderIndex = np.append(reactionOrderIndex, reactants)
                    # take nonzero value indices out for reactants
                    reactionOrder = np.append(reactionOrder, stoich_mat_pos[reac, reactants])
                    # input their reac num here.

            reactionNotReady = np.array([]).astype(int)
            reactionReady = np.array([]).astype(int)

            for timestep in range(x0.shape[0]): #find out if a reaction is ready or not
                checkReactant = 0
                for species in range(len(reactionOrderIndex)):
                    molVal = x0[timestep, reactionOrderIndex[species]]
                    reactantVal = reactionOrder[species]

                    if(molVal < reactantVal):
                        reactionNotReady = np.append(reactionNotReady, timestep)
                        checkReactant = checkReactant + 1
                        break

                if(checkReactant == 0):
                    reactionReady = np.append(reactionReady, timestep)
                    # 0 is false 1 is true for molecule having more conc than required

            xr = x0[reactionReady, ireactmat[reac, np.nonzero(ireactmat[reac, :])] - 1].astype(int)

            # get all reactants of reactions that are ready
            theOrderMatrix = np.multiply(np.ones((1, len(reactionReady))), concexp[reac, :])
            xr = np.subtract(xr, theOrderMatrix) # operands could not be broadcast together with shapes 0, 26
            concreact = np.append(concreact, np.prod(xr, axis=1))

            for ind in range(len(concreact)):
                if (concreact[ind] < 0):
                    concreact[ind] = 0

            timesHappened[reac] = timesHappened[reac] + dkdt
            timesPossible[reac] = timesPossible[reac] + np.sum(concreact)
        ##############################################################################################################################
        for j in range(len(reactRates)):  # k = (integral) dkdt
            reactRates[j] = (timesHappened[j] / timesPossible[j]) / tsteps

        for val in range(len(reactRates)):  # eliminate rare events through threshold
            if (timesHappened[val] < rarereactLimit):
                reactRates[val] = 0

        # Check for invalid reaction rates
        for k in range(len(reactRates)):  # some other error causing divide by 0
            if (reactRates[k] >= math.inf):
                print("Error. Contains NaN k")
                break
            if (sum(reactRates) > (1 / tsteps)):
                # if the time of reaction is greater than the time allocated for simulation
                print("Error. Some k values are too large.")
                break
        return reactRates
