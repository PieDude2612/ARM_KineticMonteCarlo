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
            tally = 0
            dkdt = np.array([])

            for i in range(reactperFrame.shape[0]):
                if (reactperFrame[i, 0] == reac):
                    tally = tally + 1
                    continue  # how many times it occurs

                if (tally == 0):  # if reaction given does not occur
                    break
                else:
                    dkdt = np.append(dkdt, tally)
                    continue

            reactionOrderIndex = np.nonzero(stoich_mat_pos[reac, :])
            # take nonzero value indices out for reactants
            reactionOrder = stoich_mat_pos[reac, reactionOrderIndex[0:len(reactionOrderIndex)]]
            # input their reac num here.

            # reacOrderMatrix = np.multiply(np.ones(1, (tsteps - backt)), reactionOrder)  # matrix of reactant concs

            # reactionNotReady = x0[startFrame:endFrame, reactionOrderIndex[0:len(reactionOrderIndex)]] < reacOrderMatrix[0:reacOrderMatrix.shape[0], 0:reacOrderMatrix.shape[1]] # 1D array of 0 and 1
            reactionNotReady = np.array([])

            for timestep in range(x0.shape[0]): #find out if a reaction is ready or not
                for species in range(len(reactionOrderIndex)):
                    checkConcCol = np.bool([])

                    molVal = x0[timestep, reactionOrderIndex[species]]
                    reactantVal = reactionOrder[species]
                    checkConcCol = np.append(checkConcCol, molVal < reactantVal)
                    # 0 is false 1 is true for molecule having more conc than required
                if (np.sum(checkConcCol) > 0):
                    reactionNotReady = np.append(reactionNotReady, 1)
                    continue

            yesreact = np.array([])
            reactionReady = np.where(reactionNotReady[0:len(reactionNotReady)] == 0)
            yesreact = np.append(yesreact, reactionReady)

            xr = np.array([]) #do the math to find the rates
            xr_col = np.array([])
            for timestepWhereReady in range(len(reactionReady)):
                for speciesAtTimestep in range(len(ireactmat[reac, :])):
                    xr_col = np.append(xr_col, x0[reactionReady[timestepWhereReady], ireactmat[reac, speciesAtTimestep]])
                xr = np.append(xr, xr_col, axis=0)

            # get all reactants of reactions that are ready
            print(concexp[reac, 0:np.sum(reactionOrder)])
            print(len(yesreact))
            xr = xr - (np.dot(np.ones([1, len(yesreact)]), concexp[reac, 0:np.sum(reactionOrder)]))
            concreact = np.prod(xr, axis=1)

            np.delete(yesreact, 0, 0)
            for ind in range(len(concreact)):
                if (concreact[ind] < 0):
                    concreact[ind] = 0

            timesHappened[reac] = timesHappened[reac] + dkdt[yesreact]
            lenyesr[reac] = lenyesr[reac] + len(yesreact)
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
            if (sum(reactRates) > (
                    1 / tsteps)):  # if the time of reaction is greater than the time allocated for simulation
                print("Error. Some k values are too large.")
                break
        return reactRates
