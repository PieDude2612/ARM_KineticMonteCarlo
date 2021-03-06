import numpy as np

class propensity_fcn():

    def calc_propensity(self, stoichmat, xarr, rp):
        # determine if reaction passed is uni, bi, or tri molecular.
        # determine if it is a forward or reverse reaction.

        propArr = np.double([])
        numRows = stoichmat.shape[0]
        for row in range(numRows):
            propArr = np.append(propArr, propensity_fcn.findSectionPropensity(propensity_fcn(), stoichmat[row, :],
                                                                              xarr, rp[row]))

        return propArr

    def findSectionPropensity(self, stoich, x, rateparam):
        negatives = np.array([])
        val = np.array([])
        i = 0

        # iterate through the passed vector to find any negatives (reactants) or
        # positives (products). i is a dummy variable

        while (i < len(stoich)):
            if (stoich[i] < 0):
                negatives = np.append(negatives, i)
                val = np.append(val, np.absolute(stoich[i]))
                i = i + 1
                continue
            else:
                i = i + 1
                continue

        hXt = 1
        k = rateparam

        for num in range(len(negatives)):  # any kind of reaction between diff species
            n = int(negatives[num])
            v = int(val[num])

            if (x[n] < 0):
                return 0

            if (v > 1):  # more than one of same species
                c = v
                while (c > 0):
                    hXt = hXt * (x[n] + (c - v))
                    c = c - 1
                continue
            else:  # different species
                hXt = hXt * x[n]
                continue

        if ((hXt * k) == k):  # improbable reaction
            return 0

        reaction_prob = hXt * k
        return reaction_prob
