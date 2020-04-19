import numpy as np
class matchNcreate():
    def doStringMatch(self, reactionArr, moleculeArr):
        masterSM = np.zeros((len(reactionArr), len(moleculeArr)))
        reactantIndex = np.array([]).astype(int)
        productIndex = np.array([]).astype(int)

        for reac in range(len(reactionArr)):
            wholeLine = (reactionArr[reac]).split('  => ')
            reactants = wholeLine[0].split('  + ')
            products = wholeLine[1].split('  + ')

            for rtant in range(len(reactants)):
                reactantIndex = np.append(reactantIndex, moleculeArr.index(reactants[rtant]))
            masterSM[reac, reactantIndex] = masterSM[reac, reactantIndex] - 1
            # Document the indices that are reactants and subtract 1 from them in the end

            for prod in range(len(products)):
                productIndex = np.append(productIndex, moleculeArr.index(products[prod]))
            masterSM[reac, productIndex] = masterSM[reac, productIndex] + 1
            # Document the indices that are products and add 1 to them in the end

        return masterSM