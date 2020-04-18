import numpy as np
import re
class matchNcreate():
    def doStringMatch(self, reactionArr, moleculeArr):
        masterSM = np.zeros((len(reactionArr), len(moleculeArr)))
        reactantIndex = np.array([]).astype(int)
        productIndex = np.array([]).astype(int)

        for reac in range(len(reactionArr)):
            wholeLine = np.array((reactionArr[reac]).split('  => '))
            reactants = np.array(wholeLine[0].split('  + ')).astype(str)
            products = np.array(wholeLine[1].split('  + ')).astype(str)

            for rtant in range(len(reactants)):
                for mole in range(len(moleculeArr)):
                    if(moleculeArr[mole] == reactants[rtant]):
                        reactantIndex = np.append(reactants, mole)
            masterSM[reac, reactantIndex] = masterSM[reac, reactantIndex] - 1
            # Document the indices that are reactants and subtract 1 from them in the end

            for prod in range(len(products)):
                for mole in range(len(moleculeArr)):
                    if(moleculeArr[mole] == products[prod]):
                        productIndex = np.append(productIndex, mole)
            masterSM[reac, productIndex] = masterSM[reac, productIndex] + 1
            # Document the indices that are products and add 1 to them in the end

        return masterSM