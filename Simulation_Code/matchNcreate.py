import re
import numpy as np

class matchNcreate:
    def doStringMatch(reactionArr, moleculeArr):
        masterSM = np.zeros([len(reactionArr), len(moleculeArr)])
        reactantIndex = np.array([])
        productIndex = np.array([])

        for reac in range(len(reactionArr)):
            wholeLine = reactionArr[reac, :].split('  => ')
            reactants = wholeLine[0].split('  + ')
            products = wholeLine[1].split('  + ')

            if 'reactants' not in dir():
                reactants = np.array(reactants)
            if 'products' not in dir():
                products = np.array(products)
            # Check if the variables are numpy arrays and if not, make them

            for rtant in range(len(reactants)):
                for mole in range(len(moleculeArr)):
                    if(re.match(moleculeArr[mole], reactants[rtant])):
                        reactantIndex = np.append(reactantIndex, mole)
            masterSM[reac, reactantIndex] = masterSM[reac, reactantIndex] - 1
            # Document the indices that are reactants and subtract 1 from them in the end

            for prod in range(len(products)):
                for mole in range(len(moleculeArr)):
                    if(re.match(moleculeArr[mole], products[prod])):
                        productIndex = np.append(productIndex, mole)
            masterSM[reac, productIndex] = masterSM[reac, productIndex] + 1
            # Document the indices that are products and add 1 to them in the end

        return masterSM