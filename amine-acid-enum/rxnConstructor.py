# Reaction Constructor for Amine and Carboxylic Acid Reaction Enumeration
# Implementation of the algorithm as published in https://doi.org/10.1038/s41586-020-2142-y

# BSD 3-Clause License
# Copyright (c) 2020-2021, Cernak Lab at the University of Michigan


from rdkit import Chem

class Amine(Chem.Mol):
    '''
    Class for Amine, based on rdkit.Mol object
    args:
        h - hybridization of the alpha and beta carbon
        X - reacting atom 
        Y - modification (optional)
        addBeta - add a second beta carbon (default to False)
    
    The amine molecule is used for enumerating and constructing 
    the amine-carboxylic acid reactions
    '''
    def __init__(self, h, X, Y = None, addBeta = False):
        # Assign hybridization num and construct mol
        if h in (3, 'sp3'):
            if addBeta: 
                super().__init__(Chem.MolFromSmarts('[C:3][C:2]([N:1])[C:4]'))
            else:
                super().__init__(Chem.MolFromSmarts('[C:3][C:2]([N:1])'))
            self.hybridization = 'sp3'
            
        elif h in (2, 'sp2'):
            if addBeta: 
                super().__init__(Chem.MolFromSmarts('[C:3]=[C:2]([N:1])[C:4]'))
            else:
                super().__init__(Chem.MolFromSmarts('[C:3]=[C:2]([N:1])'))
            self.hybridization = 'sp2'
            
        elif h in (2, 'aro'):
            if addBeta: 
                super().__init__(Chem.MolFromSmarts('[c:3][c:2]([N:1])[c:4]'))
            else:
                super().__init__(Chem.MolFromSmarts('[c:3][c:2]([N:1])'))
            self.hybridization = 'aro'
        else: 
            raise ValueError('Invalid input "{}" for hybridization!'.format(h)+\
                             ' Value of h must be 0, 2, 3, "aro", "sp2", or "sp3".')
        self.numBeta = int(addBeta)+1
        
        # Assign reacting atom
        if X in (1,2,3):
            self.reactAtom = X
        else:
            reactAtomDict = {'A':1, 'NH2':1, 'N':1, 'NH_2':1, 'alpha':2 , 'beta':3}
            try:
                self.reactAtom = reactAtomDict[X]
            except KeyError: 
                print('Invalid input "{}" for reacting atom!'.format(X)+\
                      ' Value of X must be 1, 2, 3 or "NH2", "alpha", "beta".')
        
        # Assign modification
        if Y not in ('-A', None):
            raise ValueError('Invalid input "{}" for modification!'.format(Y)+\
                             ' Only "-A" is allowed for amine.')
        else: self.modifArg = Y
            
        if Y and self.reactAtom == 1:
            raise ValueError('Conflict: "NH2" and "-A". Cannot remove "NH2" group while using "N" as reaction atom!')
        
    def GetEnumRxnProps(self):
        # Report assigned properties
        print("Amine\n"+\
              "Hybridization of C : \th = {}\n".format(self.hybridization)+\
              "Reaction Atom(s) : \tX = {}\n".format({1:'A[NH2]', 2:'alpha', 3:'beta'}[self.reactAtom])+\
              'Modification : \t\tY = {}\n'.format(self.modifArg)+\
              "Number of Beta C : \t{}\n".format(self.numBeta))
        
    def Modify(self):
        if self.modifArg:
            return Chem.DeleteSubstructs(self, Chem.MolFromSmiles('N'))
        else:
            # Return with no modification
            return self
        
    def GetSmarts(self):
        return Chem.MolToSmarts(self)
    
    def GetSmiles(self):
        return Chem.MolToSmiles(self)
    
class Carboxy(Chem.Mol):
    def __init__(self, h, X, Y = None, addBeta = False):
        '''
        Class for Carboxylic Acid, based on rdkit.Mol object
        args:
            h - hybridization of the alpha and beta carbon
            X - reacting atom 
            Y - modification (optional)
            addBeta - add a second beta carbon (default to False)

        The carboxylic acid molecule is used for enumerating and constructing 
        the amine-carboxylic acid reactions
        '''      
        # Assign hybridization num and construct mol
        if h in (3, 'sp3'):
            if addBeta: 
                super().__init__(Chem.MolFromSmarts('[C:8][C:7]([C:6](-[O:5])=[O:10])[C:9]'))
            else:
                super().__init__(Chem.MolFromSmarts('[C:8][C:7]([C:6](-[O:5])=[O:10])'))
            self.hybridization = 'sp3'
            
        elif h in (2, 'sp2'):
            if addBeta: 
                super().__init__(Chem.MolFromSmarts('[C:8]=[C:7]([C:6](-[O:5])=[O:10])[C:9]'))
            else:
                super().__init__(Chem.MolFromSmarts('[C:8]=[C:7]([C:6](-[O:5])=[O:10])'))
            self.hybridization = 'sp2'
            
        elif h in (0, 'aro'):
            if addBeta: 
                super().__init__(Chem.MolFromSmarts('[c:8][c:7]([C:6](-[O:5])=[O:10])[c:9]'))
            else:
                super().__init__(Chem.MolFromSmarts('[c:8][c:7]([C:6](-[O:5])=[O:10])'))
            self.hybridization = 'aro'
            
        else: 
            raise ValueError('Invalid input "{}" for hybridization!'.format(h)+\
                             ' Value of h must be 0, 2, 3, "aro", "sp2", or "sp3".')
        self.numBeta = int(addBeta)+1
        
        # Assign reacting atom
        if X in (1,2,3,4):
            self.reactAtom = X + 4
        else:
            reactAtomDict = {'B[O]':5, 'O':5, 'OH':5, 'B[C]':6, 'C':6, 'alpha':7 , 'beta':8}
            try:
                self.reactAtom = reactAtomDict[X]
            except KeyError: 
                print('Invalid input "{}" for reacting atom!'.format(X)+\
                      ' Value of X must be 1, 2, 3 or "NH2", "alpha", "beta".')
        
        # Assign modification
        if Y not in ('-B', '-OH', '+H', '+H2', None):
            raise ValueError('Invalid input "{}" for modification!'.format(Y)+\
                             ' Only "-B", "-OH", "+H", "+H2" are allowed for carboxylic acid.')
        else: self.modifArg = Y
            
        if self.reactAtom == 5:
            if Y == "-OH":
                raise ValueError('Conflict: "B[O]" with "-OH". '+\
                                 'Cannot remove "OH" group while using "O" as reaction atom!')
            if Y == "-B":
                raise ValueError('Conflict: "B[O]" with "-B". '+\
                                 'Cannot remove "COOH" group while using "O" as reaction atom!')
        if self.reactAtom == 6 and Y == "-B":
            raise ValueError('Conflict: "B[C]" with "-B". '+\
                             'Cannot remove "COOH" group while using "B[C]" as reaction atom!')
                
    def GetEnumRxnProps(self):
        # Report assigned properties
        print("Carboxylic Acid\n"+\
              "Hybridization of C : \th = {}\n".format(self.hybridization)+\
              "Reaction Atom(s) : \tX = {}\n".format({5:'B[O]', 6:'B[C]', 7:'alpha', 8:'beta'}[self.reactAtom])+\
              'Modification : \t\tY = {}\n'.format(self.modifArg)+\
              "Number of Beta C : \t{}\n".format(self.numBeta))
        
        
    def Modify(self):
        # Construct an intermediate form
        self.interm = self
        
        if self.modifArg == "-B":
        # drop the carboxylic group
            self.interm = Chem.DeleteSubstructs(self, Chem.MolFromSmiles('C(=O)O'))
        
        # Evaluate the rest modification arguments
        elif self.modifArg == "-OH":
            self.interm = Chem.ReplaceSubstructs(self, 
                                             Chem.MolFromSmiles('C(=O)O'), 
                                             Chem.MolFromSmarts('[C:6]=[O:10]'),
                                             replaceAll=True)[0]
        
        elif self.modifArg == "+H":
            if self.reactAtom == 5:
                # keep the reacting oxygen
                self.interm = Chem.ReplaceSubstructs(self, 
                                             Chem.MolFromSmiles('C(=O)O'), 
                                             Chem.MolFromSmarts('[C:6](-[O:5])[O:10]'),
                                             replaceAll=True)[0]
            
            else:
                # remove alcohol
                self.interm = Chem.ReplaceSubstructs(self, 
                                                 Chem.MolFromSmarts('C(=O)O'), 
                                                 Chem.MolFromSmarts('[C:6][O:10]'),
                                                 replaceAll=True)[0]
            
        elif self.modifArg == "+H2":
            if self.reactAtom == 5:
                # keep the reacting oxygen
                self.interm = Chem.ReplaceSubstructs(self, 
                                             Chem.MolFromSmiles('C(=O)O'), 
                                             Chem.MolFromSmarts('[C:6][O:5]'),
                                             replaceAll=True)[0]
            else:
                # fully reduce all oxygens
                self.interm = Chem.ReplaceSubstructs(self, 
                                             Chem.MolFromSmarts('C(=O)O'), 
                                             Chem.MolFromSmarts('[C:6]'),
                                             replaceAll=True)[0]
                
        elif self.modifArg == None and self.reactAtom == 6:
            # remove the alcohol if reaction atom is B[C]
            self.interm = Chem.ReplaceSubstructs(self, 
                                             Chem.MolFromSmiles('C(=O)O'), 
                                             Chem.MolFromSmarts('[C:6]=[O:10]'),
                                             replaceAll=True)[0]
                
        # return self if no modification on other reaction atoms
        return self.interm
    
    def GetSmarts(self):
        return Chem.MolToSmarts(self)
    
    def GetSmiles(self):
        return Chem.MolToSmiles(self)
    
class Product(Chem.Mol):
    '''
    Product class for the amine-acid reaction, based on rdkit.Chem.Mol
    args:
        - Amine Mol object
        - Carboxy Mol object
    Construct the product based on the modification of the reagents and form 
    the bond based on the indicated reaction atoms
    '''
    def __init__(self, amine, carboxy):
        # Construct a editable Mol from modified Amine and Carboxylic acid
        edMol = Chem.EditableMol(Chem.CombineMols(amine.Modify(),carboxy.Modify()))
        
        atomIdxDict = {}
        # Match atom list indice to atom mapping numbers
        for i, atm in enumerate(edMol.GetMol().GetAtoms()):
            atomIdxDict[atm.GetAtomMapNum()] = i
            
        # Find attaching atoms by mapping number on each constituent mol and make a single bond
        edMol.AddBond(atomIdxDict[amine.reactAtom], 
                     atomIdxDict[carboxy.reactAtom],
                     order=Chem.rdchem.BondType.SINGLE)
        #replace all single bond in the product as an arbitrary single bond (allow aromaticity)
        self.sma = Chem.MolToSmarts(edMol.GetMol()).replace("-,:", '')#.replace("-", '')
        # Construct the final product
        super().__init__(Chem.MolFromSmarts(self.sma))
        hA = "{"+amine.hybridization+"}"
        if amine.modifArg:
            aReac = "{"+{1:'A', 2:'\\alpha', 3:'\\beta'}[amine.reactAtom]+amine.modifArg+"}"
        else: 
            aReac = "{"+{1:'A', 2:'\\alpha', 3:'\\beta'}[amine.reactAtom]+"}"
        
        hB = "{"+carboxy.hybridization+"}"
        
        if carboxy.modifArg:
            bReac = "{"+{5:'B[O]', 6:'B[C]', 7:'\\alpha', 8:'\\beta'}[carboxy.reactAtom]+carboxy.modifArg+"}"
        else:
            bReac = "{"+{5:'B[O]', 6:'B[C]', 7:'\\alpha', 8:'\\beta'}[carboxy.reactAtom]+"}"
        
        latexModStr = '$^{}{}^{}/^{}{}^{}$'.format(hA, '{NH_2}', aReac,
                                                   hB, '{COOH}', bReac)
        self.SetProp('Modification', latexModStr)
        
        
    def GetSmarts(self):
        return self.sma
    
    def GetSmiles(self):
        return Chem.MolToSmiles(self)
    
class RxnConstructor:
    
    '''
    Class for constructing and storing objects for the amine-acid reaction
    The reaction object is based on rdkit.AllChem.ChemicalReaction
    args:
        - Amine object args (h, X, Y, addBeta) as a list or a tuple
        - Carboxy object args (h, X, Y, addBeta) as a list or a tuple
        
    Construct the reagents and product based upon initialization and form the reaction template
    '''
    
    def __init__(self, amineRxnProps, carboxyRxnProps):
        # Construct reagents and product from the reaction prop tuple
        self.amine = Amine(*amineRxnProps)
        self.carboxy = Carboxy(*carboxyRxnProps)
        self.prod = Product(self.amine, self.carboxy)
        self.rxnSmarts = " ".join([self.amine.GetSmarts(),'.',
                                   self.carboxy.GetSmarts(),'>>',
                                   self.prod.GetSmarts()])
        
        self.reverseSmarts = " ".join([self.prod.GetSmarts(),'>>',
                                       self.amine.GetSmarts(),'.',
                                       self.carboxy.GetSmarts()])
        # Construct the rxn
        self.rxn = Chem.rdChemReactions.ReactionFromSmarts(self.rxnSmarts)
        self.reverse = Chem.rdChemReactions.ReactionFromSmarts(self.reverseSmarts)
    
    def ReportRxnProps(self):
        self.amine.GetEnumRxnProps()
        self.carboxy.GetEnumRxnProps()
        print("RXN Smarts:", self.rxnSmarts)
        
def enumAllRxns(# Args default to all the possible properties
                # that enumerate 720 rxns
                
                # hybridization
                hA = ('sp3','sp2','aro'),
                hB = ('sp3','sp2','aro'),

                # reaction atom
                XA = ("NH2", "alpha", "beta"),
                XB = ("B[O]", "B[C]", "alpha", "beta"),

                # modification
                YA = (None, '-A'),
                YB = (None, '-B', '-OH', '+H', '+H2'),
                
                # add beta carbon
                addBetaA = (False, True),
                addBetaB = (False, True),
                
                ):
    
    '''
    Function for enumerating all the possible reactions between the amine and carboxylic acid
    based on all the available rxn property args.
    Generates all 720 rxns and their rxn property args
    '''
    
    # Create enumeration properties for both reagents
    aArgs = []
    for x in hA:
        for y in XA:
            for z in YA:
                for add in addBetaA:
                    # exclude combinations containing NH2 with -A
                    if 'NH2' in (x,y,z) and '-A' in (x,y,z):
                        pass
                    else:
                        aArgs.append((x,y,z,add))

    bArgs = []
    for x in hB:
        for y in XB:
            for z in YB:
                for add in addBetaB:
                    # exclude combinations containing B[O] with -OH or -B
                    if 'B[O]' in (x,y,z) and ('-OH'  in (x,y,z) or '-B' in (x,y,z)):
                        pass
                    # exclude combinations containing B[C] with -OH or -B
                    elif 'B[C]' in (x,y,z) and ('-OH'  in (x,y,z) or '-B' in (x,y,z)):
                        pass
                    else:
                        bArgs.append((x,y,z,add))

    print('Total enumerated rxns: {} * {} = {}'.format(len(aArgs),len(bArgs),len(aArgs)*len(bArgs)))

    # Construct Reactions
    rxns = []
    rxnConds = []
    for a in aArgs:
        for b in bArgs:
            # Append a tuple of properties expanded from both reactant
            # For creating a dataframe later
            rxnConds.append((*a,*b))
            # Construct the reaction 
            rxnConst = RxnConstructor(a, b)
            rxns.append(rxnConst)

    print(len(rxns), "reactions created")
    
    return rxns, rxnConds

if __name__ == "__main__":
    
    from pandas import DataFrame
    
    rxns, rxnConds = enumAllRxns()
    df = DataFrame(rxnConds,columns=['hybridization A', 'Rxn Atom A', 'Modification A', 'hybridization B', 'Rxn Atom B', 'Modification B'])
    df['Reactions Smarts'] = [rxnConst.rxnSmarts for rxnConst in rxns]
    df.to_csv('enumRxns720.csv', index = False)
    print('Data saved at enumRxns720.csv')