# Reaction Constructor for C-H Functionalization Reaction Enumeration

# BSD 3-Clause License
# Copyright (c) 2020-2021, Cernak Lab at the University of Michigan

from rdkit import Chem
from rdkit.Chem import Draw
from copy import copy
import itertools as itt
from reactants import Carbon, Amine, Alcohol, Carboxyl, Bromide, Boronate

# Suppress rdkit waring messeges
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

class Product:
    '''
    Product class for the C-H functionalization reaction, contains rdkit.Chem.Mol
    args:
        - mA: Carbon Mol object
        - mB: Carbon, Amine, or Alcohol Mol object
        - AMapNum: Tuple indicating reaction atom(s) of mA
        - BMapNum: Tuple indicating reaction atom(s) of mB
    Construct the product based on the modification of the reagents and form 
    the bond based on the indicated reaction atoms
    '''
    def __init__(self, mA, mB, AMapNum, BMapNum):
        
        # Construct a editable Mol from modified reactA and reactB
        edMol = Chem.EditableMol(Chem.CombineMols(mA, mB))

        atomIdxDict = {}
        # Match atom list indice to atom mapping numbers
        for i, atm in enumerate(edMol.GetMol().GetAtoms()):
            atomIdxDict[atm.GetAtomMapNum()] = i

        # Find attaching atoms by mapping number on each constituent mol and make a bond
        if len(AMapNum) == len(BMapNum) == 1:
            for i,j in zip(AMapNum, BMapNum):
                idxA = atomIdxDict[i]
                idxB = atomIdxDict[j]
                BondOrder = Chem.rdchem.BondType.SINGLE
                # Make it single to ensure the formation of biaryl
                edMol.AddBond(idxA, idxB, order=BondOrder)
                
        elif len(AMapNum) == len(BMapNum) == 2:
            for i,j in zip(AMapNum, BMapNum):
                idxA = atomIdxDict[i]
                idxB = atomIdxDict[j]
                BondOrder = self.DetermineBO(edMol.GetMol(), idxA, idxB)
                edMol.AddBond(idxA, idxB, order=BondOrder)

        else:
            for i in AMapNum:
                for j in BMapNum:
                    # For bondings of one atom from A to two atoms from B to form a ring
                        idxA = atomIdxDict[i]
                        idxB = atomIdxDict[j]
                        BondOrder = self.DetermineBO(edMol.GetMol(), idxA, idxB)
                        edMol.AddBond(idxA, idxB, order=BondOrder)
            
        self.mol = edMol.GetMol()
        self.sma = Chem.MolToSmarts(self.mol)
        
    @staticmethod
    def DetermineBO(mol, atom_idx_i, atom_idx_j):
        aro_i = mol.GetAtoms()[atom_idx_i].GetIsAromatic()
        aro_j = mol.GetAtoms()[atom_idx_j].GetIsAromatic()
        if aro_i and aro_j:
            BO = Chem.rdchem.BondType.AROMATIC
        else:
            BO = Chem.rdchem.BondType.SINGLE
        return BO
        
    def _repr_svg_(self):
        return Draw._MolsToGridSVG([self.mol])
    
    # TODO: Create LaTex strings for the reaction identifiers
    '''
    def GetLatex(self):
        hA = "{"+reactA.hybridization+"}"
        if reactA.modifArg:
            aReac = "{"+{1:'A', 2:'\\alpha', 3:'\\beta'}[reactA.reactAtom]+reactA.modifArg+"}"
        else: 
            aReac = "{"+{1:'A', 2:'\\alpha', 3:'\\beta'}[reactA.reactAtom]+"}"
        
        hB = "{"+reactB.hybridization+"}"
        
        if reactB.modifArg:
            bReac = "{"+{5:'B[O]', 6:'B[C]', 7:'\\alpha', 8:'\\beta'}[reactB.reactAtom]+reactB.modifArg+"}"
        else:
            bReac = "{"+{5:'B[O]', 6:'B[C]', 7:'\\alpha', 8:'\\beta'}[reactB.reactAtom]+"}"
        
        latexModStr = '$^{}{}^{}/^{}{}^{}$'.format(hA, '{NH_2}', aReac,
                                                   hB, '{COOH}', bReac)
        self.SetProp('Modification', latexModStr)
    '''
    
                
class Flask:
    '''
    A namespace for holding rxn, mols, and Smarts
    '''
    def __init__(self, smaA, smaB, mProd, rxnSmarts, reverseSmarts):
        self.A = smaA
        self.B = smaB
        self.prod = mProd
        
        self.rxnSmarts = rxnSmarts
        self.reverseSmarts = reverseSmarts
                                   
        self.rxn = Chem.rdChemReactions.ReactionFromSmarts(self.rxnSmarts)
        self.rrxn = Chem.rdChemReactions.ReactionFromSmarts(self.reverseSmarts)
        
    def _repr_svg_(self):
        return Draw.ReactionToImage(self.rxn, useSVG=True)
        
class GetRxn():
    
    '''
    Class for constructing and storing objects for the C-H functionalization reactions
    The reaction object is based on rdkit.AllChem.ChemicalReaction
    args:
        - reactA class (Carbon(h, X)) args should be the class object itself
        - reactB class (Carbon, Amine, or Alcohol) as class object itself
        
    Construct the reagents and product upon initialization and form the reaction template
    '''
    
    def __init__(self, reactA, reactB):
        # Check if reactant A is a carbon class
        #if reactA.__class__.__name__ != "Carbon":
            #raise TypeError("Reactant A has to be Carbon/Carbon Chain")
            
        # Check if both reactants are Carbon
        if reactA.__class__.__name__ == reactB.__class__.__name__ == "Carbon":
            # Reconstruct the second Carbon for different carbon mapping numbers
            reactB = copy(reactB)
            reactB.RemapAtoms()
            
        self.reactA = reactA
        self.reactB = reactB
        
        smaPairs = list(itt.product(reactA.smas,
                                          reactB.smas))
        
        self.rxns = []
        unique_rxns = []
        # Construct reactants and product from the A B mols
        for smaA, smaB in smaPairs:
            
            modAs = [Chem.MolFromSmarts(sma) for sma in reactA.ModifySmarts([smaA])]
            modBs = [Chem.MolFromSmarts(sma) for sma in reactB.ModifySmarts([smaB])]
            
            molPairs = list(itt.product(modAs, modBs))
            
            for mA, mB in molPairs:
                prod = Product(mA, mB, reactA.reactAtom, reactB.reactAtom)
                prod.mol = self.CorrectAromaticity(prod.mol)

                if smaA == '[c:1]'  or smaA == '[C:1]':
                # Save the single aromatic C template from being sanitized  
                # for C-C coupling on any arbituary positon of aromatic ring.
                # Otherwise, this template will be sanitized and will miss finding many possible rxns.
                    pass

                elif smaB == '[c:1]' or smaB == '[C:1]':
                # Same for the second carbon
                    pass

                else:
                    try:
                        Chem.SanitizeMol(prod.mol)
                        prod.sma = Chem.MolToSmarts(prod.mol)
                    except:
                        continue

                rxnSmarts = " ".join([smaA,'.',smaB,'>>',prod.sma])
                reverseSmarts = " ".join([prod.sma,'>>',smaA,'.',smaB])

                if rxnSmarts not in unique_rxns:
                # Save everything
                    unique_rxns.append(rxnSmarts)
                    self.rxns.append(Flask(smaA, smaB, prod.mol, rxnSmarts, reverseSmarts))
                
    def __iter__(self):  
        return iter(self.rxns)
    
    def __getitem__(self, x):
        return self.rxns[x]
    
    def __len__(self):
        return len(self.rxns)
    
    def _ipython_display_(self):
        from IPython.display import display
        self.reactA.GetEnumRxnProps()
        self.reactB.GetEnumRxnProps()
        for i,rxn in enumerate(self.rxns):
            print(i,'\t',rxn.rxnSmarts)
            display(rxn)
    
    def CorrectAromaticity(self, mol):
        # Reconstruct Mol with SMILES
        tempMol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        if tempMol:
            # Mol fraction from Smarts may not convert to valid Smiles
            bondrings = tempMol.GetRingInfo().BondRings()
            if len(bondrings) > 0:
            # First check if it's a ring
                if self.IsRingAromatic(tempMol, bondrings[0]):
                    # Save the aromatic ring instead
                    mol = tempMol 

        return mol

    @staticmethod
    def IsRingAromatic(mol, bondRing):
        ## Copied from rdkit cookbook
        for i in bondRing:
            if not mol.GetBondWithIdx(i).GetIsAromatic():
                return False
        return True
    
    def ReportRxnProps(self):
        self.reactA.GetEnumRxnProps()
        self.reactB.GetEnumRxnProps()
        for rxn in self.rxns:
            print("RXN Smarts:", rxn.rxnSmarts)
            
class Enumerator:
    def __init__(self, A, B, leastNumRingMembers = 5):
        h = ('sp3','sp2')
        # All possible X, Y options for each functional group
        ConfigDict = {'Carbon':{'X':(1, (1,2), (1,3), (2,1), (3,1))},
                      'Amine':{'X':(1, 2, 3, *itt.permutations(range(1,4),2)),
                               'Y':(None, '-B')},
                      'Alcohol':{'X':(1, 2, 3, *itt.permutations(range(1,4),2)),
                                 'Y':(None, '-B')},
                      'Carboxyl':{'X':(1,2,3,4,(1,2),(1,3),(2,3),(2,4),
                                      (3,4),(2,1),(3,1),(3,2),(4,2),(4,3)),
                                  'Y':(None, '-B', '-OH', '+H', '+H2')},
                      'Boronate':{'X':(1,2,3,4,5, *itt.permutations(range(2,6),2)),
                                  'Y':(None, '-C')},
                      'Bromide':{'X':(2, 3, (2,3), (3,2)),
                                 'Y':(None, '-B')},
                  }
        if A in ConfigDict.keys() and B in ConfigDict.keys():
            self.ClassA = eval(A)
            self.ClassB = eval(B)
        else:
            raise ValueError('Invalid functional group!\n\
            Please select from {}'.format(tuple(ConfigDict.keys())))
            
        # Generate all possible rxn configurations for each functional group
        if A == 'Carbon':
            AConfigs = list(itt.product(h,ConfigDict['Carbon']['X']))
        else:
            AConfigs = list(itt.product(h,
                                        ConfigDict[A]['X'],
                                        ConfigDict[A]['Y']))
        if B == 'Carbon':
            BConfigs = list(itt.product(h,ConfigDict['Carbon']['X']))
        else:
            BConfigs = list(itt.product(h,
                                        ConfigDict[B]['X'],
                                        ConfigDict[B]['Y']))
            
        self.AConfigs = self.ApplyConstraints(A, AConfigs)
        self.BConfigs = self.ApplyConstraints(B, BConfigs)
        
        self.rxnConfigs = list(itt.product(self.AConfigs,self.BConfigs))
        
        
        if leastNumRingMembers:
            self.rxnConfigs = self.ReduceRingMembers(self.rxnConfigs, leastNumRingMembers)
        
    def ApplyConstraints(self, funcGroupName, configs):
        if funcGroupName in ('Amine', 'Alcohol'):
            for x in copy(configs):
                # Avoid X = 1 and Y = -B
                if '-B' in x:
                    if 1 in x or (type(x[1]) is tuple and 1 in x[1]):
                        configs.remove(x)
                
        elif funcGroupName == 'Boronate':
            for x in copy(configs):
                # Avoid X = 2 or 3 and Y = -C
                if '-C' in x:
                    if 2 in x or (type(x[1]) is tuple and 2 in x[1]):
                        configs.remove(x)
                    elif 3 in x or (type(x[1]) is tuple and 3 in x[1]):
                        configs.remove(x)
                        
        elif funcGroupName == 'Carboxyl':
            for x in copy(configs):
                if '-OH' in x:
                    if (4 in x) or (type(x[1]) is tuple and 4 in x[1]):
                        configs.remove(x)
                if '-B' in x:
                    # Avoid X = 1 or 4 and Y = -B
                    if (1 in x) or (type(x[1]) is tuple and 1 in x[1]):
                        configs.remove(x)
                    if (4 in x) or (type(x[1]) is tuple and 4 in x[1]):
                        configs.remove(x)
                        
        return configs
        
    def ReduceRingMembers(self, rxnConfigs, leastMembers = 5):
        qualifiedConfigs = []
        for t in rxnConfigs:
            _XA = t[0][1]
            _XB = t[1][1]
            # exclude rings smaller than the specified number of members
            if type(_XA) != type(_XB):
                continue
            if type(_XA) == type(_XB) == tuple and (abs(_XA[0]-_XA[1])+abs(_XB[0]-_XB[1]) < leastMembers-2):
                continue
                
            qualifiedConfigs.append(t)
            
        return qualifiedConfigs
    
    def GetAllRxns(self):
        allRxns = []
        for A, B in self.rxnConfigs:
            rxns = GetRxn(self.ClassA(*A), self.ClassB(*B)).rxns
            allRxns.append((A,B,rxns))
            
        return allRxns
    
    def MakeDataFrame(self):
        import pandas as pd
        all_entries = []
        for A, B in self.rxnConfigs:
            for rxn in GetRxn(self.ClassA(*A), self.ClassB(*B)):
                all_entries.append(((A, B), rxn.A, rxn.B, rxn.rxnSmarts, rxn.prod, rxn.rxn))

        return pd.DataFrame(all_entries, columns=['A B Configs','A','B','rxn SMARTS', 
                                                  'Product', 'rxn'])
                                                  