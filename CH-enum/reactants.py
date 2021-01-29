# Reactant classes for C-H Functionalization Reaction Enumeration

# BSD 3-Clause License
# Copyright (c) 2020-2021, Cernak Lab at the University of Michigan

from rdkit import Chem
from rdkit.Chem import Draw
import itertools

class Carbon:
    '''
    Namespace for Carbon, constructs SMARTS and rdkit.Mol objects
    args:
        h - hybridization of the alpha and beta carbon
        X - reacting atom 
        
    The Carbon molecules created by this class is used for 
    enumerating and constructing the C-H functionalization reactions
    '''
    def __init__(self, h, X):
        # Assign hybridization num and construct mol
        if h in (3, 'sp3'):
            self.smas = ['[C:1]']
            self.hybridization = 'sp3'
            
        elif h in (2, 'sp2'):
            self.smas = ['[C:1]=[*:2]', '[c:1]']
            self.hybridization = 'sp2'
            
        else: 
            raise ValueError('Invalid input "{}" for hybridization!'.format(h)+\
                             ' Value of h must be 2, 3, "sp2", or "sp3".')
        
        # Assign reacting atom
        # X has to be 1, (1,2) or (1,3)
        if X == 1 or X == (1,): 
            self.reactAtom = (1,)
        elif X in ((1, 2), (2, 1)):
            self.reactAtom = X
            if self.hybridization == 'sp2':
                self.smas = ['[C:1]=[C:2]', '[c:1][c:2]']
            else:  # sp3
                self.smas = ['[C:1]-[C:2]','[C:1]-[c:2]']
        elif X in ((1, 3), (3, 1)):
            self.reactAtom = X
            if self.hybridization == 'sp2':
                self.smas = ['[C:1]=[C:2]-[C:3]', '[c:1][c:2][c:3]']
            else: # sp3
                self.smas = ['[C:1]-[C:2]-[C:3]','[C:1]-[c:2][c:3]']
        else:
            raise ValueError('Invalid input "{}" for reacting atom!'.format(X)+\
                  ' Value of X must be 1, (1,2) or (1,3).')
            
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
       
    def __copy__(self):
        return self.__class__(self.hybridization, self.reactAtom)
    
    def _repr_svg_(self):
        labels = [str(i+1) for i in range(len(self.smas))]
        return Draw._MolsToGridSVG(self.mols,legends=labels)
    
    def GetEnumRxnProps(self):
        # Report assigned properties
        print(self.__class__.__name__+'\n'+\
              "Hybridization of C1 : \th = {}\n".format(self.hybridization)+\
              "Reaction Atom(s) : \tX = {}\n".format(self.reactAtom)
             )
        
    def ModifySmarts(self,smas):
        # No modification allowed
        return smas
    
    def Modify(self):
        # No modification allowed
        return self

    def RemapAtoms(self):
        # remap with prefix 9 to avoid conflict with the first carbon
        self.smas = [sma.replace('1','91').replace('2','92').replace('3','93')
                    for sma in self.smas]
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
        
        if self.reactAtom == (1,):
            self.reactAtom = (91,)
        else: 
            self.reactAtom = (self.reactAtom[0]+90, self.reactAtom[1]+90)
        return self
        
class Amine:
    '''
    Namespace for Amine, constructs SMARTS and rdkit.Mol objects
    args:
        h - hybridization of the alpha and beta carbon
        X - reacting atom 
        Y - modification (optional)
    
    The amine molecule created by this class is used for 
    enumerating and constructing the C-H Amine functionalization reactions
    '''
    def __init__(self, h, X, Y = None):
        # Assign hybridization num and construct SMARTS
        if h in (3, 'sp3'):
            self.smas = ['[C:13]-[C:12]-[N:11]']
            self.hybridization = 'sp3'
            
        elif h in (2, 'sp2'):
            self.smas = ['[C:13]=[C:12]-[N:11]','[c:13][c:12]-[N:11]']
            self.hybridization = 'sp2'
            
        else: 
            raise ValueError('Invalid input "{}" for hybridization!'.format(h)+\
                             ' Value of h must be 2, 3, "sp2", or "sp3".')
        
        # Assign reacting atom
        if X in (1, 2, 3, (1,2), (1,3), (2,3), (2,1), (3,1), (3,2)):
            # Add 10 to each number to match atom mapping number 
            # and avoid conflicts with other class mols when adding bonds
            if type(X) is int:
                self.reactAtom = (X+10,)
            else:
                self.reactAtom = (X[0]+10, X[1]+10)
        else:
            # Check if they are string inputs
            reactAtomDict = {'NH2':11, 'N':11, 'NH':11, 'alpha':12, 'beta':13}
            
            try:
                if type(X) is str:
                    self.reactAtom = (reactAtomDict[X],)
                else:
                    self.reactAtom = (reactAtomDict[X[0]], 
                                      reactAtomDict[X[1]])
            except (KeyError, TypeError): 
                raise ValueError(
                    'Invalid input "{}" for reacting atom!\n'.format(X)+\
                    'Value of X must be numerical values: 1, 2, 3, (1, 2), (1, 3), (2, 3)\n'+\
                    'Or strings: "N", "alpha", "beta".\n'+\
                    'Revering tuple order is allowed'
                )
        
        # Assign modification
        if Y not in ('-N', '-B', None):
            raise ValueError('Invalid input "{}" for modification!'.format(Y)+\
                             ' Only "-N" or "-B" is allowed for amine.')
            
        else: self.modifArg = Y
        
        # Checking args conflict
        if self.modifArg and 11 in self.reactAtom:
            raise ValueError('Conflict: X = "N" and Y = "-N". Cannot remove "N" while using it as reaction atom!')
            
        # Construct the molecule
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
    
    def _repr_svg_(self):
        labels = [str(i+1) for i in range(len(self.smas))]
        return Draw._MolsToGridSVG(self.mols,legends=labels)
    
    def GetEnumRxnProps(self):
        # Report assigned properties
        print(self.__class__.__name__+'\n'+\
              "Hybridization of C : \th = {}\n".format(self.hybridization)+\
              "Reaction Atom(s) : \tX = {}\n".format(self.reactAtom)+\
              'Modification : \t\tY = {}\n'.format(self.modifArg)
              )
        
    def ModifySmarts(self, smalist):
        smalist = list(smalist)
        # This is needed for separate treatments on individual SMARTS later
        if self.modifArg in ('-N', '-B'):
            smalist = [sma.replace('-[N:11]','') for sma in smalist]
            
        arolist = [sma.replace('[C:12]-[N:11]','[C:12]-[n:11]')\
                   .replace('[c:12]-[N:11]','[c:12][n:11]')
                   for sma in smalist]
        smalist.extend(arolist)
        smalist = list(set(smalist))
        return smalist
        
    def Modify(self):
        self.smas = self.ModifySmarts(self.smas)
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
        return self
        
class Carboxyl:
    def __init__(self, h, X, Y = None):
        '''
        Namespace for Carboxylic Acid, construct SMARTS and rdkit.Mol objects
        args:
            h - hybridization of the alpha and beta carbon
            X - reacting atom 
            Y - modification (optional)

        The carboxylic acid molecule created by this class is used for 
        enumerating and constructing the C-H functionalization reactions
        '''      
        # Assign hybridization num and construct mol
        if h in (3, 'sp3'):
            self.smas = ['[C:23]-[C:22]-[C:21](-[O:24])=[O:25]']
            self.hybridization = 'sp3'
            
        elif h in (2, 'sp2'):
            self.smas = ['[C:23]=[C:22]-[C:21](-[O:24])=[O:25]','[c:23][c:22]-[C:21](-[O:24])=[O:25]']
            self.hybridization = 'sp2'
            
        else: 
            raise ValueError('Invalid input "{}" for hybridization!'.format(h)+\
                             ' Value of h must be 2, 3, "sp2", or "sp3".')
        
        # Assign reacting atom
        if X in (1,2,3,4,(1,2),(1,3),(2,3),(2,4),(3,4),(2,1),(3,1),(3,2),(4,2),(4,3)):
            if type(X) is int:
                self.reactAtom = (X+20,)
            else:
                self.reactAtom = (X[0]+20, X[1]+20)
                
        else:
            reactAtomDict = {'B[C]':21, 'C':21, 'alpha':22 , 'beta':23, 
                             'B[O]':24, 'O':24, 'OH':24}
            try:
                if type(X) is str:
                    self.reactAtom = (reactAtomDict[X],)
                else:
                    self.reactAtom = (reactAtomDict[X[0]], 
                                      reactAtomDict[X[1]])
            except (KeyError, TypeError): 
                raise ValueError(
                    'Invalid input "{}" for reacting atom!'.format(X)+\
                    ' Value of X must be 1, 2, 3, 4 or "B[C]", "alpha", "beta", "B[O]".'
                )
        
        # Assign modification
        if Y not in ('-B', '-OH', '+H', '+H2', None):
            raise ValueError('Invalid input "{}" for modification!'.format(Y)+\
                             ' Only "-B", "-OH", "+H", "+H2" are allowed for carboxylic acid.')
        else: self.modifArg = Y
            
        if 24 in self.reactAtom:
            if Y == "-OH":
                raise ValueError('Conflict: "B[O]" with "-OH". '+\
                                 'Cannot remove "OH" group while using "O" as reaction atom!')
            if Y == "-B":
                raise ValueError('Conflict: "B[O]" with "-B". '+\
                                 'Cannot remove "COOH" group while using "O" as reaction atom!')
                
        if 21 in self.reactAtom and Y == "-B":
            raise ValueError('Conflict: "B[C]" with "-B". '+\
                             'Cannot remove "COOH" group while using "B[C]" as reaction atom!')
        
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
        
    def _repr_svg_(self):
        labels = [str(i+1) for i in range(len(self.smas))]
        return Draw._MolsToGridSVG(self.mols,legends=labels)
               
    def GetEnumRxnProps(self):
        # Report assigned properties
        print(self.__class__.__name__+"\n"+\
              "Hybridization of C : \th = {}\n".format(self.hybridization)+\
              "Reaction Atom(s) : \tX = {}\n".format(self.reactAtom)+\
              'Modification : \t\tY = {}\n'.format(self.modifArg))
        
        
    def ModifySmarts(self, smalist):
        smalist = list(smalist)
        if self.modifArg == "-B":
            # drop the carboxylic group
            smalist = [sma.replace('-[C:21](-[O:24])=[O:25]','') for sma in smalist]
            
        # Evaluate the rest modification arguments
        elif self.modifArg == "-OH":
            smalist = [sma.replace('[C:21](-[O:24])=[O:25]','[C:21]=[O:25]') for sma in smalist]
        
        elif self.modifArg == "+H":
            if 24 in self.reactAtom:
                # keep the reacting oxygen
                smalist = [sma.replace('[C:21](-[O:24])=[O:25]','[C:21](-[O:24])-[O:25]') for sma in smalist]
            
            else:
                # remove alcohol
                smalist = [sma.replace('[C:21](-[O:24])=[O:25]','[C:21]-[O:25]') for sma in smalist]
            
        elif self.modifArg == "+H2":
            if 24 in self.reactAtom:
                # keep the reacting oxygen
                smalist = [sma.replace('[C:21](-[O:24])=[O:25]','[C:21]-[O:24]') for sma in smalist]
            else:
                # fully reduce all oxygens
                smalist = [sma.replace('[C:21](-[O:24])=[O:25]','[C:21]') for sma in smalist]
                
        elif self.modifArg == None and 21 in self.reactAtom:
            # remove the alcohol if reaction atom is B[C] (to make space for bonding)
            smalist = [sma.replace('([C:21](-[O:24])=[O:25])','[C:21]=[O:25]') for sma in smalist]
            
        aro_smas = [sma.replace('[c:23][c:22]-[C:21]','[c:23][c:22][c:21]') for sma in smalist]
        aro_smas.extend([sma.replace('[c:21]-[O:25]','[c:21][o:25]').replace('[c:21]-[O:24]','[c:21][o:24]')
                         for sma in aro_smas])
        smalist.extend(aro_smas) 
        smalist = list(set(smalist))
        return smalist
    
    def Modify(self):
        self.smas = self.ModifySmarts(self.smas)
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
        return self
    
class Boronate:
    def __init__(self, h, X, Y = None):
        '''
        Namespace for Boronate, construct SMARTS and rdkit.Mol objects
        args:
            h - hybridization of the alpha and beta carbon
            X - reacting atom 
            Y - Modification (optional)
            
        The boronate molecule created by this class is used for 
        enumerating and constructing the C-H functionalization reactions
        '''      
        # Assign hybridization num and construct mol
        if h in (3, 'sp3'):
            self.smas = ['[C:53]-[C:52]-[B:51](-[O:54])-[O:55]']
            self.hybridization = 'sp3'
            
        elif h in (2, 'sp2'):
            self.smas = ['[C:53]=[C:52]-[B:51](-[O:54])-[O:55]',
                         '[c:53][c:52]-[B:51](-[O:54])-[O:55]']
            self.hybridization = 'sp2'
            
        else: 
            raise ValueError('Invalid input "{}" for hybridization!'.format(h)+\
                             ' Value of h must be 2, 3, "sp2", or "sp3".')
        
        # Assign reacting atom
        if X in (1,2,3,4,5,
                 *itertools.permutations(range(2,6),2)):
                 #(1,2),(1,3),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5),
                 #(2,1),(3,1),(3,2),(4,2),(5,2),(4,3),(5,3),(5,4)):
            if type(X) is int:
                self.reactAtom = (X+50,)
            else:
                self.reactAtom = (X[0]+50, X[1]+50)
                
        else:
            reactAtomDict = {'B':51, 'alpha':52 , 'beta':53, 
                             'B[O1]':54, 'B[O2]':55,}
            try:
                if type(X) is str:
                    self.reactAtom = (reactAtomDict[X],)
                else:
                    self.reactAtom = (reactAtomDict[X[0]], 
                                      reactAtomDict[X[1]])
            except (KeyError, TypeError): 
                raise ValueError(
                    'Invalid input "{}" for reacting atom!'.format(X)+\
                    ' Value of X must be 1, 2, 3, 4, 5 or a two-number combination from 2-5,\
                    or using strings "B", "alpha", "beta", "B[O1]", "B[O2]".'
                )
        # Assign modification
        if Y not in ('-C', None):
            raise ValueError('Invalid input "{}" for modification!'.format(Y)+\
                             ' Only "-C" (removing alkyl/allyl/aryl group) is allowed for boronate.')
            
        if (52 in self.reactAtom or 53 in self.reactAtom) and Y == "-C":
            raise ValueError('Conflict: "alpha/beta" with "-C". '+\
                             'Cannot remove the alpha or beta carbon while using it as reaction atom!')
        else: self.modifArg = Y
            
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
        
    def _repr_svg_(self):
        labels = [str(i+1) for i in range(len(self.smas))]
        return Draw._MolsToGridSVG(self.mols,legends=labels)
               
    def GetEnumRxnProps(self):
        # Report assigned properties
        print(self.__class__.__name__+"\n"+\
              "Hybridization of C : \th = {}\n".format(self.hybridization)+\
              "Reaction Atom(s) : \tX = {}\n".format(self.reactAtom)+\
              'Modification :\t\tY = {}\n'.format(self.modifArg))
        
        
    def ModifySmarts(self, smalist):
        if self.modifArg:
            smalist = [sma.replace('[C:53]-[C:52]-','').replace('[C:53]=[C:52]-','')\
                       .replace('[c:53][c:52]-','')
                       for sma in smalist]
            
        if 51 in self.reactAtom:
            smalist = [sma.replace('[B:51](-[O:54])-[O:55]','[B+:51](-[O:54])-[O:55]')
                       for sma in smalist]
            #smalist = ['[B+:51](-[O:54])-[O:55]']
        
        return smalist
    
    def Modify(self):
        return self
    
class Alcohol:
    '''
    Namespace for Alcohol, constructs SMARTS and rdkit.Mol objects
    args:
        h - hybridization of the alpha and beta carbon
        X - reacting atom 
        Y - modification (optional)
    
    The alcohol molecule created by this class is used for 
    enumerating and constructing the C-H functionalization reactions
    '''
    def __init__(self, h, X, Y = None):
        # Assign hybridization num and construct SMARTS
        if h in (3, 'sp3'):
            self.smas = ['[C:33]-[C:32]-[O:31]']
            self.hybridization = 'sp3'
            
        elif h in (2, 'sp2'):
            self.smas = ['[C:33]=[C:32]-[O:31]', '[c:33][c:32]-[O:31]']
            self.hybridization = 'sp2'
            
        else: 
            raise ValueError('Invalid input "{}" for hybridization!'.format(h)+\
                             ' Value of h must be 2, 3, "sp2", or "sp3".')
        
        # Assign reacting atom
        if X in (1, 2, 3, (1,2), (1,3), (2,3), (2,1), (3,1), (3,2)):
            # Add 30 to each number to match atom mapping number 
            # and avoid conflicts with other class mols when adding bonds
            if type(X) is int:
                self.reactAtom = (X+30,)
            else:
                self.reactAtom = (X[0]+30, X[1]+30)
        else:
            # Check if they are string inputs
            reactAtomDict = {'OH':31, 'O':31, 'alpha':32, 'beta':33}
            try:
                if type(X) is str:
                    self.reactAtom = (reactAtomDict[X],)
                else:
                    self.reactAtom = (reactAtomDict[X[0]], 
                                      reactAtomDict[X[1]])
            except (KeyError, TypeError): 
                raise ValueError(
                    'Invalid input "{}" for reacting atom!\n'.format(X)+\
                    'Value of X must be numerical values: 1, 2, 3, (1, 2), (1, 3), (2, 3)\n'+\
                    'Or strings: "O", "alpha", "beta"'
                )
                    
        # Assign modification
        if Y not in ('-O', '-B', None):
            raise ValueError('Invalid input "{}" for modification!'.format(Y)+\
                             ' Only "-O" or "-B" is allowed for alcohol.')
            
        # Below is not implemented, redundancies with C-C coupling are allowed for right now.
        #elif self.reactAtom in (22, 23, (22,23), (23,22)):
            # Removed OH to avoid redundancies with C-C coupling in enumerating rxns 
            #self.modifArg = '-O'   
        
        else: self.modifArg = Y
        
        # Checking args conflict
        if self.modifArg and 31 in self.reactAtom:
            raise ValueError('Conflict: X = "O" and Y = "-O". Cannot remove "O" while using it as reaction atom!')
            
        # Construct the molecule
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
        
    def _repr_svg_(self):
        labels = [str(i+1) for i in range(len(self.smas))]
        return Draw._MolsToGridSVG(self.mols,legends=labels)
    
    def GetEnumRxnProps(self):
        # Report assigned properties
        print(self.__class__.__name__+'\n'+\
              "Hybridization of C : \th = {}\n".format(self.hybridization)+\
              "Reaction Atom(s) : \tX = {}\n".format(self.reactAtom)+\
              'Modification : \t\tY = {}\n'.format(self.modifArg)
              )
    
    def ModifySmarts(self, smalist):
        # This is needed for separate treatments on individual SMARTS later
        smalist = list(smalist)
        if self.modifArg in ('-O', '-B'):
            smalist = [sma.replace('-[O:31]','') for sma in smalist]
    
        arolist = [sma.replace('[C:32]-[O:31]','[C:32]-[o:31]')\
                   .replace('[c:32]-[O:31]','[c:32][o:31]') 
                   for sma in smalist]
        smalist.extend(arolist)
        smalist = list(set(smalist))
        return smalist
        
    def Modify(self):
        self.smas = self.ModifySmarts(self.smas)
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
        return self
    
class Bromide:
    '''
    Namespace for Bromide, constructs SMARTS and rdkit.Mol objects
    args:
        h - hybridization of the alpha and beta carbon
        X - reacting atom 
        Y - modification (optional)
    
    The bromide molecule created by this class is used for 
    enumerating and constructing the C-H functionalization reactions
    '''
    def __init__(self, h, X, Y = None):
        # Assign hybridization num and construct SMARTS
        if h in (3, 'sp3'):
            self.smas = ['[C:43]-[C:42]-[Br:41]']
            self.hybridization = 'sp3'
            
        elif h in (2, 'sp2'):
            self.smas = ['[C:43]=[C:42]-[Br:41]', '[c:43][c:42]-[Br:41]']
            self.hybridization = 'sp2'
            
        else: 
            raise ValueError('Invalid input "{}" for hybridization!'.format(h)+\
                             ' Value of h must be 2, 3, "sp2", or "sp3".')
        
        # Assign reacting atom
        
        if X in (2, 3, (2,3), (3,2)):
            # Add 40 to each number to match atom mapping number 
            # and avoid conflicts with other class mols when adding bonds
            if type(X) is int:
                self.reactAtom = (X+40,)
            else:
                self.reactAtom = (X[0]+40, X[1]+40)
        else:
            # Check if they are string inputs
            reactAtomDict = {'alpha':42, 'beta':43}
            try:
                if type(X) is str:
                    self.reactAtom = (reactAtomDict[X],)
                else:
                    self.reactAtom = (reactAtomDict[X[0]], 
                                      reactAtomDict[X[1]])
            except (KeyError, TypeError): 
                raise ValueError(
                    'Invalid input "{}" for reacting atom!\n'.format(X)+\
                    'Value of X must be numerical values: 2, 3, (2, 3), (3, 2)\n'+\
                    'Or strings: "alpha", "beta" (Bromine itself has to be removed during rxn)'
                )
        
        # Assign modification
        # Y is actually an optional arg here 
        # as Bromine has to be removed even if it's not specified
        if Y not in ('-Br', '-B', None):
            raise ValueError('Invalid input "{}" for modification!'.format(Y)+\
                             ' Only "-Br" or "-B" is allowed for bromide.')
            
        else: self.modifArg = Y
        
        # Construct the molecule
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
        
    def _repr_svg_(self):
        labels = [str(i+1) for i in range(len(self.smas))]
        return Draw._MolsToGridSVG(self.mols,legends=labels)
    
    def GetEnumRxnProps(self):
        # Report assigned properties
        print(self.__class__.__name__+'\n'+\
              "Hybridization of C : \th = {}\n".format(self.hybridization)+\
              "Reaction Atom(s) : \tX = {}\n".format(self.reactAtom)+\
              'Modification : \t\tY = {}\n'.format(self.modifArg)
              )
    
    def ModifySmarts(self, smalist):
        if 42 in self.reactAtom or self.modifArg:
            # Remove Br
            smalist = list(smalist)
            smalist = [sma.replace('-[Br:41]','') for sma in smalist]
        return smalist
    
    def Modify(self):
        self.smas = self.ModifySmarts(self.smas)
        self.mols = [Chem.MolFromSmarts(sma) for sma in self.smas]
        return self
     