# Script for finding substructure match using kekulized smiles (ignore aromatic bond matching)
# Reproduces results as seen in https://doi.org/10.1038/s41586-020-2142-y

# BSD 3-Clause License
# Copyright (c) 2020-2021, Cernak Lab at the University of Michigan

from rxnConstructor import RxnConstructor, enumAllRxns
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd

# Suppress waring messeges
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def create_df():
    rxns, rxnConds = enumAllRxns(

                    # hybridization
                    hA = ('sp3', 'sp2'),
                    hB = ('sp3', 'sp2'),

                    # reaction atom
                    XA = ("NH2", "alpha", "beta"),
                    XB = ("B[O]", "B[C]", "alpha", "beta"),

                    # modification
                    YA = (None, '-A'),
                    YB = (None, '-B', '-OH', '+H', '+H2'),

                    # add an aditional beta carbon
                    addBetaA = (False,),
                    addBetaB = (False,),
    )


    df = pd.DataFrame(rxnConds,columns=['Hybridization A', 'Rxn Atom A', 
                                        'Modification A', 'Add Beta A',
                                        'Hybridization B', 'Rxn Atom B', 
                                        'Modification B', 'Add Beta B'])
    df['Reactions Smarts'] = [rxnConst.rxnSmarts for rxnConst in rxns]
    #df.to_csv('enumRxns320.csv', index = False)

    # Add mol and reaction objects
    df['Product'] = [rxnConst.prod for rxnConst in rxns]
    df['Reactions'] = [rxnConst.rxn for rxnConst in rxns]
    
    return df
    
def subMatchSMI(df, drug_smi, img_file = 'drugSubMatches.png'):
    drug_mol = Chem.MolFromSmiles(drug_smi)

    Chem.Kekulize(drug_mol, clearAromaticFlags=True)

    matches = []
    rxn_idx = []
    for i,prod in enumerate(df['Product']):
        match = drug_mol.GetSubstructMatches(prod)
        if len(match) > 0:
            matches.append(match[0])
            rxn_idx.append(i)

    print('Matched:', len(matches))
    
    labels = ['RXN '+str(i) for i in rxn_idx]
    
    mols = [drug_mol for m in matches]
    pic = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(200, 200),
                               legends = labels,
                               #useSVG = True, 
                               highlightAtomLists = matches)

    pic.save(img_file)

if __name__ == '__main__':
    
    noscapine = "[H][C@]1([C@]2([H])C(C=CC(OC)=C3OC)=C3C(O2)=O)C4=C(OC)C5=C(OCO5)C=C4CCN1C"
    quinine = "C=C[C@@H]1[C@@H]2CC([C@@H](C3=C(C=C(OC)C=C4)C4=NC=C3)O)[N@](C1)CC2"
    sitagliptin = "FC1=CC(C[C@H](N)CC(N2CCN3C(C2)=NN=C3C(F)(F)F)=O)=C(F)C(F)=C1"
    drugs = [(noscapine,'noscapine.png'),
             (quinine,'quinine.png'),
             (sitagliptin,'sitagliptin.png')]
    
    df = create_df()
    for drug in drugs:
        subMatchSMI(df,*drug)
    