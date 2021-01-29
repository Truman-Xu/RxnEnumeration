from rdkit import Chem
from rdkit.Chem import Draw, AllChem

# Suppress waring messeges
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

def enumProds(reacts, rxns, drop_duplicates = False, report = True):
    products_lists = [rxn.RunReactants(reacts) for rxn in rxns]

    all_smi = []
    total = 0
    for l in products_lists:
        smi_list = []
        for p in l:
            try:
                Chem.SanitizeMol(p[0]) 
                smi = Chem.MolToSmiles(p[0])
                mol = Chem.MolFromSmiles(smi)
                smi_list.append(Chem.MolToSmiles(mol))
            except:
                pass
        all_smi.append(smi_list)
        total += len(smi_list)
    if drop_duplicates:
        all_smi = [list(set(l)) for l in all_smi]
    if report:
        print('Total products:',total)
    return all_smi

def drawProds(reacts, prod_smi_lists):
    
    AllChem.SetPreferCoordGen(True)
    AllChem.Compute2DCoords(reacts[1])
    AllChem.Compute2DCoords(reacts[0])

    pics = []
        
    for i, smi_list in enumerate(prod_smi_lists):
        if len(smi_list) > 50:
            smi_list = smi_list[:50]
            print("Truncated list %i to 50 mols"%i)
            
        elif len(smi_list) == 0:
            continue
        
        mols = [Chem.MolFromSmiles(m) for m in smi_list]
        matches = []

        for m in mols:
            try:
                match = list(m.GetSubstructMatches(reacts[1])[0])
                matches.append(match)
            except IndexError:
                pass
            
            AllChem.GenerateDepictionMatching2DStructure(m,reacts[0])

        labels = [str(s) for s in range(1, len(mols)+1)]
        pic = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(300, 300),
                                   legends = labels,
                                   #useSVG = True, 
                                   highlightAtomLists = matches)
        pics.append(pic)
        
    return pics
