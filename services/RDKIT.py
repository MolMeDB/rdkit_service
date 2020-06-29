from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem import rdDepictor as Depict

class RDKIT:
    # Canonize SMILES
    def canonizeSmiles(self, params = {}):
        if not "smi" in params:
            print("Parameter 'smi' not found in argument list")
            return []
        
        req_smiles = params["smi"]
        canonized_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(req_smiles))
        return [canonized_smiles]

    # Makes 3D structure
    def make3Dstructure(self, params = {}):
        if not "smi" in params:
            print("Parameter 'smi' not found in argument list. Please, fill some SMILES.")
            return []
        
        req_smiles = params["smi"]
        mol = Chem.MolFromSmiles(req_smiles)

        Depict.Compute2DCoords(mol)

        if not mol:
            raise Exception("Wrong SMILES INPUT")

        return [Chem.MolToMolBlock(mol)]

    # Generate InChIKey
    def makeInchi(self, params = {}):
        if not "smi" in params:
            print("Parameter 'smi' not found in argument list")
            return []
        
        req_smiles = params["smi"]
        inchi = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(req_smiles))
        return [inchi]

    # Compute LogP
    def getGeneralInfo(self, params = {}):
        if not "smi" in params:
            print("Parameter 'smi' not found in argument list")
            return []
        
        req_smiles = params["smi"]
        mol = Chem.MolFromSmiles(req_smiles)

        canonized_smiles = Chem.MolToSmiles(mol)
        inchi = Chem.inchi.MolToInchi(mol)
        inchiKey = Chem.inchi.InchiToInchiKey(inchi)
        MW = Descriptors.MolWt(mol)
        LogP = Crippen.MolLogP(mol)

        return {
            "canonized_smiles": canonized_smiles,
            "inchi": inchiKey,
            "MW": MW,
            "LogP": LogP
        }

    # Fingerprint similarity