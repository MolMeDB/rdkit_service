from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem import rdDepictor as Depict
from rdkit.Chem import AllChem
import services.mmpa.rfrag as rfrag

class RDKIT:
    # Fragment molecule
    def mmpaFragment(self, params = {}):
        if "mol" not in params.keys():
            return None

        smiles = params["mol"]
        id = params["id"] if "id" in params.keys() else 1

        # First, canonize smiles
        canonized_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))

        if not canonized_smiles:
            return None

        fragments = rfrag.fragment_mol(canonized_smiles, id)

        result = []

        for i in fragments:
            t = i.split(",")
            if(len(t) < 4):
                continue
            f = t[3].split(".")
            f_res = []
            for ft in f:
                f_res.append(ft)
            result.append({
                "input": t[0],
                "identifier": t[1],
                "core": t[2],
                "fragments": f_res,
            })
        return result


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

        m3 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(m3)

        structure = Chem.MolToMolBlock(m3)

        return [structure]

    # Makes 2D structure
    def make2Dstructure(self, params = {}):
        if not "smi" in params:
            print("Parameter 'smi' not found in argument list. Please, fill some SMILES.")
            return []
        
        req_smiles = params["smi"]
        mol = Chem.MolFromSmiles(req_smiles)

        Depict.Compute2DCoords(mol)

        if not mol:
            raise Exception("Wrong SMILES INPUT")

        structure = Chem.MolToMolBlock(mol)

        return [structure]

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