from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem import rdDepictor as Depict
from rdkit.Chem import AllChem
import services.mmpa.rfrag as rfrag
import json
import re
import services.gen_confomers as cnf

class RDKIT:
    # Fragment molecule
    def mmpaFragment(self, params = {}, parent = None):
        if "mol" not in params.keys():
            return None

        smiles = params["mol"]
        id = params["id"] if "id" in params.keys() else 1
        print_result = True if "asHTML" in params.keys() else False

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
                f_res.append({
                    "smiles": ft,
                    "similarity": self.similarity(canonized_smiles, ft)
                })
            
            if(t[2]):
                f_res.append({
                    "smiles": t[2],
                    "similarity": self.similarity(canonized_smiles, t[2])
                })

            result.append({
                "input": t[0],
                "identifier": t[1],
                "fragments": f_res,
            })

        return True, result, print_result

    # Returns mol fingerprints
    def getFingerprint(self, params = {}):
        if "mol" not in params:
            return False

        s = params["mol"]

        mol = Chem.MolFromSmiles(s)
        fp = Chem.RDKFingerprint(mol)

        return {
            "smiles": s,
            "fingeprint": fp.ToBitString()
        }

    def similarity(self, smiles1, smiles2):
        if not smiles1 or not smiles2:
            return None

        smiles1 = re.sub(r'\[\*:[0-9]*\]', "[H]", smiles1)
        smiles2 = re.sub(r'\[\*:[0-9]*\]', "[H]", smiles2)

        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)

        fp1 = Chem.RDKFingerprint(mol1)
        fp2 = Chem.RDKFingerprint(mol2)

        return {
            "Tanimoto": DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.TanimotoSimilarity),
            "Dice": DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.DiceSimilarity),
            "Cosine": DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.CosineSimilarity),
            "Sokal": DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.SokalSimilarity),
            "Russel": DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.RusselSimilarity),
            "Kulczynski": DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.KulczynskiSimilarity),
            "McConnaughey": DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.McConnaugheySimilarity),
        }


    # Compute tanimoto betwween molecules
    def computeSimilarity(self, params = {}):
        if "smi1" not in params.keys() or "smi2" not in params.keys():
            return None

        smiles1 = params["smi1"]
        smiles2 = params["smi2"]

        return True, {
            "smiles1": smiles1,
            "smiles2": smiles2,
            "similarity": self.similarity(smiles1, smiles2)
        }


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

    # Returns all charge states for given molecule smiles
    def getAllChargeSmiles(self, params = {}):
        if params is None or "smi" not in params or params["smi"] is None:
            print("Missing SMI parameter")
            return []
        
        if "limit" not in params:
            limit = 20
        else:
            limit = int(params["limit"])

        smiles = params["smi"]
        # Check smiles validity
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            print("Invalid SMILES:", smiles)
            return []

        # Add library
        from dimorphite_dl import DimorphiteDL

        phLimits = [
            (0,14),
            (1,13),
            (2,11),
            (3,10),
            (4,9),
            (5,8),
            (6,8),
            (6.8,7.5),
            (7,7.5)
        ]

        for start, end in phLimits:
            dimorphite_dl = DimorphiteDL(
                min_ph=start, # Whole pH range
                max_ph=end,
                max_variants=128,
                label_states=False,
                pka_precision=1.0
            )

            # Get all protonated/deprotonated states
            prot_smiles_list = dimorphite_dl.protonate(smiles)

            if len(prot_smiles_list) < limit:
                break

        return {"pH_start": start, "pH_end": end, "molecules": prot_smiles_list}

        
    # Returns COSMO conformers
    def COSMO_conformers(self, params = {}):
        if params is None or "smi" not in params or params["smi"] is None:
            print("Missing `smi` parameter.")
            return []

        if "name" not in params or not params["name"]:
            print("Invalid structure name.")
            return []

        smiles = params["smi"]
        name = params["name"]
        # Check smiles validity
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            print("Invalid SMILES:", smiles)
            return []

        instance = cnf.Conformers()
        return instance.run(mol, name)

        

        