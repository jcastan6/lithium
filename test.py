import psycopg2
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolFromSmiles

search_conn = psycopg2.connect('postgresql:///molecules', connect_timeout=1)
search_cursor = search_conn.cursor()
mol = MolFromSmiles("CC2CCC(C1CCCC1)C2")
mol= AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512)

search_cursor.execute("SELECT fingerprints.id, tanimoto_sml(ecfp4, %s) as similarity, mol_to_smiles(substance.mol) as smiles FROM fingerprints \
                        LEFT JOIN substance ON fingerprints.id = substance.id \
                        WHERE fingerprints.id >= %s AND fingerprints.id < %s ORDER BY tanimoto_sml(ecfp4, %s) DESC LIMIT 1", (str(mol), 0*10000000, (0+1)*10000000, str(mol)))
res = search_cursor.fetchall()
print(res)
search_cursor.close()
search_conn.close()

