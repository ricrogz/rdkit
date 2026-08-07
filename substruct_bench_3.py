from rdkit import Chem

# this is from https://github.com/rdkit/rdkit/issues/8719

# Target molecule (a large, symmetric molecule)
mol_smi = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCCCCC(CCCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)(CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)COCCCCOCCOc1ccc(C[O:30][C@H:29]2[C@@H:20]([O:21][CH2:22][c:23]3[cH:24][cH:25][cH:26][cH:27][cH:28]3)[C@@H:10]([CH2:11][O:12][CH2:13][c:14]3[cH:15][cH:16][cH:17][cH:18][cH:19]3)[O:9][C@@H:8]([O:7][c:6]3[cH:5][cH:4][c:3]([O:2][CH3:1])[cH:41][cH:40]3)[C@@H:31]2[O:32][CH2:33][c:34]2[cH:35][cH:36][cH:37][cH:38][cH:39]2)cc1"
mol = Chem.MolFromSmiles(mol_smi)

# Query molecule (also contains symmetric segments, like the perfluorinated chains)
query_smi = "C(OC[C@H]1O[C@@H](Oc2ccc(OC)cc2)[C@H](OCc2ccccc2)[C@@H](OCc2ccc(OCCOCCCCOCC(CCCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)(CCCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)F)cc2)[C@H]1OCc1ccccc1)c1ccccc1"
query_mol = Chem.MolFromSmiles(query_smi)

# This call hangs / takes an extremely long time
for match in mol.GetSubstructMatches(query_mol):
    print(match)
