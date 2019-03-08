from contextlib import contextmanager
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess as sp

def generate_xyz(smiles, xtb):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    print(Chem.MolToMolBlock(mol), file=open('mol','w+'))

    p = sp.Popen(['babel', 'mol', 'xyz'], stdout=sp.PIPE)
    o, e = p.communicate()


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
