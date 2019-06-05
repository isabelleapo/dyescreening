import rdkit, rdkit.Chem as rdkit
import stk
import subprocess as sp
import itertools
from .utilities import *
from joblib import Parallel, delayed
import re
import os

class DyeMaker:

    '''
    Builds dye structures and calculates IP, EA, optical gap and oscillator strength of first
    excitation energy using xTB.

    Parameters
    ----------
    path_to_A_smiles:  :class:'str'
        Path to text file containing list of smiles for monomer A

    path_to_B_smiles:  :class:'str'
        Path to text file containing list of smiles for monomer B

    path_to_C_smiles:  :class:'str'
        Path to text file containing list of smiles for monomer C

    nconfs:  :class:'int'
        Number of conformers to embed within conformer search

    solvent:  :class:'str'
        Solvent to be used in xTB (available solvents: toluene,thf, methanol, h2o, ether,
        chcl3, acetonitrile, acetone, cs2)

    min_osc_strength:  :class:'float'
        Minimum oscillator strength of the excitation selected to be the optical gap


    Methods
    -------
    dye_screen:
        Builds dyes made up of all possible combinations of monomers from libraries for A, B and C
        in the sequence ABCBA, performs a conformer search and calculates properties of the lowest
        energy conformers of these dyes.


    Returns
    -------
    .mol file containing lowest energy conformer

    Text file containing calculated properties

    '''

    def __init__(self,
                path_to_A_smiles,
                path_to_B_smiles,
                path_to_C_smiles,
                nconfs = 30,
                solvent = None,
                min_osc_strength = 0.0
                ):

        self.path_to_A_smiles = path_to_A_smiles
        self.path_to_B_smiles = path_to_B_smiles
        self.path_to_C_smiles = path_to_C_smiles
        self.nconfs = nconfs
        self.solvent = solvent
        self.min_osc_strength = min_osc_strength



    def _create_smiles_list(self):
        with open(self.path_to_A_smiles, 'r') as f:
            f2 = f.read()
            A_list = f2.split()

        with open(self.path_to_B_smiles, 'r') as f3:
            f4 = f3.read()
            B_list = f4.split()

        with open(self.path_to_C_smiles, 'r') as f5:
            f6 = f5.read()
            C_list = f6.split()

        return A_list, B_list, C_list

    def _combine_to_dyes(self, i, smile_tuple):
        A_smiles = smile_tuple[0]
        B_smiles = smile_tuple[1]
        C_smiles = smile_tuple[2]

        Amol = rdkit.AddHs(rdkit.MolFromSmiles(A_smiles))
        Bmol = rdkit.AddHs(rdkit.MolFromSmiles(B_smiles))
        Cmol = rdkit.AddHs(rdkit.MolFromSmiles(C_smiles))

        A = stk.StructUnit2.rdkit_init(Amol, 'bromine')
        B = stk.StructUnit2.rdkit_init(Bmol, 'bromine')
        C = stk.StructUnit2.rdkit_init(Cmol, 'bromine')

        try:
            dye1 = stk.Polymer([B, C, B], stk.Linear('ABA', [0,0,1], n=1, ends='fg'))
            X = stk.StructUnit2.rdkit_init(dye1.mol, 'bromine')
            dye = stk.Polymer([A, X, A], stk.Linear('ABA', [0,0,1], n=1))
            dyemol = dye.mol
            stk.rdkit_ETKDG(dye)
            id = '{:08d}'.format(i)
        except Exception as e:
            print(e, 'Build error')
        return A_smiles, B_smiles, C_smiles, dye, id

    def _conformer_search(self, dye, id):
        dye = dye.mol
        dye = rdkit.AddHs(dye)
        confs = rdkit.AllChem.EmbedMultipleConfs(dye, self.nconfs, rdkit.AllChem.ETKDG())
        rdkit.SanitizeMol(dye)

        lowest_energy = 10**10
        for conf in confs:
            ff = rdkit.AllChem.MMFFGetMoleculeForceField(dye, rdkit.AllChem.MMFFGetMoleculeProperties(dye), confId=conf)
            ff.Initialize()
            energy = ff.CalcEnergy()

            if energy < lowest_energy:
                lowest_energy = energy
                lowest_conf = conf

        rdkit.MolToMolFile(dye, f'lowest_conformers/{id}_lowest-conformer.mol', confId=lowest_conf)


    def _create_prop_files(self):
        os.mkdir('lowest_conformers')
        with open('props.txt', 'w') as f:
            f.write('Dye_ID\txTB_IP\txTB_EA\txTB_Opticalgap\txTB_oscillator_strength\tSmiles\tA\tB\tC\n')

    def _properties(self, A_smiles, B_smiles, C_smiles, id):
        with open(f'lowest_conformers/{id}_lowest-conformer.mol', 'r') as f:
            molfile = f.read()
        molfile = rdkit.MolFromMolBlock(molfile)
        molfile = rdkit.AddHs(molfile)
        rdkit.AllChem.EmbedMolecule(molfile, rdkit.AllChem.ETKDG())
        dirname = str(id)+'_properties'
        os.mkdir(dirname)
        os.chdir(dirname)
        xyz = generate_xyz(molfile)

        if self.solvent != None:
            p = sp.Popen(['xtb','xyz','-opt','-gbsa', self.solvent], stdout=sp.PIPE)
            output, _ = p.communicate()

        else:
            p = sp.Popen(['xtb','xyz','-opt'], stdout=sp.PIPE)
            output, _ = p.communicate()

        if self.solvent != None:
            p = sp.Popen(['xtb','xtbopt.xyz', '-vip', '-gbsa', self.solvent], stdout=sp.PIPE)
            output, _ = p.communicate()

        else:
            p = sp.Popen(['xtb','xtbopt.xyz', '-vip'], stdout=sp.PIPE)
            output, _ = p.communicate()

        with open('temp.txt', 'wb') as f:
            f.write(output)

        with open('temp.txt', 'rb') as f:
            temp = f.read()
            temp = str(temp)

        pattern = re.compile('(?<=delta SCC IP [(]eV[)].).*\d\.\d{4}')
        ip = pattern.findall(temp)

        rt = sp.Popen(['rm','temp.txt'], stdout=sp.PIPE)
        o, e = rt.communicate()

        if self.solvent != None:
            p = sp.Popen(['xtb','xtbopt.xyz', '-vea', '-gbsa', self.solvent], stdout=sp.PIPE)
            output, _ = p.communicate()

        else:
            p = sp.Popen(['xtb','xtbopt.xyz', '-vea'], stdout=sp.PIPE)
            output, _ = p.communicate()

        with open('temp.txt', 'wb') as f:
            f.write(output)

        with open('temp.txt', 'rb') as f:
            temp = f.read()
            temp = str(temp)

        pattern = re.compile('(?<=delta SCC EA [(]eV[)].).*\d\.\d{4}')
        ea = pattern.findall(temp)

        rt = sp.Popen(['rm','temp.txt'], stdout=sp.PIPE)
        o, e = rt.communicate()

        p = sp.Popen(['xtb','xtbopt.xyz'], stdout=sp.PIPE)
        output, _ = p.communicate()

        p = sp.Popen(['stda','-xtb', '-e', '8'], stdout=sp.PIPE)
        output, _ = p.communicate()
        with open('temp.txt', 'wb') as f:
            f.write(output)

        with open('temp.txt', 'rb') as f:
            temp = f.read()
            temp = str(temp)

        pattern = re.compile(r'\\n\s+\d+\s+(\d\.\d+)\s+[\d.-]+\s+(\d+\.\d+)[>\(\)\s\d.-]*')
        og = pattern.findall(temp)

        rt = sp.Popen(['rm','temp.txt'], stdout=sp.PIPE)
        o, e = rt.communicate()

        smilenoH = rdkit.RemoveHs(molfile)
        smiles = rdkit.MolToSmiles(smilenoH, canonical=True)

        string = str(id)+'\t'
        for match in ip:
            match = match.strip()
            string += match+'\t'
        for match in ea:
            match = match.strip()
            string += match+'\t'
        opgap = []
        osc = []
        for match in og:
            if float(match[1]) > self.min_osc_strength:
                opgap.append(match[0])
                osc.append(match[1])
                break
        string += opgap[0]+'\t'
        string += osc[0]+'\t'
        string += smiles+'\t'
        string += A_smiles+'\t'
        string += B_smiles+'\t'
        string += C_smiles+'\n'

        with open('props.txt', 'w') as f:
            f.write(string)
        os.chdir('../')


    def _write_to_props(self):
       dirs = str(os.listdir('.'))
       pattern = re.compile('\d+_properties')
       match = pattern.findall(dirs)
       for m in match:
           try:
               with open(m+'/props.txt', 'r') as p:
                   prop = p.read()
               with open('props.txt', 'a') as f:
                   f.write(prop)
           except Exception as e:
               print(e, 'Error in writing to file')

    def _dye_screen(self, i, item):
        try:
            A_smiles, B_smiles, C_smiles, dye, id = self._combine_to_dyes(i, item)
            self._conformer_search(dye, id)
        except Exception as e:
            print(e, 'Error in dye formation or conformer search')
        try:
            self._properties(A_smiles, B_smiles, C_smiles, id)
        except Exception as e:
            print(e, 'Error in property calculation')
            os.chdir('../')

    def screen(self, nprocs=1):
        self._create_prop_files()
        A_list, B_list, C_list = self._create_smiles_list()
        results = Parallel(n_jobs=nprocs)(delayed(self._dye_screen)(i, item) for i, item in enumerate(itertools.product(A_list, B_list, C_list, repeat=1)))
        self._write_to_props()
