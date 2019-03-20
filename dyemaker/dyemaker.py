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
                solvent = None
                ):

        self.path_to_A_smiles = path_to_A_smiles
        self.path_to_B_smiles = path_to_B_smiles
        self.path_to_C_smiles = path_to_C_smiles
        self.nconfs = nconfs
        self.solvent = solvent



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



    def _dye_smile_to_mol(self, smiles):
        molobject = rdkit.MolFromSmiles(smiles)
        rdkit.AddHs(molobject)
        rdkit.AllChem.EmbedMolecule(molobject, rdkit.AllChem.ETKDG())
        return molobject



    def _make_dye_unit_a(self, molobject):
        A = stk.StructUnit2.rdkit_init(molobject, 'bromine')
        return A



    def _make_dye_unit(self, molobject):
        B = stk.StructUnit2.rdkit_init(molobject, 'bromine')
        return B



    def _mol_lists(self, A_list, B_list, C_list):
        Amols = []
        Bmols = []
        Cmols = []
        for item in A_list:
            A = self._dye_smile_to_mol(item)
            Amols.append(A)

        for item in B_list:
            B = self._dye_smile_to_mol(item)
            Bmols.append(B)

        for item in C_list:
            C = self._dye_smile_to_mol(item)
            Cmols.append(C)

        return Amols, Bmols, Cmols



    def _combiner_dyes(self, Amols, BMols, Cmols):
        dye_combo = itertools.product(Amols, BMols, Cmols, repeat=1)
        dye_combo = list(dye_combo)
        dyes = []
        ids = []
        Amol = [item[0] for item in dye_combo]
        Bmol = [item[1] for item in dye_combo]
        Cmol = [item[2] for item in dye_combo]
        for i, value in enumerate(dye_combo):
            #A = self._make_dye_unit_a(value[0])
            #B = self._make_dye_unit(value[1])
            #C = self._make_dye_unit(value[2])

            A = self._make_dye_unit_a(Amol[i])
            B = self._make_dye_unit(Bmol[i])
            C = self._make_dye_unit(Cmol[i])

            try:
                dye = stk.Polymer([A, B, C], stk.Linear('ABCBA', [0,0,0,1,1], n=1))
                stk.rdkit_ETKDG(dye)
                dyes.append(dye)
                ids.append('{:08d}'.format(i))
            except Exception as e:
                print('Conformer Error')

        n_dyes = zip(dyes, ids)
        return n_dyes



    def _conformer_search(self, dye, id):
        dye = dye.mol
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
            f.write('Dye ID\tIP (eV)\tEA (eV)\tOptical gap (eV)\tOscillator strength\tSmiles\n')
 
    def _properties(self, molfile, id):
        molfile = molfile.mol
        dirname = str(id)+'_properties'
        os.mkdir(dirname)
        os.chdir(dirname)
        xyz = generate_xyz(molfile)
        smiles = rdkit.MolToSmiles(molfile)

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

        pattern = re.compile(r'Rv\(corr\)\\n\s\s\s\s1\s*([\d\.-]+)\s+[\d\.-]+\s+([\d\.-]+)')
        og = pattern.findall(temp)

        rt = sp.Popen(['rm','temp.txt'], stdout=sp.PIPE)
        o, e = rt.communicate()


        with open('props.txt', 'w') as f:
            string = str(id)+'\t'
            for match in ip:
                match = match.strip()
                string += match+'\t'
            for match in ea:
                match = match.strip()
                string += match+'\t'
            for match in og:
                string += match[0]+'\t'
                string += match[1]+'\t'
            string += smiles+'\n'
            f.write(string)
        os.chdir('../')

    def _write_to_props(self):
        dirs = str(os.listdir('.'))
        pattern = re.compile('\d+_properties')
        match = pattern.findall(dirs)
        for m in match:
            with open(m+'/props.txt', 'r') as p:
                prop = p.read()
            with open('props.txt', 'a') as f:
                f.write(prop)


    def _compute_properties(self, dye, id):
        try:
            self._conformer_search(dye, id)
            self._properties(dye, id)
        except Exception as e:
            print(e, 'ERROR')


    def dye_screen(self, nprocs=1):
        
        A_list, B_list, C_list = self._create_smiles_list()
        Amols, Bmols, Cmols = self._mol_lists(A_list, B_list, C_list)
        n_dyes = self._combiner_dyes(Amols, Bmols, Cmols)
        self._create_prop_files()

        results = Parallel(n_jobs=nprocs)(delayed(self._compute_properties)(dye, id) for dye, id in n_dyes)
        self._write_to_props()
