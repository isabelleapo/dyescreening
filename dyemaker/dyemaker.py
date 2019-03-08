import rdkit, rdkit.Chem as rdkit
import stk
import subprocess as sp
import itertools
from .utilities import *

class DyeMaker:

    def __init__(self,
                path_to_A_smiles,
                path_to_B_smiles,
                path_to_C_smiles,
                nconfs = 30,
                #solvent = None
                ):

        self.path_to_A_smiles = path_to_A_smiles
        self.path_to_B_smiles = path_to_B_smiles
        self.path_to_C_smiles = path_to_C_smiles
        self.nconfs = nconfs
        #self.solvent = solvent

    def CreateSmilesList(self):
        with open(self.path_to_A_smiles + '.txt', 'r') as f:
            f2 = f.read()
            A_list = f2.split()

        with open(self.path_to_B_smiles + '.txt', 'r') as f3:
            f4 = f3.read()
            B_list = f4.split()

        with open(self.path_to_C_smiles + '.txt', 'r') as f5:
            f6 = f5.read()
            C_list = f6.split()

        return A_list, B_list, C_list

    # now we have a three lists containing all the strings


    def DyeSmiletoMol(self, smiles):
        molobject = rdkit.MolFromSmiles(smiles)
        rdkit.AddHs(molobject)
        rdkit.AllChem.EmbedMolecule(molobject, rdkit.AllChem.ETKDG())
        return molobject
# this turns smiles into mol object and we will feed in the lists at the end to this function (for i in Alist: DyeSmiletoMol(i))

    def MakeDyeUnitA(self, molobject):
        A = stk.StructUnit2.rdkit_init(molobject, 'bromine')
        return A

    def MakeDyeUnit(self, molobject):
        B = stk.StructUnit2.rdkit_init(molobject, 'bromine')
        return B

#workflow


#make 3 lists of smiles

#now make lists of mol objects
    def MolLists(self, A_list, B_list, C_list):
        Amols = []
        Bmols = []
        Cmols = []
        for item in A_list:
            A = self.DyeSmiletoMol(item)
            Amols.append(A)

        for item in B_list:
            B = self.DyeSmiletoMol(item)
            Bmols.append(B)

        for item in C_list:
            C = self.DyeSmiletoMol(item)
            Cmols.append(C)

        return Amols, Bmols, Cmols

    def CombinerDyes(self, Amols, BMols, Cmols):
        dye_combo = itertools.product(Amols, BMols, Cmols, repeat=1)
        dye_combo = list(dye_combo)
        dyes = []
        ids = []
        Amol = [item[0] for item in dye_combo]
        Bmol = [item[1] for item in dye_combo]
        Cmol = [item[2] for item in dye_combo]
        for i, value in enumerate(dye_combo):
            A = self.MakeDyeUnitA(Amol[i])
            B = self.MakeDyeUnit(Bmol[i])
            C = self.MakeDyeUnit(Cmol[i])

            dye = stk.Polymer([A, B, C], stk.Linear('ABCBA', [0,0,0,1,1], n=1))
            stk.rdkit_ETKDG(dye)
            dyes.append(dye)

            ids.append('{:04d}'.format(i))
        return dyes, ids


    def conformer_search(self, dyes, ids):
        ndyes = []
        for i, dye in enumerate(dyes):
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
                ndyes.append(dye)
            print(ndyes)
            print(ids)
            return ndyes, ids
            rdkit.MolToMolFile(dye, f'{ids[i]}_lowest-conformer.mol', confId=lowest_conf)

#IN PROGRESS!!!!!!!!!!-----------------------------------------------------------------------------------------------------
    def _properties(self, molfile, id):

        dirname = str(id)+'_properties'
        p = sp.Popen(['mkdir', dirname], stdin=sp.PIPE)
        p, o = p.communicate()
        cd(dirname)
        xyz = generate_xyz(molfile)
#---OPTIMISE

        p = sp.Popen(['xtb','xyz.xyz','-opt'], stdout=sp.PIPE)
        output, _ = p.communicate()

#---CALCULATE IP


        p = sp.Popen(['xtb','xtbopt.xyz', '-vip'], stdout=sp.PIPE)
        output, _ = p.communicate()
        with open('temp.txt', 'wb') as f:
            f.write(output)
            f.close()
            with open('temp.txt', 'rb') as f:
                temp = f.read()
                temp = str(temp)

        with open('f.txt', 'w') as f:
            f.write('Smiles     IP (eV)     EA (eV)     Optical gap (eV)     Oscillator strength\n')

        pattern = re.compile('(?<=delta SCC IP [(]eV[)].).*\d\.\d{4}')
        ip = pattern.findall(temp)
        with open('f.txt', 'a') as f:
            for match in ip:
                match = match.strip()
                f.write(match)

        rt = sp.Popen(['rm','temp.txt'], stdout=sp.PIPE)
        o, e = rt.communicate()

#---CALCULATE EA

        p = sp.Popen(['xtb','xtbopt.xyz', '-vea'], stdout=sp.PIPE)
        output, _ = p.communicate()
        with open('temp.txt', 'wb') as f:
            f.write(output)
            f.close()
            with open('temp.txt', 'rb') as f:
                temp = f.read()
                temp = str(temp)

        pattern = re.compile('(?<=delta SCC EA [(]eV[)].).*\d\.\d{4}')
        with open('f.txt', 'a') as f:
            f.write('        ')
            for match in ip:
                match = match.strip()
                f.write(match)

        rt = sp.Popen(['rm','temp.txt'], stdout=sp.PIPE)
        o, e = rt.communicate()

#---CALCULATE OPTICAL GAP

        p = sp.Popen(['xtb','xyzopt.xyz'], stdout=sp.PIPE)
        output, _ = p.communicate()

        p = sp.Popen(['stda','-xtb', '-e', '8'], stdout=sp.PIPE)
        output, _ = p.communicate()
        with open('temp.txt', 'wb') as f:
            f.write(output)
            f.close()
            with open('temp.txt', 'rb') as f:
                temp = f.read()
                temp = str(temp)

        pattern = re.compile('(?<=Rv[(]corr[)]\n    1).*\d\.\d{3}(?<=\d).*\d\.\d{3}')
        ip = pattern.findall(temp)
        with open('f.txt', 'a') as f:
            f.write('        ')
            for match in ip:
                match = match.strip()
                f.write(match)

        rt = sp.Popen(['rm','temp.txt'], stdout=sp.PIPE)
        o, e = rt.communicate()
        os.chdir('../')

#IN PROGRESS!!!!!!!!!!-----------------------------------------------------------------------------------------------------~

    def compute_properties(self, ndyes, ids):
        for dye in ndyes:
            self._properties(dye, ids)

    def DyeScreen(self):
        A_list, B_list, C_list = self.CreateSmilesList()
        Amols, Bmols, Cmols = self.MolLists(A_list, B_list, C_list)
        dyes, ids = self.CombinerDyes(Amols, Bmols, Cmols)
        ndyes = self.conformer_search(dyes, ids)
        self.compute_properties(ndyes, ids)
