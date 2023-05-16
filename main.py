from tkinter import *
from tkinter import filedialog
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from openbabel import openbabel as ob
import numpy

class Interface:
    def __init__(self, master):
        self.master = master
        master.title("Ma fenêtre principale")
        
        self.canvas_image = Canvas(master, width=300, height=200)
        self.canvas_image.pack()

        self.canvas_resultat = Canvas(master, width=300, height=300)
        self.canvas_resultat.pack()

        self.textvar_energie = StringVar()
        self.textvar_energie.set("Energie")

        self.label_energie = Label(master, textvariable=self.textvar_energie)
        self.label_energie.pack()
        
        self.button_ajout_molecule = Button(master, text="Ajouter une molecule", command=lambda: load_structure(self))
        self.button_ajout_molecule.pack()
    
    def show_molecule(self, mol):
        Draw.MolToFile(mol, "mol.png") # afficher la molécule dans un fichier image
        self.image_molecule = PhotoImage(file="mol.png")
        self.canvas_image.create_image(10,10,anchor=NW,image=self.image_molecule)
        self.textvar_energie.set("{}".format(getEigs_second_voisin(mol)[0]))

        self.master.update()
        return
    
    def show_results(self, mol):
        self.master.update()
        return

# fonction pour lire un fichier XYZ et afficher la structure moléculaire
def load_structure_(interface):
    filename = filedialog.askopenfilename(filetypes=[("Fichiers XYZ", "*.xyz")]) # ouvrir une boîte de dialogue pour choisir un fichier
    if filename:
        with open(filename, 'r') as f:
            lines = f.readlines()
        num_atoms = int(lines[0]) # nombre d'atomes dans la molécule
        atoms = []
        coords = []
        for i in range(2, num_atoms+2):
            tokens = lines[i].split()
            atoms.append(tokens[0])
            coords.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
        str = "{}\n\n".format(num_atoms)
        str = str + '\n'.join([f"{atoms[i]} {coords[i][0]} {coords[i][1]} {coords[i][2]}" for i in range(num_atoms)])
        mol = Chem.MolFromXYZBlock(str) # créer une molécule à partir des coordonnées
        mol.AddBond(mol.GetAtomWithIdx(0).GetIdx(), mol.GetAtomWithIdx(1).GetIdx(), Chem.rdchem.BondType.SINGLE)
        if not(mol):
            print("Erreur lors de la lecture du fichier")
    interface.show_molecule(mol)
    return

def load_structure(interface):
    filename = filedialog.askopenfilename(filetypes=[("Fichiers XYZ", "*.xyz")]) # ouvrir une boîte de dialogue pour choisir un fichier
    # créer un objet OBMol vide
    mol_ob = ob.OBMol()

    # lire le fichier xyz
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "mol")
    obConversion.ReadFile(mol_ob, filename)
    mol_rdkit = Chem.MolFromMolBlock(obConversion.WriteString(mol_ob))
    interface.show_molecule(mol_rdkit)
    return

def getEigs_premier_voisin(mol):
    mat = Chem.rdmolops.GetAdjacencyMatrix(mol).astype(numpy.float32)
    print(mat)
    eigvals, eigvecs = numpy.linalg.eig(mat)
    idx = numpy.argsort(-eigvals) #signe moins pour ordonner par ordre decroissant (beta<0)
    eigvals_sorted = eigvals[idx]
    eigvecs_sorted = eigvecs[:, idx]
    dim = len(eigvals_sorted)
    energy = 2 * sum(e for e in eigvals_sorted[:int(dim/2)])
    return energy, eigvals_sorted, eigvecs_sorted

def getVoisins(i, mol):
    mat = Chem.rdmolops.GetAdjacencyMatrix(mol).astype(numpy.float32)

    res=[]
    for j in range(len(mat[0])):
        if mat[i,j]==1:
            res.append(j)
    return res

def getVoisins_2(i, mol):
#    mat = Chem.rdmolops.GetAdjacencyMatrix(mol).astype(numpy.float32)

    res=[]
    first_Voisins = getVoisins(i, mol)
    for j in first_Voisins:
        res = res + getVoisins(j, mol)
    return res

def getEigs_second_voisin(mol):
    mat = Chem.rdmolops.GetAdjacencyMatrix(mol).astype(numpy.float32)

    for i in range(len(mat)):
        for j in getVoisins_2(i, mol):
            mat[i][j] = .2
#    print(mat)

    eigvals, eigvecs = numpy.linalg.eig(mat)
    idx = numpy.argsort(-eigvals) #signe moins pour ordonner par ordre decroissant (beta<0)
    eigvals_sorted = eigvals[idx]
    eigvecs_sorted = eigvecs[:, idx]
    dim = len(eigvals_sorted)
    energy = 2 * sum(e for e in eigvals_sorted[:int(dim/2)])
    return energy, eigvals_sorted, eigvecs_sorted

def main():
    root = Tk()
    app = Interface(root)
    root.mainloop()

if __name__ == "__main__":
    main()
