from rdkit import Chem
import pandas as pd





def noH_atom(mol):
    '''Define a function to find a C atom that is not connected to hydrogen,
    that is, a carbon atom that cannot be bonded to another carbon atom.'''
    noH_atom_list = []

    atoms = mol.GetAtoms()
    for atom in atoms:
        if atom.IsInRing() == True:
            CnoH = atom.GetTotalNumHs()
            if CnoH == 0:
                noH_atom_list.append(atom.GetIdx())
    return noH_atom_list

def hasH_atom(mol):
    '''Define a function to find the C atom connected with hydrogen,
    that is, the carbon atom that can form a bond with another carbon atom.'''
    noH_atom_list = []

    atoms = mol.GetAtoms()
    for atom in atoms:
        if atom.IsInRing() == True:
            CnoH = atom.GetTotalNumHs()
            if CnoH != 0:
                noH_atom_list.append(atom.GetIdx())
    return noH_atom_list

def between_noH(noH_idx_list,mol):
    '''Define a function to find whether there are atoms between two atoms without hydrogen connection.'''
    nei_list = []

    for noH_idx in noH_idx_list:
        nei_atoms = mol.GetAtomWithIdx(noH_idx).GetNeighbors()
        for nei_atom in nei_atoms:
            if nei_atom.IsInRing() == True and nei_atom.GetIdx() not in noH_idx_list:
                nei_list.append(nei_atom.GetIdx())

    dup = list(set([x for x in nei_list if nei_list.count(x) > 1]))

    return dup

def atom_in_ring(mol):
    _list1 = []
    ssr = Chem.GetSymmSSSR(mol)
    for i in ssr:
        for a in i:
            if a not in _list1:
                _list1.append(a)

    return _list1












