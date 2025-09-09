from rdkit import Chem




def compose_ring_atoms(smiles):
    ring_list = []
    mol = Chem.MolFromSmiles(smiles)
    ssr = Chem.GetSymmSSSR(mol)
    for ring in ssr:
        ring_list.append(list(ring))


    def ring_comatoms(ring_list,ring_num1,ring_num2):
        del_list = []
        for num_index in ring_list[ring_num2]:
            if num_index in ring_list[ring_num1]:
                del_list.append(num_index)
        for num_index in del_list:
            ring_list[ring_num2].remove(num_index)
        return ring_list


    ring_num = len(ring_list)
    if ring_num == 2:
        ring_comatoms(ring_list,0,1)
    elif ring_num == 3:
        ring_comatoms(ring_list,0,1)
        ring_comatoms(ring_list,0,2)
        ring_comatoms(ring_list,1,2)

    return ring_list


def share_ring_atoms(smiles):
    '''Calculate the shared atoms between rings in the molecule (handle up to 3 rings)'''
    ring_list  = []
    dir = {}
    mol = Chem.MolFromSmiles(smiles)
    ssr = Chem.GetSymmSSSR(mol)

    def ring_curatoms(ring_list, ring_num1, ring_num2):
        '''Auxiliary function: returns the shared atomic index list of two rings'''
        current_list = []
        '''Check whether the atoms of ring_num2 exist in ring_num1'''
        for num_index in ring_list[ring_num2]:
            if num_index in ring_list[ring_num1]:
                current_list.append(num_index)
        return current_list

    for ring in ssr:
        ring_list.append(list(ring))
        ring_num = len(ring_list)
        if ring_num == 1:
            dir['curR1R2'] = []
            dir['curR1R3'] = []
            dir['curR2R3'] = []
        if ring_num == 2:
            dir['curR1R2'] = ring_curatoms(ring_list,0,1)
            dir['curR1R3'] = []
            dir['curR2R3'] = []
        elif ring_num == 3:
            dir['curR1R2'] = ring_curatoms(ring_list,0,1)
            dir['curR1R3'] = ring_curatoms(ring_list,0,2)
            dir['curR2R3'] = ring_curatoms(ring_list,1,2)

    return dir


def IsRingConnect(smiles):
    ring_list = []
    share_atoms_R1R2 = []
    share_atoms_R1R3 = []
    share_atoms_R2R3 = []
    dir = {}
    mol = Chem.MolFromSmiles(smiles)
    ssr = Chem.GetSymmSSSR(mol)
    for i in ssr:
        ring_list.append(list(i))

    ring_num = len(ring_list)
    if ring_num == 2:
        for i in ring_list[1]:
            if i in ring_list[0]:
                share_atoms_R1R2.append(i)
            if len(share_atoms_R1R2) != 0:
                dir['R1R2'] = True
                dir['R2R3'] = False
            else:
                dir['R1R2'] = False
                dir['R2R3'] = False
    if ring_num == 3:
        for i in ring_list[1]:
            if i in ring_list[0]:
                share_atoms_R1R2.append(i)
        if len(share_atoms_R1R2) != 0:
            dir['R1R2'] = True
        else:
            dir['R1R2'] = False

        for i in ring_list[2]:
            if i in ring_list[0]:
                share_atoms_R1R3.append(i)
        if len(share_atoms_R1R3) != 0:
            dir['R1R3'] = True
        else:
            dir['R1R3'] = False

        for i in ring_list[2]:
            if i in ring_list[1]:
                share_atoms_R2R3.append(i)
        if len(share_atoms_R2R3) != 0:
            dir['R2R3'] = True
        else:
            dir['R2R3'] = False

    return dir

