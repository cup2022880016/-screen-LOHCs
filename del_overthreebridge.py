from rdkit import Chem
from conring import share_ring_atoms
import pandas as pd



'''This code is to eliminate molecules containing atoms bridging more than three rings.'''

def del_overthreebridge(df,atom_element=None):
    '''Read the excel file obtained from the previous code. Only the code containing carbon is shown here.
    If it contains heteroatoms, please change the number after range and read the file name of Excel.'''
    df_C = df
    df_C = df_C.reset_index(drop=True)
    under3curatom_list = []
    df1 = pd.DataFrame()


    '''Determine whether the number of shared atoms in the ring does not exceed 2 by entering the smiles formula'''
    for index,row in df_C.iterrows():
        smi = row['IsomericSMILES']
        curatom_dir = share_ring_atoms(smi)

        if len(curatom_dir['curR1R2']) <= 2 and len(curatom_dir['curR1R3']) <= 2 and len(curatom_dir['curR2R3']) <= 2:
            under3curatom_list.append({'IsomericSMILES': row['IsomericSMILES']})

    new_df = pd.DataFrame(under3curatom_list)
    df1 = pd.concat([df1, new_df], ignore_index=True)

    return df1