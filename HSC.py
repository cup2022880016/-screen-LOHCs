from rdkit import Chem
import pandas as pd
from query_H_atom import noH_atom,between_noH
from conring import share_ring_atoms,IsRingConnect
import SelCarbon as sc



'''This code is based on smiles formula to obtain the ratio of hydrogen obtained after complete 
dehydrogenation of the molecule to molecular weight, i.e. weight percent'''


def cal_Hweight(df,atom_element=None):
    if atom_element == 'O':
        df_O = df

        H2 = 2
        wt_list = []
        wt_df = pd.DataFrame(columns=['IsomericSMILES', 'wt', 'molecule', 'H2num'])
        for index, row in df_O.iterrows():
            threering_list = []
            ring_atom_list = []

            noH_atom_list = []
            final_H2 = 0
            smi = row['IsomericSMILES']
            mol = Chem.MolFromSmiles(smi)
            ssr = Chem.GetSymmSSSR(mol)
            total_atom = mol.GetAtoms()
            noH_atom_idx = noH_atom(mol)
            dups = between_noH(noH_atom_idx, mol)

            H_num = 2 * len(total_atom)
            one_cyc = row['IsomericSMILES'].find('1')
            count_1 = row['IsomericSMILES'].count('1')
            two_cyc = row['IsomericSMILES'].find('2')
            three_cyc = row['IsomericSMILES'].find('3')
            if one_cyc != -1 and two_cyc != -1 and three_cyc != -1:
                decrease_H = 2 * 3
            elif one_cyc != -1 and two_cyc != -1 and three_cyc == -1:
                decrease_H = 2 * 2
            elif one_cyc != -1 and two_cyc == -1 and three_cyc == -1:
                if count_1 == 4:
                    decrease_H = 2 * 2
                elif count_1 == 2:
                    decrease_H = 2 * 1
                elif count_1 == 6:
                    decrease_H = 2 * 3
            total_H = H_num - decrease_H
            molecule = f'C{len(total_atom) - 1}H{total_H}O'

            if len(ssr) == 1:
                for ring in ssr:
                    atom_num = len(list(ring)) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = (atom_num - remainder) / 2
                    print(final_H2)
            if len(ssr) == 2:
                ringcon_dir = IsRingConnect(smi)
                if ringcon_dir['R1R2'] == True:
                    shareatom_dir = share_ring_atoms(smi)
                    for ring in ssr:
                        for i in ring:
                            ring_atom_list.append(i)
                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx) - len(shareatom_dir['curR1R2'])
                    remainder = atom_num % H2
                    final_H2 = (atom_num - remainder) / 2
                    print(final_H2)
                else:
                    for ring in ssr:
                        ring_list = list(ring)
                        for i in noH_atom_idx:
                            if i in ring_list:
                                ring_list.remove(i)
                        for j in dups:
                            if j in ring_list:
                                ring_list.remove(j)
                        atom_num = len(ring_list)
                        remainder = atom_num % H2
                        cur_H2 = (atom_num - remainder) / 2
                        final_H2 += cur_H2
                    print(final_H2)
            if len(ssr) == 3:
                ringcon_dir = IsRingConnect(smi)
                if ringcon_dir['R1R2'] == False and ringcon_dir['R2R3'] == False and ringcon_dir['R1R3'] == False:
                    for ring in ssr:
                        ring_list = list(ring)
                        for i in noH_atom_idx:
                            if i in ring_list:
                                ring_list.remove(i)
                        for j in dups:
                            if j in ring_list:
                                ring_list.remove(j)
                        atom_num = len(ring_list)
                        remainder = atom_num % H2
                        cur_H2 = (atom_num - remainder) / 2
                        final_H2 += cur_H2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == True and ringcon_dir['R2R3'] == False and ringcon_dir['R1R3'] == False:
                    shareatom_dir = share_ring_atoms(smi)
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for i in range(0, 2):
                        for j in threering_list[i]:
                            ring_atom_list.append(j)
                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx) - len(shareatom_dir['curR1R2'])
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2

                    for i in noH_atom_idx:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)
                    for j in dups:
                        if j in threering_list[2]:
                            threering_list[2].remove(j)
                    atom_num = len(threering_list[2])
                    remainder = atom_num % H2
                    cur_H2 = (atom_num - remainder) / 2
                    final_H2 += cur_H2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == True and ringcon_dir['R2R3'] == True and ringcon_dir['R1R3'] == True:
                    shareatom_dir = share_ring_atoms(smi)
                    ring3_shareatom = []
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for atomR1R2 in shareatom_dir['curR1R2']:
                        if atomR1R2 in threering_list[1]:
                            threering_list[1].remove(atomR1R2)
                    for atomR1R3 in shareatom_dir['curR1R3']:
                        if atomR1R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR1R3)
                    for atomR2R3 in shareatom_dir['curR2R3']:
                        if atomR2R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR2R3)
                    for i in ring3_shareatom:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)

                    for j in threering_list:
                        for n in j:
                            ring_atom_list.append(n)

                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == True and ringcon_dir['R2R3'] == True and ringcon_dir['R1R3'] == False:
                    shareatom_dir = share_ring_atoms(smi)
                    ring3_shareatom = []
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for atomR1R2 in shareatom_dir['curR1R2']:
                        if atomR1R2 in threering_list[1]:
                            threering_list[1].remove(atomR1R2)
                    for atomR1R3 in shareatom_dir['curR1R3']:
                        if atomR1R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR1R3)
                    for atomR2R3 in shareatom_dir['curR2R3']:
                        if atomR2R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR2R3)
                    for i in ring3_shareatom:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)

                    for j in threering_list:
                        for n in j:
                            ring_atom_list.append(n)

                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == False and ringcon_dir['R2R3'] == True and ringcon_dir['R1R3'] == True:
                    shareatom_dir = share_ring_atoms(smi)
                    ring3_shareatom = []
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for atomR1R2 in shareatom_dir['curR1R2']:
                        if atomR1R2 in threering_list[1]:
                            threering_list[1].remove(atomR1R2)
                    for atomR1R3 in shareatom_dir['curR1R3']:
                        if atomR1R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR1R3)
                    for atomR2R3 in shareatom_dir['curR2R3']:
                        if atomR2R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR2R3)
                    for i in ring3_shareatom:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)

                    for j in threering_list:
                        for n in j:
                            ring_atom_list.append(n)

                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2
                    print(final_H2)

            if molecule == 'C5H10O':
                wt = round((2.016 * float(final_H2) / 86.13) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C6H12O':
                wt = round((2.016 * float(final_H2) / 100.16) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C6H10O':
                wt = round((2.016 * float(final_H2) / 98.14) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C7H14O':
                wt = round((2.016 * float(final_H2) / 114.19) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C7H12O':
                wt = round((2.016 * float(final_H2) / 112.17) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C8H16O':
                wt = round((2.016 * float(final_H2) / 128.21) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C8H14O':
                wt = round((2.016 * float(final_H2) / 126.2) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C8H12O':
                wt = round((2.016 * float(final_H2) / 124.18) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C9H18O':
                wt = round((2.016 * float(final_H2) / 142.24) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C9H16O':
                wt = round((2.016 * float(final_H2) / 140.22) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C9H14O':
                wt = round((2.016 * float(final_H2) / 138.21) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C10H20O':
                wt = round((2.016 * float(final_H2) / 156.26) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C10H18O':
                wt = round((2.016 * float(final_H2) / 154.25) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C10H16O':
                wt = round((2.016 * float(final_H2) / 152.23) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C11H22O':
                wt = round((2.016 * float(final_H2) / 170.29) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C11H20O':
                wt = round((2.016 * float(final_H2) / 168.28) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C11H18O':
                wt = round((2.016 * float(final_H2) / 166.26) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C12H24O':
                wt = round((2.016 * float(final_H2) / 184.32) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C12H22O':
                wt = round((2.016 * float(final_H2) / 182.3) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C12H20O':
                wt = round((2.016 * float(final_H2) / 180.29) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C13H26O':
                wt = round((2.016 * float(final_H2) / 198.34) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C13H24O':
                wt = round((2.016 * float(final_H2) / 196.33) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C13H22O':
                wt = round((2.016 * float(final_H2) / 194.31) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C14H28O':
                wt = round((2.016 * float(final_H2) / 212.37) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C14H26O':
                wt = round((2.016 * float(final_H2) / 210.36) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C14H24O':
                wt = round((2.016 * float(final_H2) / 208.34) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)

        new_df = pd.DataFrame(wt_list)
        wt_df = pd.concat([wt_df, new_df], ignore_index=True)
        print(wt_df)
        return wt_df
        #wt_df.to_excel(f"C{atom_element}-storage.xlsx", sheet_name='Sheet1', index=False)



    elif atom_element == 'N':
        df_N = df

        H2 = 2
        wt_list = []
        wt_df = pd.DataFrame(columns=['IsomericSMILES', 'wt', 'molecule', 'H2num'])
        for index, row in df_N.iterrows():
            threering_list = []
            ring_atom_list = []

            noH_atom_list = []
            final_H2 = 0
            smi = row['IsomericSMILES']
            mol = Chem.MolFromSmiles(smi)
            ssr = Chem.GetSymmSSSR(mol)
            total_atom = mol.GetAtoms()
            noH_atom_idx = noH_atom(mol)
            dups = between_noH(noH_atom_idx, mol)

            H_num = 2 * (len(total_atom) - 1) + 3
            one_cyc = row['IsomericSMILES'].find('1')
            count_1 = row['IsomericSMILES'].count('1')
            two_cyc = row['IsomericSMILES'].find('2')
            three_cyc = row['IsomericSMILES'].find('3')
            if one_cyc != -1 and two_cyc != -1 and three_cyc != -1:
                decrease_H = 2 * 3
            elif one_cyc != -1 and two_cyc != -1 and three_cyc == -1:
                decrease_H = 2 * 2
            elif one_cyc != -1 and two_cyc == -1 and three_cyc == -1:
                if count_1 == 4:
                    decrease_H = 2 * 2
                elif count_1 == 2:
                    decrease_H = 2 * 1
                elif count_1 == 6:
                    decrease_H = 2 * 3
            total_H = H_num - decrease_H
            molecule = f'C{len(total_atom) - 1}H{total_H}N'

            if len(ssr) == 1:
                for ring in ssr:
                    atom_num = len(list(ring)) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = (atom_num - remainder) / 2
                    print(final_H2)
            if len(ssr) == 2:
                ringcon_dir = IsRingConnect(smi)
                if ringcon_dir['R1R2'] == True:
                    shareatom_dir = share_ring_atoms(smi)
                    for ring in ssr:
                        for i in ring:
                            ring_atom_list.append(i)
                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx) - len(shareatom_dir['curR1R2'])
                    remainder = atom_num % H2
                    final_H2 = (atom_num - remainder) / 2
                    print(final_H2)
                else:
                    for ring in ssr:
                        ring_list = list(ring)
                        for i in noH_atom_idx:
                            if i in ring_list:
                                ring_list.remove(i)
                        for j in dups:
                            if j in ring_list:
                                ring_list.remove(j)
                        atom_num = len(ring_list)
                        remainder = atom_num % H2
                        cur_H2 = (atom_num - remainder) / 2
                        final_H2 += cur_H2
                    print(final_H2)
            if len(ssr) == 3:
                ringcon_dir = IsRingConnect(smi)
                if ringcon_dir['R1R2'] == False and ringcon_dir['R2R3'] == False and ringcon_dir['R1R3'] == False:
                    for ring in ssr:
                        ring_list = list(ring)
                        for i in noH_atom_idx:
                            if i in ring_list:
                                ring_list.remove(i)
                        for j in dups:
                            if j in ring_list:
                                ring_list.remove(j)
                        atom_num = len(ring_list)
                        remainder = atom_num % H2
                        cur_H2 = (atom_num - remainder) / 2
                        final_H2 += cur_H2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == True and ringcon_dir['R2R3'] == False and ringcon_dir['R1R3'] == False:
                    shareatom_dir = share_ring_atoms(smi)
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for i in range(0, 2):
                        for j in threering_list[i]:
                            ring_atom_list.append(j)
                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx) - len(shareatom_dir['curR1R2'])
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2

                    for i in noH_atom_idx:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)
                    for j in dups:
                        if j in threering_list[2]:
                            threering_list[2].remove(j)
                    atom_num = len(threering_list[2])
                    remainder = atom_num % H2
                    cur_H2 = (atom_num - remainder) / 2
                    final_H2 += cur_H2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == True and ringcon_dir['R2R3'] == True and ringcon_dir['R1R3'] == True:
                    shareatom_dir = share_ring_atoms(smi)
                    ring3_shareatom = []
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for atomR1R2 in shareatom_dir['curR1R2']:
                        if atomR1R2 in threering_list[1]:
                            threering_list[1].remove(atomR1R2)
                    for atomR1R3 in shareatom_dir['curR1R3']:
                        if atomR1R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR1R3)
                    for atomR2R3 in shareatom_dir['curR2R3']:
                        if atomR2R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR2R3)
                    for i in ring3_shareatom:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)

                    for j in threering_list:
                        for n in j:
                            ring_atom_list.append(n)

                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == True and ringcon_dir['R2R3'] == True and ringcon_dir['R1R3'] == False:
                    shareatom_dir = share_ring_atoms(smi)
                    ring3_shareatom = []
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for atomR1R2 in shareatom_dir['curR1R2']:
                        if atomR1R2 in threering_list[1]:
                            threering_list[1].remove(atomR1R2)
                    for atomR1R3 in shareatom_dir['curR1R3']:
                        if atomR1R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR1R3)
                    for atomR2R3 in shareatom_dir['curR2R3']:
                        if atomR2R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR2R3)
                    for i in ring3_shareatom:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)

                    for j in threering_list:
                        for n in j:
                            ring_atom_list.append(n)

                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == False and ringcon_dir['R2R3'] == True and ringcon_dir['R1R3'] == True:
                    shareatom_dir = share_ring_atoms(smi)
                    ring3_shareatom = []
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for atomR1R2 in shareatom_dir['curR1R2']:
                        if atomR1R2 in threering_list[1]:
                            threering_list[1].remove(atomR1R2)
                    for atomR1R3 in shareatom_dir['curR1R3']:
                        if atomR1R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR1R3)
                    for atomR2R3 in shareatom_dir['curR2R3']:
                        if atomR2R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR2R3)
                    for i in ring3_shareatom:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)

                    for j in threering_list:
                        for n in j:
                            ring_atom_list.append(n)

                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2
                    print(final_H2)

            if molecule == 'C5H11N':
                wt = round((2.016 * float(final_H2) / 85.15) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C6H13N':
                wt = round((2.016 * float(final_H2) / 99.17) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C6H11N':
                wt = round((2.016 * float(final_H2) / 97.16) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C7H15N':
                wt = round((2.016 * float(final_H2) / 113.2) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C7H13N':
                wt = round((2.016 * float(final_H2) / 111.18) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C8H17N':
                wt = round((2.016 * float(final_H2) / 127.23) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C8H15N':
                wt = round((2.016 * float(final_H2) / 125.21) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C8H13N':
                wt = round((2.016 * float(final_H2) / 123.2) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C9H19N':
                wt = round((2.016 * float(final_H2) / 141.25) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C9H17N':
                wt = round((2.016 * float(final_H2) / 139.24) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C9H15N':
                wt = round((2.016 * float(final_H2) / 137.22) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C10H21N':
                wt = round((2.016 * float(final_H2) / 155.28) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C10H19N':
                wt = round((2.016 * float(final_H2) / 153.26) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C10H17N':
                wt = round((2.016 * float(final_H2) / 151.25) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C11H23N':
                wt = round((2.016 * float(final_H2) / 169.31) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C11H21N':
                wt = round((2.016 * float(final_H2) / 167.29) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C11H19N':
                wt = round((2.016 * float(final_H2) / 165.27) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C12H25N':
                wt = round((2.016 * float(final_H2) / 183.33) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C12H23N':
                wt = round((2.016 * float(final_H2) / 181.32) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C12H21N':
                wt = round((2.016 * float(final_H2) / 179.3) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C13H27N':
                wt = round((2.016 * float(final_H2) / 197.36) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C13H25N':
                wt = round((2.016 * float(final_H2) / 195.34) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C13H23N':
                wt = round((2.016 * float(final_H2) / 193.33) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C14H29N':
                wt = round((2.016 * float(final_H2) / 211.39) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C14H27N':
                wt = round((2.016 * float(final_H2) / 209.37) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)
            elif molecule == 'C14H25N':
                wt = round((2.016 * float(final_H2) / 207.35) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule, 'H2num': final_H2}
                wt_list.append(c2)

        new_df = pd.DataFrame(wt_list)
        wt_df = pd.concat([wt_df, new_df], ignore_index=True)
        print(wt_df)
        return wt_df
        #wt_df.to_excel(f"C{atom_element}-storage.xlsx", sheet_name='Sheet1', index=False)


    elif atom_element == None:
        df_C = df

        H2 = 2
        wt_list = []
        wt_df = pd.DataFrame(columns=['IsomericSMILES', 'wt','molecule','H2num'])
        for index,row in df_C.iterrows():
            threering_list = []
            ring_atom_list = []

            noH_atom_list = []
            final_H2 = 0
            smi = row['IsomericSMILES']
            mol = Chem.MolFromSmiles(smi)
            ssr = Chem.GetSymmSSSR(mol)
            total_atom = mol.GetAtoms()
            noH_atom_idx = noH_atom(mol)
            dups = between_noH(noH_atom_idx,mol)

            '''Calculate the theoretical number of hydrogen atoms'''
            H_num = 2 * (len(total_atom) + 1)
            one_cyc = smi.find('1')
            two_cyc = smi.find('2')
            three_cyc = smi.find('3')
            '''Adjust the number of hydrogen atoms according to the number of rings'''
            if one_cyc != -1 and two_cyc != -1 and three_cyc != -1:
                decrease_H = 2 * 3
            elif one_cyc != -1 and two_cyc != -1 and three_cyc == -1:
                decrease_H = 2 * 2
            elif one_cyc != -1 and two_cyc == -1 and three_cyc == -1:
                decrease_H = 2 * 1
            total_H = H_num - decrease_H
            molecule = f'C{len(total_atom)}H{total_H}'

            '''===Calculate the number of H2 according to the number of rings==='''
            if len(ssr) == 1:
                for ring in ssr:
                    '''Calculate the number of atoms on the ring (minus non hydrogen atoms and repeated atoms)'''
                    atom_num = len(list(ring)) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = (atom_num - remainder) / 2
                    print(final_H2)
            if len(ssr) == 2:
                '''Determine whether the ring is connected'''
                ringcon_dir = IsRingConnect(smi)
                if ringcon_dir['R1R2'] == True:
                    shareatom_dir = share_ring_atoms(smi)
                    for ring in ssr:
                        for i in ring:
                            ring_atom_list.append(i)
                    '''Calculate the number of effective atoms (minus non hydrogen, duplicate, and shared atoms)'''
                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx) - len(shareatom_dir['curR1R2'])
                    remainder = atom_num % H2
                    final_H2 = (atom_num - remainder) / 2
                    print(final_H2)
                else:
                    '''Ring not connected'''
                    '''Calculate the amount of Dehydrogenation on each ring separately and then accumulate'''
                    for ring in ssr:
                        ring_list = list(ring)
                        for i in noH_atom_idx:
                            if i in ring_list:
                                ring_list.remove(i)
                        for j in dups:
                            if j in ring_list:
                                ring_list.remove(j)
                        atom_num = len(ring_list)
                        remainder = atom_num % H2
                        cur_H2 = (atom_num - remainder) / 2
                        final_H2 += cur_H2
                    print(final_H2)
            if len(ssr) == 3:
                ringcon_dir = IsRingConnect(smi)
                if ringcon_dir['R1R2'] == False and ringcon_dir['R2R3'] == False and ringcon_dir['R1R3'] == False:
                    for ring in ssr:
                        ring_list = list(ring)
                        for i in noH_atom_idx:
                            if i in ring_list:
                                ring_list.remove(i)
                        for j in dups:
                            if j in ring_list:
                                ring_list.remove(j)
                        atom_num = len(ring_list)
                        remainder = atom_num % H2
                        cur_H2 = (atom_num - remainder) / 2
                        final_H2 += cur_H2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == True and ringcon_dir['R2R3'] == False and ringcon_dir['R1R3'] == False:
                    shareatom_dir = share_ring_atoms(smi)
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for i in range(0,2):
                        for j in threering_list[i]:
                            ring_atom_list.append(j)
                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx) - len(shareatom_dir['curR1R2'])
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2

                    for i in noH_atom_idx:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)
                    for j in dups:
                        if j in threering_list[2]:
                            threering_list[2].remove(j)
                    atom_num = len(threering_list[2])
                    remainder = atom_num % H2
                    cur_H2 = (atom_num - remainder) / 2
                    final_H2 += cur_H2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == True and ringcon_dir['R2R3'] == True and ringcon_dir['R1R3'] == True:
                    shareatom_dir = share_ring_atoms(smi)
                    ring3_shareatom = []
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for atomR1R2 in shareatom_dir['curR1R2']:
                        if atomR1R2 in threering_list[1]:
                            threering_list[1].remove(atomR1R2)
                    for atomR1R3 in shareatom_dir['curR1R3']:
                        if atomR1R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR1R3)
                    for atomR2R3 in shareatom_dir['curR2R3']:
                        if atomR2R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR2R3)
                    for i in ring3_shareatom:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)

                    for j in threering_list:
                        for n in j:
                            ring_atom_list.append(n)

                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == True and ringcon_dir['R2R3'] == True and ringcon_dir['R1R3'] == False:
                    shareatom_dir = share_ring_atoms(smi)
                    ring3_shareatom = []
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for atomR1R2 in shareatom_dir['curR1R2']:
                        if atomR1R2 in threering_list[1]:
                            threering_list[1].remove(atomR1R2)
                    for atomR1R3 in shareatom_dir['curR1R3']:
                        if atomR1R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR1R3)
                    for atomR2R3 in shareatom_dir['curR2R3']:
                        if atomR2R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR2R3)
                    for i in ring3_shareatom:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)

                    for j in threering_list:
                        for n in j:
                            ring_atom_list.append(n)

                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2
                    print(final_H2)
                elif ringcon_dir['R1R2'] == False and ringcon_dir['R2R3'] == True and ringcon_dir['R1R3'] == True:
                    shareatom_dir = share_ring_atoms(smi)
                    ring3_shareatom = []
                    for ring in ssr:
                        threering_list.append(list(ring))
                    for atomR1R2 in shareatom_dir['curR1R2']:
                        if atomR1R2 in threering_list[1]:
                            threering_list[1].remove(atomR1R2)
                    for atomR1R3 in shareatom_dir['curR1R3']:
                        if atomR1R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR1R3)
                    for atomR2R3 in shareatom_dir['curR2R3']:
                        if atomR2R3 not in ring3_shareatom:
                            ring3_shareatom.append(atomR2R3)
                    for i in ring3_shareatom:
                        if i in threering_list[2]:
                            threering_list[2].remove(i)

                    for j in threering_list:
                        for n in j:
                            ring_atom_list.append(n)

                    atom_num = len(ring_atom_list) - len(dups) - len(noH_atom_idx)
                    remainder = atom_num % H2
                    final_H2 = final_H2 + (atom_num - remainder) / 2
                    print(final_H2)

            '''Used to calculate WT and write into the table'''
            if molecule == 'C6H12':
                wt = round((2.016 * float(final_H2) / 84.16) * 100,2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C7H14':
                wt = round((2.016 * float(final_H2) / 98.19) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C7H12':
                wt = round((2.016 * float(final_H2) / 96.17) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C8H16':
                wt = round((2.016 * float(final_H2) / 112.21) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C8H14':
                wt = round((2.016 * float(final_H2) / 110.2) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C9H18':
                wt = round((2.016 * float(final_H2) / 126.24) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C9H16':
                wt = round((2.016 * float(final_H2) / 124.22) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C9H14':
                wt = round((2.016 * float(final_H2) / 122.21) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C10H20':
                wt = round((2.016 * float(final_H2) / 140.27) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C10H18':
                wt = round((2.016 * float(final_H2) / 138.25) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C10H16':
                wt = round((2.016 * float(final_H2) / 136.23) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C11H22':
                wt = round((2.016 * float(final_H2) / 154.29) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C11H20':
                wt = round((2.016 * float(final_H2) / 152.28) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C11H18':
                wt = round((2.016 * float(final_H2) / 150.26) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C12H24':
                wt = round((2.016 * float(final_H2) / 168.32) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C12H22':
                wt = round((2.016 * float(final_H2) / 166.3) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C12H20':
                wt = round((2.016 * float(final_H2) / 164.29) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C13H26':
                wt = round((2.016 * float(final_H2) / 182.35) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C13H24':
                wt = round((2.016 * float(final_H2) / 180.33) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C13H22':
                wt = round((2.016 * float(final_H2) / 178.31) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C14H28':
                wt = round((2.016 * float(final_H2) / 196.37) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C14H26':
                wt = round((2.016 * float(final_H2) / 194.36) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C14H24':
                wt = round((2.016 * float(final_H2) / 192.34) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C15H30':
                wt = round((2.016 * float(final_H2) / 210.4) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C15H28':
                wt = round((2.016 * float(final_H2) / 208.38) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)
            elif molecule == 'C15H26':
                wt = round((2.016 * float(final_H2) / 206.37) * 100, 2)
                c2 = {'IsomericSMILES': row['IsomericSMILES'], 'wt': wt, 'molecule': molecule,'H2num':final_H2}
                wt_list.append(c2)

        new_df = pd.DataFrame(wt_list)
        wt_df = pd.concat([wt_df, new_df], ignore_index=True)
        print(wt_df)
        return wt_df
        #wt_df.to_excel(f"C-storage.xlsx", sheet_name='Sheet1', index=False)


def HSCover6(df):
    data_frame = pd.DataFrame(columns=['IsomericSMILES', 'wt','molecule','H2num'])
    over_six = []
    df1 = df
    for index,row in df1.iterrows():
        if row[1] > 6:
            c2 = {'IsomericSMILES': row[0], 'wt': row[1], 'molecule': row[2],'H2num':row[3]}
            over_six.append(c2)
    new_df = pd.DataFrame(over_six)
    data_frame = pd.concat([data_frame, new_df], ignore_index=True)
    #data_frame.to_excel("storageover6%.xlsx", sheet_name='Sheet1', index=False)

    return data_frame



def underMCH(df_smi,df_bde):
    smiles_dataframe = pd.DataFrame()
    df_smi1 = sc.smitoalfbetsmi(df_smi)
    smiles_list = []
    for smi in df_smi1:
        result = df_bde[df_bde['molecule'] == smi]
        aa_list = []
        for index,row in result.iterrows():
            aa_list.append(row['bde_pred'])
        sums = sum(aa_list)
        bb = sums / len(aa_list)
        if bb <= 96.35616302490234:
            smiles_list.append(smi)


    new_df = pd.DataFrame(smiles_list)
    smiles_dataframe = pd.concat([smiles_dataframe, new_df], ignore_index=True)
    smiles_dataframe.to_excel(f"underMCH.xlsx", sheet_name='Sheet1', index=False)
    return smiles_dataframe