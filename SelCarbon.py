import pubchempy as pcp
import pandas as pd
from rdkit import Chem
from alfabet import model


class SiftCarbonRings:
    """Create a class for filtering fully saturated carbocycles."""

    def __init__(self, need_C_num,type_heteroatom=None):
        """An external interface that collects the required number of carbon atoms through this function."""
        self.need_C_num = need_C_num
        self.type_heteroatom = type_heteroatom

    def find_CH(self):
        """Define a function to obtain the molecular formula."""
        C_formula_list = []
        """Create an empty list to store hydrocarbon molecular formulas."""
        C_num = self.need_C_num
        if self.type_heteroatom == None:
            """Obtain the input carbon atoms num."""
            H_num = 2
            """Define the initial number of hydrogen atoms."""
            while H_num <= 2 * C_num:
                C_formula = f'C{C_num}H{H_num}'
                C_formula_list.append(C_formula)
                H_num = H_num + 2
                """Obtain the molecular formula through a loop, 
                where the maximum number of hydrogen atoms is twice the number of carbon atoms, 
                and increase the number of hydrogen atoms by two each time the loop iterates."""
            return C_formula_list
        else:
            heteroatom = self.type_heteroatom[0].upper()
            if heteroatom == 'O' or heteroatom == 'S':
                H_num = 2 * C_num - 4
                while H_num <= 2 * C_num:
                    hetero_ring = f'C{C_num}H{H_num}{heteroatom[0]}'
                    C_formula_list.append(hetero_ring)
                    H_num = H_num + 2
                for None_ring in C_formula_list:
                    if None_ring == 'C13H4O':
                        '''Remove the molecular formula C13H4O, 
                        as substances with this molecular formula are not listed in PubChem, 
                        which may cause bugs.'''
                        C_formula_list.remove(None_ring)
                return C_formula_list
            elif heteroatom == 'N' or heteroatom == 'P':
                H_num = C_num * 2 - 3
                '''If the heteroatom is a nitrogen atom, change the initial number of hydrogen atoms to three.'''
                while H_num <= 2 * C_num + 1:
                    hetero_ring = f'C{C_num}H{H_num}{heteroatom[0]}'
                    C_formula_list.append(hetero_ring)
                    H_num = H_num + 2
                for None_ring in C_formula_list:
                    if None_ring == 'C12H5N' or None_ring == 'C13H3N' or None_ring == 'C15H3N':
                        C_formula_list.remove(None_ring)
                return C_formula_list


    def choosen_CH(self,carbon_list):
        """Define a function to filter carbocycles."""
        carbon_dataframe = pd.DataFrame(columns=['CID', 'IsomericSMILES', 'num'])
        """Create an empty DataFrame containing the keywords CID, SMILES, and num, and store the final results in a DataFrame."""
        C_formula_list = carbon_list
        for one_of_Carbon in C_formula_list:
            Carbon_chains = pcp.get_properties(["isomeric_smiles"], one_of_Carbon, "formula")
            #print(Carbon_chains)
            """Obtain data from the PubChem website, with the obtained data labeled as SMILES, complexity, and CID."""
            str_Brackets = '['
            str_Break = '.'
            str_double_bond = '='
            str_triple_bond = '#'
            """Select the required keywords, remove double (=) and triple (#) bonds, isotopes, radicals ([), and disconnected substances (.) ."""
            isotope_break_removal_list = []
            cycle_list = []
            cycle_list2 = []
            """Create an empty list to store the filtered data."""

            for Carbon_chain in Carbon_chains:
                """Traverse the data obtained from the PubChem website."""
                Isotope_Brackets = Carbon_chain['SMILES'].find(str_Brackets)
                Break_material = Carbon_chain['SMILES'].find(str_Break)
                Double_bond = Carbon_chain['SMILES'].find(str_double_bond)
                Triple_bond = Carbon_chain['SMILES'].find(str_triple_bond)
                """Use the find function to check if the obtained SMILES equation has "[" (representing isotopes, radicals), 
                "." (representing unconnected substances), 
                "=" (representing double bonds), 
                "#" (representing triple bonds). If not, the output value is -1."""
                if Isotope_Brackets == -1 and Break_material == -1 and Double_bond == -1 and Triple_bond == -1:
                    isotope_break_removal_list.append(Carbon_chain)
                """When it is determined that there are no of the above four situations in this SMILES equation, 
                add this SMILES equation to a list."""

            for choose_Carbon_chain in isotope_break_removal_list:
                """Traverse the remaining substances in the previous cycle."""
                carbon_mol_1 = Chem.MolFromSmiles(choose_Carbon_chain['SMILES'])
                """Use RDkit to convert Smiles style to Mol style and perform the following operations."""
                carbon_ssr_1 = Chem.GetSymmSSSR(carbon_mol_1)
                """This function can return a sequence of rings containing atomic IDs."""
                num_ring = len(carbon_ssr_1)
                """Obtain the length of this sequence, which is the number of loops."""
                if num_ring <= 3:
                    cycle_list.append(choose_Carbon_chain)
                """If the number of rings is less than or equal to three, add this SMILES equation to the new list."""

            for cycle_list_casual in cycle_list:
                carbon_mol_2 = Chem.MolFromSmiles(cycle_list_casual['SMILES'])
                carbon_ssr_2 = Chem.GetSymmSSSR(carbon_mol_2)
                ring_atom_num = 0
                """Define the initial ring coefficient as zero."""
                for rings_1 in carbon_ssr_2:
                    """Traverse the loop sequence to ensure that only one atom on the loop is counted at a time."""
                    carbon_rings_1 = len(list(rings_1))
                    """Obtain the number of atoms on the ring."""
                    if carbon_rings_1 > 8 or carbon_rings_1 < 5:
                        ring_atom_num = ring_atom_num + 1
                        """If the number of atoms on the ring is greater than eight or less than five, increase the ring coefficient by one."""
                if ring_atom_num == 0:
                    modified_dict = {
                        'CID': cycle_list_casual['CID'],
                        'IsomericSMILES': cycle_list_casual['SMILES'],
                        'num': cycle_list_casual.get('num', 0)  # 使用get方法避免KeyError
                    }
                    cycle_list2.append(modified_dict)
                    """If the ring coefficient remains zero after the above operation, 
                    it can be determined that the number of atoms on the ring is between five and eight, 
                    and the SMILES equation can be added to the list."""

            if cycle_list2:
                new_df = pd.DataFrame(cycle_list2)
                carbon_dataframe = pd.concat([carbon_dataframe, new_df], ignore_index=True)
                #carbon_dataframe = carbon_dataframe.append(cycle_list2, ignore_index=True)
                """If there is data in the list after the above operation, add the data in this list to the initially created dataframe."""
                print(carbon_dataframe)
        return carbon_dataframe



    def smitoalfbetsmi(self,carbon_dataframe):
        """Convert the SMILES equation obtained from pubchem to the SMILES equation used by alfabet."""
        alfabetsmi_list = []
        dataframe = carbon_dataframe
        for index, row in dataframe.iterrows():
            '''Traverse the data in the columns of the dataframe.'''
            pred_BDE = model.predict([row['IsomericSMILES']])
            '''Use the module provided by alfabet to predict BDE, and use the SMILES formula provided by pubchem.'''
            line = pred_BDE.shape[0]
            '''Determine the number of BDE given.'''
            if line != 1:
                pred_BDE = pred_BDE.drop(index=range(1, line))
                '''If the quantity is greater than one, remove the excess data and only retain one BDE value for future use.'''
            for index2, row2 in pred_BDE.iterrows():
                '''Traverse the data frame for storing and predicting BDE.'''
                carbon = row2['molecule']
                '''Extract the keyword 'molecule' in the form of SMILES, which can correctly obtain the data of alfabet.'''
                alfabetsmi_list.append(carbon)
                '''Put the obtained SMILES into a list.'''
        return alfabetsmi_list


    def carbon_bde(self,alfabetsmi_list):
        """Define a function for obtaining and processing predicted BDE."""
        void_dataframe = pd.DataFrame()
        '''Create an empty dataframe to store BDE.'''

        for smiles in alfabetsmi_list:
            '''Traverse the list of SMILES stored as input strings.'''
            pred_BDE = model.predict([smiles])
            '''Predicting BDE.'''
            pred_BDE_modify = pred_BDE.drop(labels=['bdfe_pred', 'bdfe'], axis=1)
            '''Delete unnecessary data and streamline the structure.'''
            void_dataframe = pd.concat([void_dataframe, pred_BDE_modify], ignore_index=True)
            '''Store data in the created empty dataframe.'''
        pred_BDE_line = void_dataframe.shape[0]
        '''Obtain the number of data.'''
        final_pred_BDE = void_dataframe

        for one_of_predict_bde_line in range(pred_BDE_line):
            '''Traverse each row of data.'''
            predict_bde_series = void_dataframe.iloc[one_of_predict_bde_line]
            predict_bde_keyword = predict_bde_series['bond_type']
            '''Get filtered keywords.'''
            if predict_bde_keyword == 'C-C' or predict_bde_keyword == 'C-N' or predict_bde_keyword == 'H-N' \
                        or predict_bde_keyword == 'C-O' or predict_bde_keyword == 'H-O' or predict_bde_keyword == 'C-P' or predict_bde_keyword == 'H-P' or predict_bde_keyword == 'C-S' or predict_bde_keyword == 'H-S':
                final_pred_BDE = final_pred_BDE.drop(labels=[one_of_predict_bde_line], axis=0)
                '''Determine the keyword, if it is C-C, delete this line.'''
        final_pred_BDE = final_pred_BDE.sort_values(by=['bde_pred'], ascending=True)
        final_pred_BDE = final_pred_BDE.reset_index(drop=True)
        '''Sort the predicted BDE in ascending order and initialize the index with sequence numbers.'''

        return final_pred_BDE

    def none_rings_CH(self,final_pred_BDE):
        """Define a function to remove C-H bonds that are not on the ring."""
        frames = final_pred_BDE
        frames2 = frames
        for index, framerow in frames.iterrows():
            smiles = [framerow['fragment2']]
            '''Traverse the list to obtain column keywords.'''
            mol = Chem.MolFromSmiles(smiles[0])
            carbon_ssr_1 = Chem.GetSymmSSSR(mol)
            num = 0
            '''Set a coefficient equal to zero.'''
            atom_indexs = []
            for atom in mol.GetAtoms():
                atom_indexs.append([(atom.GetSmarts(), atom.GetIdx())])
                '''Getting atoms and their indexes is a nested list.'''
            for atom_index in atom_indexs:
                for one_atom in atom_index:
                    if one_atom[0] == '[C]' or one_atom[0] == '[CH]':
                        for rings_1 in carbon_ssr_1:
                            carbon_rings_1 = list(rings_1)
                            if one_atom[1] in carbon_rings_1:
                                num = num + 1
                                '''Determine whether the atomic ID of the free radical in the given SMILES equation is within
                                 the atomic ID of the ring, and if so, add one coefficient.'''

            if num == 0:
                msk = frames[frames.fragment2 == smiles[0]].index.tolist()
                frames2 = frames2.drop(index=msk, axis=0)
                '''If the coefficient is equal to zero, it means that the free radical position of the SMILES equation in 
                this loop is outside the loop, and the row where this SMILES equation is located is deleted from the dataframe.'''

        return frames2


class SiftCarbonHeteroatom:

    def __init__(self, need_C_num, *type_heteroatom):
        self.need_C_num = need_C_num
        self.heteroatom = type_heteroatom

    def find_C_hetero(self):
        hetero_ring_list = []
        C_num = self.need_C_num
        heteroatom = self.heteroatom
        heteroatom = heteroatom[0].upper()

        if heteroatom == 'O' or heteroatom == 'S':
            H_num = 2 * C_num - 4
            while H_num <= 2 * C_num:
                hetero_ring = f'C{C_num}H{H_num}{heteroatom[0]}'
                hetero_ring_list.append(hetero_ring)
                H_num = H_num + 2
            for None_ring in hetero_ring_list:
                if None_ring == 'C13H4O':
                    hetero_ring_list.remove(None_ring)
            return hetero_ring_list
        elif heteroatom == 'N' or heteroatom == 'P':
            H_num = C_num * 2 - 3
            while H_num <= 2 * C_num + 1:
                hetero_ring = f'C{C_num}H{H_num}{heteroatom[0]}'
                hetero_ring_list.append(hetero_ring)
                H_num = H_num + 2
            for None_ring in hetero_ring_list:
                if None_ring == 'C12H5N' or None_ring == 'C13H3N' or None_ring == 'C15H3N':
                    hetero_ring_list.remove(None_ring)

            return hetero_ring_list

    def choosen_heteroatom_ring(self, hetero_ring_list):
        carbon_dataframe = pd.DataFrame(columns=['CID', 'IsomericSMILES', 'num'])
        C_formula_list = hetero_ring_list

        for one_of_Carbon in C_formula_list:
            Carbon_chains = pcp.get_properties(["isomeric_smiles", "complexity"], one_of_Carbon, "formula")
            str_Brackets = '['
            str_Break = '.'
            str_double_bond = '='
            str_triple_bond = '#'
            isotope_break_removal_list = []
            cycle_list = []
            cycle_list2 = []

            for Carbon_chain in Carbon_chains:
                Isotope_Brackets = Carbon_chain['SMILES'].find(str_Brackets)
                Break_material = Carbon_chain['SMILES'].find(str_Break)
                Double_bond = Carbon_chain['SMILES'].find(str_double_bond)
                Triple_bond = Carbon_chain['SMILES'].find(str_triple_bond)
                if Isotope_Brackets == -1 and Break_material == -1 and Double_bond == -1 and Triple_bond == -1:
                    isotope_break_removal_list.append(Carbon_chain)

            for choose_Carbon_chain in isotope_break_removal_list:
                carbon_mol_1 = Chem.MolFromSmiles(choose_Carbon_chain['SMILES'])
                carbon_ssr_1 = Chem.GetSymmSSSR(carbon_mol_1)
                num_ring = len(carbon_ssr_1)
                if num_ring <= 3:
                    cycle_list.append(choose_Carbon_chain)

            for cycle_list_casual in cycle_list:
                carbon_mol_2 = Chem.MolFromSmiles(cycle_list_casual['SMILES'])
                carbon_ssr_2 = Chem.GetSymmSSSR(carbon_mol_2)
                ring_atom_num = 0
                for rings_1 in carbon_ssr_2:
                    carbon_rings_1 = len(list(rings_1))
                    if carbon_rings_1 > 8 or carbon_rings_1 < 5:
                        ring_atom_num = ring_atom_num + 1
                if ring_atom_num == 0:
                    modified_dict = {
                        'CID': cycle_list_casual['CID'],
                        'IsomericSMILES': cycle_list_casual['SMILES'],
                        'num': cycle_list_casual.get('num', 0)  # 使用get方法避免KeyError
                    }
                    cycle_list2.append(modified_dict)

            print(cycle_list2)

            if cycle_list:
                new_df = pd.DataFrame(cycle_list2)
                carbon_dataframe = pd.concat([carbon_dataframe, new_df], ignore_index=True)
        return carbon_dataframe



    def smitoalfbetsmi(self, carbon_dataframe):
        dataframe = carbon_dataframe
        alfabetsmi_list = []

        for index, row in dataframe.iterrows():
            pred_BDE = model.predict([row['IsomericSMILES']])
            line = pred_BDE.shape[0]
            if line != 1:
                pred_BDE = pred_BDE.drop(index=range(1, line))
            for index2, row2 in pred_BDE.iterrows():
                carbon = row2['molecule']
                alfabetsmi_list.append(carbon)

        return alfabetsmi_list

    def carbon_bde(self, Carbon_list):
        void_dataframe = pd.DataFrame()

        for smiles in Carbon_list:
            pred_BDE = model.predict([smiles])
            pred_BDE_modify = pred_BDE.drop(labels=['bdfe_pred', 'bdfe'], axis=1)
            void_dataframe = pd.concat([void_dataframe, pred_BDE_modify], ignore_index=True)

        pred_BDE_line = void_dataframe.shape[0]
        final_pred_BDE = void_dataframe

        for one_of_predict_bde_line in range(pred_BDE_line):
            predict_bde_series = void_dataframe.iloc[one_of_predict_bde_line]
            predict_bde_keyword = predict_bde_series['bond_type']
            if predict_bde_keyword == 'C-C' or predict_bde_keyword == 'C-N' or predict_bde_keyword == 'H-N' \
                    or predict_bde_keyword == 'C-O' or predict_bde_keyword == 'H-O' or predict_bde_keyword == 'C-P' or predict_bde_keyword == 'H-P' or predict_bde_keyword == 'C-S' or predict_bde_keyword == 'H-S':
                final_pred_BDE = final_pred_BDE.drop(labels=[one_of_predict_bde_line], axis=0)

        final_pred_BDE = final_pred_BDE.sort_values(by=['bde_pred'], ascending=True)
        final_pred_BDE = final_pred_BDE.reset_index(drop=True)

        return final_pred_BDE

    def none_rings_CH(self):
        frames = pd.read_excel(f"C:\\Users\\cwn\\Desktop\\{self.need_C_num}C{self.heteroatom[0]}-bde测试4.xlsx", index_col=None)
        frames2 = frames
        for index, framerow in frames.iterrows():
            smiles = [framerow['fragment2']]
            mol = Chem.MolFromSmiles(smiles[0])
            carbon_ssr_1 = Chem.GetSymmSSSR(mol)
            num = 0
            atom_indexs = []
            for atom in mol.GetAtoms():
                atom_indexs.append([(atom.GetSmarts(), atom.GetIdx())])
            for atom_index in atom_indexs:
                for one_atom in atom_index:
                    if one_atom[0] == '[C]' or one_atom[0] == '[CH]':
                        for rings_1 in carbon_ssr_1:
                            carbon_rings_1 = list(rings_1)
                            if one_atom[1] in carbon_rings_1:
                                num = num + 1

            if num == 0:
                msk = frames[frames.fragment2 == smiles[0]].index.tolist()
                frames2 = frames2.drop(index=msk, axis=0)

        return frames2

def smitoalfbetsmi(carbon_dataframe):
    """Convert the SMILES equation obtained from pubchem to the SMILES equation used by alfabet."""
    alfabetsmi_list = []
    dataframe = carbon_dataframe
    for index, row in dataframe.iterrows():
        '''Traverse the data in the columns of the dataframe.'''
        pred_BDE = model.predict([row['IsomericSMILES']])
        '''Use the module provided by alfabet to predict BDE, and use the SMILES formula provided by pubchem.'''
        line = pred_BDE.shape[0]
        '''Determine the number of BDE given.'''
        if line != 1:
            pred_BDE = pred_BDE.drop(index=range(1, line))
            '''If the quantity is greater than one, remove the excess data and only retain one BDE value for future use.'''
        for index2, row2 in pred_BDE.iterrows():
            '''Traverse the data frame for storing and predicting BDE.'''
            carbon = row2['molecule']
            '''Extract the keyword 'molecule' in the form of SMILES, which can correctly obtain the data of alfabet.'''
            alfabetsmi_list.append(carbon)
            '''Put the obtained SMILES into a list.'''
    return alfabetsmi_list