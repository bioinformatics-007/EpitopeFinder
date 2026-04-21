# modules/vaccine_assembly.py
import itertools
import pandas as pd

order_dict = {1: [0, 1, 2], 2: [0, 2, 1], 3: [2, 0, 1], 4: [2, 1, 0], 5: [1, 2, 0], 6: [1, 0, 2]}

def generate_combinations(strings, num):
    return list(itertools.permutations(strings, num))

def final_components(seqs, linker):
    final_comp = []
    for combo in seqs:
        str_temp = linker.join(combo)
        final_comp.append(linker + str_temp)
    return final_comp

def epitope_combinations(sequence, list_epitopes):
    return [sequence + i for i in list_epitopes]

def combination_function(combination, epitopes):
    new_list_with_all_epitopes = []
    for new_epitope in epitopes:
        temp_list = list(combination) if isinstance(combination, (tuple, list)) else [combination]
        temp_list.append(new_epitope)
        new_list_with_all_epitopes.append(temp_list)
    return new_list_with_all_epitopes

def epitopes_on_vaccine(epitopes_all, order, nums):
    all_epitopes_combined = []
    temp_nums = list(nums)
    for i in range(3):
        if len(epitopes_all[i]) == 0: temp_nums[i] = 0

    for j in range(len(order)):
        category_idx = order[j]
        epitopes = epitopes_all[category_idx]
        if not epitopes: continue
        if not all_epitopes_combined:
            all_epitopes_combined = [e for e in epitopes]
        else:
            new_list = []
            for existing_combo in all_epitopes_combined:
                new_list.extend(combination_function(existing_combo, epitopes))
            all_epitopes_combined = new_list
    
    names = ['B-cell', 'CTL', 'HTL']
    df_data = {n: [] for n in names}
    for row in all_epitopes_combined:
        curr_idx = 0
        for cat_idx in order:
            name = names[cat_idx]
            if temp_nums[cat_idx] > 0:
                val = row[curr_idx] if isinstance(row[curr_idx], (list, tuple)) else [row[curr_idx]]
                df_data[name].append(list(val))
                curr_idx += 1
            else:
                df_data[name].append("N/A")
    return df_data

def multiepitope_vaccine(components, order, start='', end=''):
    subunit = []
    for category_idx in order:
        comp_list = components[category_idx]
        if not comp_list: continue
        if not subunit:
            subunit = [start + c for c in comp_list]
        else:
            new_ls = []
            for vaccine_base in subunit:
                new_ls.extend(epitope_combinations(vaccine_base, comp_list))
            subunit = new_ls
            
    # Add C-terminal (His-tag)
    return [s + end for s in subunit] if subunit else [start + end]

def run_assembly(dfs, num_epitopes, linkers, order_opt, n_term, c_term):
    b_df, c_df, h_df = dfs
    components = [[], [], []]
    comb_seqs = [[], [], []]

    for i, (df, linker, num) in enumerate(zip([b_df, c_df, h_df], linkers, num_epitopes)):
        if df is not None and not df.empty and num > 0:
            col = next((c for c in df.columns if c.lower() in ['sequences', 'peptide', 'sequence', 'matched sequence']), None)
            if col:
                seqs = df[col].dropna().unique()[:6] 
                actual_num = min(num, len(seqs))
                if actual_num > 0:
                    comb = generate_combinations(seqs, actual_num)
                    comb_seqs[i] = comb
                    components[i] = final_components(comb, linker)
    
    order = order_dict.get(int(order_opt), [0, 1, 2])
    vacc_sequences = multiepitope_vaccine(components, order, n_term, c_term)
    epit_data = epitopes_on_vaccine(comb_seqs, order, num_epitopes)
    
    final_df = pd.DataFrame({'Vaccine_Candidate': vacc_sequences})
    final_df['B_cell_Epitopes'] = epit_data['B-cell']
    final_df['CTL_Epitopes'] = epit_data['CTL']
    final_df['HTL_Epitopes'] = epit_data['HTL']
    return final_df
