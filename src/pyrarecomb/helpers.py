import numpy as np
import pandas as pd
from collections import Counter
from scipy.stats import binomtest
from statsmodels.stats.power import tt_ind_solve_power
from statsmodels.stats.proportion import proportion_effectsize


def preprocess_boolean(boolean_input_df, input_format, output_format, min_indv_threshold, max_freq_threshold):
    """
    This function parses and filters user-given boolean dataframe to rarecomb
    1) It defines input and output columns
    2) It creates case and control boolean matrix for data mining algorithm
    3) It filters the matrix to only keep relevant items based on user defined conditions
    """
    # Identify all the input and output variables
    input_colname_list = [col for col in boolean_input_df.columns if col.startswith(input_format)]
    output_column = [col for col in boolean_input_df.columns if col.startswith(output_format)][0]
    # Make cases and controls apriori input df
    apriori_input_cases_df = boolean_input_df.loc[boolean_input_df[output_column]==1, input_colname_list].applymap(int)
    apriori_input_controls_df = boolean_input_df.loc[boolean_input_df[output_column]==0, input_colname_list].applymap(int)
    # Get the number of cases from input param - max freq thresh
    number_of_cases = apriori_input_cases_df.shape[0]
    max_instances = round(number_of_cases * max_freq_threshold)
    # Filter case and control df to only include gene cols that satisfy input criterion
    apriori_input_cases_df = apriori_input_cases_df.loc[:, (apriori_input_cases_df.sum() >= min_indv_threshold) 
                                                        & (apriori_input_controls_df.sum() >= 1) & (apriori_input_cases_df.sum() < max_instances)]
    apriori_input_controls_df = apriori_input_controls_df.loc[:, (apriori_input_cases_df.sum() >= min_indv_threshold) 
                                                              & (apriori_input_controls_df.sum() >= 1) & (apriori_input_cases_df.sum() < max_instances)]
    # Select columns that remain
    sel_input_colname_list = [col for col in apriori_input_cases_df.columns if col.startswith(input_format)]
    return apriori_input_cases_df, apriori_input_controls_df, sel_input_colname_list, output_column, number_of_cases 

def get_counts(uniq_items, input_df):
    """
    This function gets the counts of combos from an input boolean df
    """
    query = " & ".join([f"({i} == 1)" for i in uniq_items.split("|")])
    return len(input_df.query(query))

def refine_control_frequencies(all_sig_case_cont_freqitems_df, apriori_input_controls_df,number_of_controls,check_enrichment=True):
    # Check if zero frequency cases exist in controls
    cont_combos_w_zero_freq_df = all_sig_case_cont_freqitems_df.loc[all_sig_case_cont_freqitems_df["Cont_Obs_Count_Combo"] == 0]
    zero_freq_combo_count = cont_combos_w_zero_freq_df.shape[0]
    print('Number of combinations with support less than 2 in controls:', zero_freq_combo_count)
    if zero_freq_combo_count > 0:
        # REFINE CONTROL FREQUENCIES
        # for the zero frequency combos in controls, get their actual combo size
        cont_combos_w_zero_freq_df["Cont_Obs_Count_Combo"] = cont_combos_w_zero_freq_df.uniq_items.apply(get_counts, args=(apriori_input_controls_df, ))
        cont_combos_w_zero_freq_df["Cont_Obs_Prob_Combo"] = cont_combos_w_zero_freq_df['Cont_Obs_Count_Combo'] / number_of_controls
        if check_enrichment:
            cont_combos_w_zero_freq_df['Cont_pvalue_more'] = cont_combos_w_zero_freq_df.apply(lambda row: binomtest(int(row['Cont_Obs_Count_Combo']), number_of_controls, row['Cont_Exp_Prob_Combo'], alternative='greater').pvalue, axis=1)
        else:
            # check for depletion
            cont_combos_w_zero_freq_df['Cont_pvalue_less'] = cont_combos_w_zero_freq_df.apply(lambda row: binomtest(int(row['Cont_Obs_Count_Combo']), number_of_controls, row['Cont_Exp_Prob_Combo'], alternative='less').pvalue, axis=1)
        # rewrite the items in the main df
        all_sig_case_cont_freqitems_df.update(cont_combos_w_zero_freq_df)
    return all_sig_case_cont_freqitems_df

def calculate_power(all_sig_case_cont_freqitems_df, number_of_cases, number_of_controls):
    all_sig_case_cont_freqitems_df['Case_Exp_Count_Combo'] = round((all_sig_case_cont_freqitems_df['Case_Exp_Prob_Combo'] * number_of_cases), 2)
    all_sig_case_cont_freqitems_df['Cont_Exp_Count_Combo'] = round((all_sig_case_cont_freqitems_df['Cont_Exp_Prob_Combo'] * number_of_controls), 2)
    all_sig_case_cont_freqitems_df['Effect_Size'] = all_sig_case_cont_freqitems_df.apply(
        lambda row: proportion_effectsize(row['Case_Obs_Prob_Combo'], row['Cont_Obs_Prob_Combo']), axis=1)
    all_sig_case_cont_freqitems_df['Power_One_Pct'] = round(all_sig_case_cont_freqitems_df.apply(
        lambda row: tt_ind_solve_power(effect_size=row['Effect_Size'], nobs1=number_of_cases, ratio=number_of_controls/number_of_cases, alpha=0.01), axis=1), 3)
    all_sig_case_cont_freqitems_df['Power_Five_Pct'] = round(all_sig_case_cont_freqitems_df.apply(
        lambda row: tt_ind_solve_power(effect_size=row['Effect_Size'], nobs1=number_of_cases, ratio=number_of_controls/number_of_cases, alpha=0.05),axis=1), 3)
    return all_sig_case_cont_freqitems_df

def get_samples(row, samples_df, output_column):
    """
    Helper function for adding sample information
    """
    items = row.uniq_items.split("|")
    count_cases = Counter(samples_df.loc[(samples_df.Items.isin(items))&(samples_df[output_column]==1)].Sample_Name)
    case_samples = [s for s,ns in count_cases.items() if ns==len(items)]
    count_controls = Counter(samples_df.loc[(samples_df.Items.isin(items))&(samples_df[output_column]==0)].Sample_Name)
    control_samples = [s for s,ns in count_controls.items() if ns==len(items)]
    return pd.Series({"Case_Samples": "|".join(case_samples), "Control_Samples": "|".join(control_samples)})

def add_sample_info(boolean_input_df, output_sig_case_cont_freqitems_df, output_column):
    # add case and control samples for each combo
    samples_df = boolean_input_df.set_index(["Sample_Name", output_column])
    samples_df = samples_df.mask(samples_df == 0).stack().reset_index().drop(0, axis=1).rename(columns={"level_2": "Items"})
    output_sig_case_cont_freqitems_df = output_sig_case_cont_freqitems_df.merge(output_sig_case_cont_freqitems_df.apply(get_samples, args=(samples_df, output_column), axis=1), left_index=True, right_index=True)
    return output_sig_case_cont_freqitems_df
