import numpy as np
import pandas as pd
from collections import Counter
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

from .utils.run_apriori_freqitems import run_apriori_freqitems
from .utils import helpers as hp


def compare_enrichment_depletion(
        boolean_input_df, combo_length, min_indv_threshold, max_freq_threshold, 
        input_format="Input_", output_format="Output_", pval_filter_threshold=0.05, 
        adj_pval_type="BH", min_power_threshold=0.7, sample_names_ind="Y", method="fpgrowth"
        ):
    ##########
    # Filter #
    ##########
    apriori_input_cases_df, apriori_input_controls_df, sel_input_colname_list, output_column, number_of_cases = hp.preprocess_boolean(
        boolean_input_df, input_format, output_format, min_indv_threshold, max_freq_threshold
    )
    num_cases, num_controls, num_genes = len(apriori_input_cases_df), len(apriori_input_controls_df), len(sel_input_colname_list)
    # debugging
    print(f"Number of cases remaining after filtration: {num_cases}")
    print(f"Number of controls remaining after filtration: {num_controls}")
    print(f"Number of items remaining after filtration: {num_genes}")

    if min(num_cases, num_controls, num_genes)==0:
        raise ValueError("No samples/items detected: Relax your thresholds")

    ############################
    # CASES / SEVERE Phenotype #
    ############################
    # Introduce a support threshold
    support_threshold = min_indv_threshold / apriori_input_cases_df.shape[0]
    case_freqitems_df, case_freqitems_size1_df = run_apriori_freqitems(apriori_input_cases_df, combo_length, support_threshold, method=method)
    # set the number of frequent items column name 
    case_freqitems_df = case_freqitems_df.rename(columns={"Obs_Count_Combo": "Case_Obs_Count_Combo"})
    # get the number of unique items forming combinations
    uniq_combo_items = set(case_freqitems_df.loc[:, [f"Item_{i}" for i in range(1, combo_length+1)]].values.flatten())
    # Store the counts as a dictionary for each item
    case_freqitems_countdict = dict(zip(case_freqitems_size1_df.Item_1, case_freqitems_size1_df.Obs_Count_Combo.astype(int)))
    # Get the observed count for each item
    for i in range(1, combo_length+1):
        case_freqitems_df[f"Case_Obs_Count_I{i}"] = case_freqitems_df[f"Item_{i}"].map(case_freqitems_countdict)
    # Get the expected probability of observing the combos
    case_freqitems_df["Case_Exp_Prob_Combo"] = case_freqitems_df.loc[:, [f"Case_Obs_Count_I{i}" for i in range(1, combo_length+1)]].prod(axis=1)/(number_of_cases**combo_length)
    # Get the observed probability of observing the combos
    case_freqitems_df['Case_Obs_Prob_Combo'] = case_freqitems_df['Case_Obs_Count_Combo'] / number_of_cases
    # Using bionomial test, calculate p-value
    case_freqitems_df['Case_pvalue_more'] = case_freqitems_df.apply(lambda row: binomtest(int(row['Case_Obs_Count_Combo']), number_of_cases, row['Case_Exp_Prob_Combo'], alternative='greater').pvalue, axis=1)
    # debugging
    print(f'Number of initial combinations identified for cases: {case_freqitems_df.shape[0]}')
    print(f'Number of unique items in cases: {len(uniq_combo_items)}')

    #############################
    # CONTROLS / MILD Phenotype #
    #############################
    # Get the control profile of the combo items
    apriori_input_controls_df = apriori_input_controls_df.loc[:, list(uniq_combo_items)]
    number_of_controls = apriori_input_controls_df.shape[0]
    # define support threshold for controls
    support_threshold = 2 / number_of_controls
    # get the frequently mutated genes in controls using apriori
    cont_freqitems_df, cont_freqitems_size1_df = run_apriori_freqitems(apriori_input_controls_df, combo_length, support_threshold, method=method)
    # set the number of frequent items column name 
    cont_freqitems_df = cont_freqitems_df.rename(columns={"Obs_Count_Combo": "Cont_Obs_Count_Combo"})
    # Store the counts as a dictionary for each item
    cont_freqitems_countdict = dict(zip(cont_freqitems_size1_df.Item_1, cont_freqitems_size1_df.Obs_Count_Combo.astype(int)))
    # Keep combos found in case only
    case_cont_freqitems_df = case_freqitems_df.merge(cont_freqitems_df, left_on="uniq_items", right_on="uniq_items", how="left", suffixes=('', '_cont')).drop(columns=[f"Item_{i}_cont" for i in range(1, combo_length + 1)]).fillna(0.)
    # Get the observed count in controls for each item
    for i in range(1, combo_length+1):
        case_cont_freqitems_df[f"Cont_Obs_Count_I{i}"] = case_cont_freqitems_df[f"Item_{i}"].map(cont_freqitems_countdict)
    # Get the expected probability of observing the combos in controls
    case_cont_freqitems_df["Cont_Exp_Prob_Combo"] = case_cont_freqitems_df.loc[:, [f"Cont_Obs_Count_I{i}" for i in range(1, combo_length+1)]].prod(axis=1)/(number_of_controls**combo_length)
    # Get the observed probability of observing the combos
    case_cont_freqitems_df['Cont_Obs_Prob_Combo'] = case_cont_freqitems_df['Cont_Obs_Count_Combo'] / number_of_controls
    # Using bionomial test, calculate p-value
    case_cont_freqitems_df['Cont_pvalue_less'] = case_cont_freqitems_df.apply(lambda row: binomtest(int(row['Cont_Obs_Count_Combo']), number_of_controls, row['Cont_Exp_Prob_Combo'], alternative='less').pvalue, axis=1)
    # debugging
    print(f"Number of controls: {number_of_controls}")
    print(f"Number of combinations with support of at least 2 in controls: {cont_freqitems_df.shape[0]}")

    ########################
    # Nominal significance #
    ########################
    # TODO: This step is not omitted from compare_enrichment_depletion since we
    # TODO: consider all for multiple testing.
    sel_case_cont_freqitems_df = case_cont_freqitems_df
    # debugging
    print(f"Number of combinations considered for multiple testing correction: {sel_case_cont_freqitems_df.shape[0]}")

    ####################
    # Multiple testing #
    ####################
    # Create variable for number of tests done
    number_of_tests = sel_case_cont_freqitems_df.shape[0]
    # multiple test BH and Bonferroni - round to 3 places of decimal will change in later versions
    sel_case_cont_freqitems_df['Case_Adj_Pval_bonf'] = np.round(multipletests(sel_case_cont_freqitems_df['Case_pvalue_more'].values, method='bonferroni')[1], 3)
    sel_case_cont_freqitems_df['Case_Adj_Pval_BH'] = np.round(multipletests(sel_case_cont_freqitems_df['Case_pvalue_more'].values, method='fdr_bh')[1], 3)
    # add a column for number of tests done
    sel_case_cont_freqitems_df['Num_tests'] = number_of_tests
    # filter significant items
    if adj_pval_type == 'BH':
        all_sig_case_cont_freqitems_df = sel_case_cont_freqitems_df[
            (sel_case_cont_freqitems_df['Case_Adj_Pval_BH'] < pval_filter_threshold) &
            (sel_case_cont_freqitems_df['Cont_pvalue_less'] < pval_filter_threshold)
        ]
    elif adj_pval_type == 'bonferroni':
        all_sig_case_cont_freqitems_df = sel_case_cont_freqitems_df[
            (sel_case_cont_freqitems_df['Case_Adj_Pval_bonf'] < pval_filter_threshold) &
            (sel_case_cont_freqitems_df['Cont_pvalue_less'] < pval_filter_threshold)
        ]
    multtest_sig_comb_count = all_sig_case_cont_freqitems_df.shape[0]
    # debugging
    print(f"Number of combinations that are significant after multiple testing correction: {multtest_sig_comb_count}")

    ###################
    # Post processing #
    ###################
    # Check if there is at least a single significant combination after multiple testing correction
    if multtest_sig_comb_count > 0:
        all_sig_case_cont_freqitems_df = hp.refine_control_frequencies(all_sig_case_cont_freqitems_df, apriori_input_controls_df,number_of_controls,check_enrichment=False)

        ######################
        # POWER CALCULATIONS #
        ######################
        all_sig_case_cont_freqitems_df = hp.calculate_power(all_sig_case_cont_freqitems_df, number_of_cases, number_of_controls)
        output_sig_case_cont_freqitems_df = all_sig_case_cont_freqitems_df.loc[all_sig_case_cont_freqitems_df.Power_Five_Pct >= min_power_threshold]
        
        #####################
        # SAMPLES DETECTION #
        #####################
        if len(output_sig_case_cont_freqitems_df)>0:
            print(f"Number of significant combinations that meet the power threshold is {len(output_sig_case_cont_freqitems_df)}")
            if sample_names_ind == "Y":
                output_sig_case_cont_freqitems_df = hp.add_sample_info(boolean_input_df, output_sig_case_cont_freqitems_df, output_column)

        
        else:
            print("No significant combinations that meet the specified power threshold")
            print("Returning ONLY the non-significant combinations")
            output_sig_case_cont_freqitems_df = sel_case_cont_freqitems_df

    else:
        print("No significant combinations were found after multiple testing correction")
        print("Returning ONLY the non-significant combinations")
        output_sig_case_cont_freqitems_df = sel_case_cont_freqitems_df

    return output_sig_case_cont_freqitems_df


if __name__ == "__main__":
    # load the input df
    boolean_input_df = pd.read_csv("/data6/deepro/computational_pipelines/pyrarecomb/test/input/test_input.csv")
    # define all other params
    combo_length = 2
    min_indv_threshold = 5
    max_freq_threshold = 0.25

    compare_enrichment_depletion(
        boolean_input_df, combo_length, min_indv_threshold, max_freq_threshold
        )
