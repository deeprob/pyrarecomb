from functools import reduce
import pandas as pd
from mlxtend.frequent_patterns import apriori, fpgrowth, association_rules


def get_freq_items_combo(frequent_itemsets, combo_length):
    # get items of specific combo size
    frequent_itemsets = frequent_itemsets.loc[frequent_itemsets.length==combo_length].drop(columns=["length"])
    # Split itemsets into separate columns
    frequent_itemsets = frequent_itemsets.merge(frequent_itemsets['itemsets'].apply(lambda x: pd.Series(list(x))), left_index=True, right_index=True).drop(columns=["support", "itemsets"])
    # rename columns
    old_colnames = range(combo_length)
    new_colnames = [f"Item_{i}" for i in range(1, combo_length+1)]
    frequent_itemsets = frequent_itemsets.rename(columns=dict(zip(old_colnames, new_colnames)))
    return frequent_itemsets.loc[:, new_colnames + ["uniq_items", "Obs_Count_Combo"]].reset_index(drop=True)

def add_frozensets(a, b):
    return a.union(b)  

def run_apriori_freqitems(apriori_input_df, combo_length, support_threshold, primary_entities=None, method="fpgrowth"):
    method_dict = {"apriori": apriori, "fpgrowth": fpgrowth}
    frequent_itemsets = method_dict[method](
        apriori_input_df.astype(bool), min_support=support_threshold, 
        use_colnames=True, max_len=combo_length
        )
    if primary_entities is not None:
        # run association rules
        assoc_df = association_rules(frequent_itemsets, metric="confidence", min_threshold=0.0)
        # the primary entities must be in the consequents only
        assoc_df = assoc_df.loc[
            (assoc_df.consequents.apply(lambda x: True if len(x.intersection(primary_entities))>0 else False)) &
            (assoc_df.antecedents.apply(lambda x: True if len(x.intersection(primary_entities))==0 else False))
            ]
        # convert assoc df to frequent items
        assoc_df["itemsets"] = assoc_df.apply(lambda row: row["antecedents"].union(row["consequents"]), axis=1)
        # get the individual items from frequent items
        frequent_itemsets = frequent_itemsets.loc[frequent_itemsets.itemsets.apply(lambda x: len(x)==1)]
        # add the filtered items obtained from association mining
        frequent_itemsets = pd.concat((frequent_itemsets, assoc_df.loc[:, ["itemsets", "support"]]))

    frequent_itemsets['count'] = frequent_itemsets['support'] * len(apriori_input_df)
    frequent_itemsets['Obs_Count_Combo'] = frequent_itemsets.pop('count')
    frequent_itemsets['length'] = frequent_itemsets['itemsets'].apply(lambda x: len(x))
    frequent_itemsets['uniq_items'] = frequent_itemsets['itemsets'].apply(lambda x: "|".join(sorted(x)))
    frequent_itemsets_len_combo = get_freq_items_combo(frequent_itemsets, combo_length)
    frequent_itemsets_len_1 = get_freq_items_combo(frequent_itemsets, 1)
    return frequent_itemsets_len_combo, frequent_itemsets_len_1
