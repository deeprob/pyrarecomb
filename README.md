# Pythonic version of RareComb
RareComb is a tool to find oligogenic combinations of genes with rare variants that are enriched in individuals with a specific phenotype. RareComb was orginally developed in R (https://github.com/girirajanlab/RareComb). Here we provide a pythonic version of RareComb with some additional utilities.

# Installation
```bash
$ pip install pyrarecomb
```

# User interface 
The pythonic version of RareComb currently has 3 user facing functions:

1. **compare_enrichment**: Checks for oligogenic combinations of rare genetic variants that are enriched in cases but not in controls.

2. **compare_enrichment_depletion**:  Checks for oligogenic combinations of rare genetic variants that are enriched in cases but depleted in controls.

3. **compare_enrichment_modifiers**: Checks for oligogenic combinations of rare genetic variants that are enriched in cases but not in controls where one of the items in a combination must be within an user-defined set of genes.

All these functions have the following **required** arguments:

- *boolean_input_df*: A dataframe where rows are the number of samples and columns include sample ids (represented by the column name: "*Sample_Name*") along with one hot encoded information about the sample genotype (presence or absence rare deleterious mutation within a gene; these columns should start with the prefix "*Input_*") and phenotype (presence or absence of a phenotype; this column should start with the prefix "*Output_*"). Example dataframe is as follows:

Sample_Name | Input_GeneA | Input_GeneB | Input_GeneC | ... | Output_phenotype
--- | --- | --- | --- | --- | --- 
Sample_1111 | 0 | 1 | 1 | ... | 1
Sample_2198 | 0 | 1 | 0 | ... | 0
... 
Sample_N | 0 | 0 | 1 | ... | 0

- *combo_length*: The number of items to mine within a combination.
- *min_indv_threshold*: The minimum number of individuals to consider that must possess a combination before checking for enrichment.
- *max_freq_threshold*: The maximum fraction of the cohort size that possess a combination (to filter out highly frequent combinations).

Along with the other required arguments, **compare_enrichment_modifiers** has an additional required argument:

- *primary_input_entities*: List of genes that must be part of the enriched combinations

All these functions have the following **optional** arguments:

- *input_format*: The prefix of the input columns in the boolean matrix; default="Input_"
- *output_format*: The prefix of the output column in the boolean matrix; default="Output_"
- *pval_filter_threshold*: The p-value significance threshold that the combinations must satisfy; default=0.05
- *adj_pval_type*: The adjusted p-value method to run for multiple testing, one of bonferroni/BH; default="BH"
- *min_power_threshold*: The minimum power threhsold that the significant combinations must satisfy; default=0.7
- *sample_names_ind*: Add samples who possess each combo, one of "Y"/"N"; default="Y"
- *method*: The frequent itemset mining method, one of "fpgrowth"/"apriori"; default="fpgrowth"

# Usage examples
Please refer to the notebooks dir in repo.

# Citation
1. Pounraja VK, Girirajan S. A general framework for identifying oligogenic combinations of rare variants in complex disorders. Genome Res. 2022 May;32(5):904-915. doi: 10.1101/gr.276348.121. Epub 2022 Mar 17. PMID: 35301265; PMCID: PMC9104696.

# Modifications in v0.1.0
## Major
1. Options between apriori and fpgrowth algorithms for frequent itemsets mining

## Minor
1. After filter, raise ValueError check introduced if there is no data)
2. Optional arguments bug fixed

# Possible modifications for v0.2.0
1. Refining control frequencies step may not be required
2. Create function for getting exp and obs prob for combos
3. Create function for calculating p values
4. Discuss the nominal significance filtration strategy
5. Create multiple testing function
6. Rounding adjusted p-values to 3 digits not a good idea
7. compare enrichment modifiers why are we checking for primary entities only as consequents?

# Internal use
## Package creation
```bash
$ python3 -m pip install --upgrade pip
$ python3 -m pip install --upgrade build
$ python3 -m pip install --upgrade twine
$ python -m build
$ python3 -m twine upload dist/*
