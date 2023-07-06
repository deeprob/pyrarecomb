# Pythonic version of RareComb
RareComb is a tool to find oligogenic combinations of genes with rare variants that are enriched in individuals with a specific phenotype. RareComb is orginally developed in R (https://github.com/girirajanlab/RareComb). Here we provide a pythonic version of RareComb with some additional utilities.

# Installation
```bash
$ pip install pyrarecomb
```

# User interface 
The pythonic version of RareComb currently has 3 user facing functions:

1. **compare_enrichment**: Checks for oligogenic combinations of rare genetic variants that are enriched in cases but not in controls.

2. **compare_enrichment_depletion**:  Checks for oligogenic combinations of rare genetic variants that are enriched in cases but depleted in controls.

3. **compare_enrichment_modifiers**: Checks for oligogenic combinations of rare genetic variants that are enriched in cases but not in controls where one of the items in a combination must be within an user-defined set of genes.

All these functions have the following required arguments:

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

# Usage examples
Please refer to the notebooks dir in repo.

# Possible modifications for v0.1.0
1. Refining control frequencies step may not be required
2. After filter, introduce raise ValueError step if there is no data
3. Create function for getting exp and obs prob for combos
4. Create function for calculating p values
5. Discuss the nominal significance filtration strategy
6. Create multiple testing function
7. Rounding adjusted p-values to 3 digits not a good idea
8. compare enrichment modifiers why are we checking for primary entities only as consequents?

# Internal use
## Package creation
```bash
$ python3 -m pip install --upgrade pip
$ python3 -m pip install --upgrade build
$ python3 -m pip install --upgrade twine
$ python -m build
$ python3 -m twine upload dist/*
