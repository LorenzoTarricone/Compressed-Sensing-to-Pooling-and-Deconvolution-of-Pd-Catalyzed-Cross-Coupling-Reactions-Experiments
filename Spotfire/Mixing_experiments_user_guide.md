# Guided workflow for mixed experiments in Spotfire
This small guide aims to give all the instructions to set up the experiment and later choose the parameters to perform the analysis of mixed experiments. Each input parameter will be discussed in its importance for determining the output, and some practical suggestions will be given to facilitate the choice of the hyperparameter
## Pooling Algorithm / Plate Builder
**Input:**
 - `cats_list` (column of Table): the list of catalysts we want to test for in our mixed experiment

**Output:**
 - `pooling_df` (Table): data frame containing the instructions for the pooling. Each row is a pool to use in the experiment, and the seven entries in the columns correspond to the different catalysts to use in that mixture

**REMARK:** Not any arbitrary number of catalysts can be tested in a mixed experiment. To use a (partial or full) Kirkman matrix as a pooling matrix, the number of catalysts to be tested needs to be either a **multiple of 5 (up to 35)** or a **multiple of 7 (up to 70)**. An error message should come out in the debugging tool if the length of the input catalyst list is not correct

## Run the experiment
All the details to run your mixed experiments or replicate ours are available in our paper _Compressed Sensing Approach to Pooling and Deconvolution of Pd-Catalyzed Cross-Coupling Reactions_ by L. Tarricone et Al.

## Deconvolution Algorithm / Results Analysis
**Input:** (for standard workflow)
- `data_raw` (Table): data frame containing the result of the mixed experiments. It is assumed that the data in this table come from the same ELN. To be used correctly by the algorithm. This has to have the following columns:
     1. ELN_ID: (String) Specific ELN identification number
     2. SAMPLE_ID: (String) Unique identifier of the sample (vial)
     3. PLATENUMBER: (Integer) Plate number inside this ELN
     4. SAMPLENUMBER: (Integer) Sample number (sample timepoint)
     5. PLATEROW: (String) Capital letter indicating the row (from A to H)
     6. PLATECOLUMN: (Integer) Number indicating the column (from 1 to 12)
     7. COMPOUND:  (String) Type of compound considered (e.g. "Product", "lim SM", "unknown", etc...)
     8. ColumnLabel: (String) Name of the solvent and base combination used throughout that row.
     9. AREA_TOTAL_REDUCED: (Real) The area percentage that we have discarded non-interesting elements. We calculate it from AREA_PERCENT as `[AREA_TOTAL] * (100 / (100 - sum(If(([COMPOUND]="Internal Standard") Or ([COMPOUND]="Solvent") Or ([COMPOUND]="Ligand") Or ([COMPOUND]="ignore Peak"),[AREA_TOTAL],0)) OVER ([SAMPLE_ID],[SAMPLENUMBER])))`
- `cats_names`: (Column)  	Cats names needed to find again the matrix that was used for the pooling
- `Time_correction`: (Boolean) Switch that chooses if applying time correction in the algorithm. Time correction assumes that the product is relatively stable throughout sampling points. If you suspect that the product considered was not stable, you should set this to False
-  `SM_weights`: (Boolean) Switch that chooses if applying starting material weighting in the algorithm. SM weighting assumes that more importance should be given to sampling points where a higher fraction of vials have some limiting starting material inside (independently of the quantity). It is suggested to look qualitatively at the result of the mixed experiment and see if there are sampling points with a lot of Product localised in a few experiments (good signal/noise ratio) and no limiting starting material at all. In that case, it might be wise to turn this switch off to not have a small weight on a very informative time point.
-  `toll_exp`: When selecting the chemical conditions (solvent base combinations /rows ) to be considered for the average, we want to retain the signal that comes from conditions with enough non-zero entries. discard the conditions where there are fewer than min_non_zero_exp non-zero experiments. It is suggested to decrease this value (therefore being more generous with the solvent/base conditions to consider) when, maybe in hard reactions, a lot of mixed experiments are zero
-  `plate_number`: (Integer) Switch that depends on the Document property that filters the ELN for the plate to be analysed with the mixed experiment
-  `chemical_prior`: (Boolean) Switch that allows you to choose if you want to use the chemical prior (embeddings) for the reconstruction. It is suggested to set this to true when a lot of product is found across the plate.

**Input** (needed when using the chemical prior)
- `M_E_embeddings_neutral`: (Table) Chemical prior input (corrected for typos and missing bits). Pandas data frame of dimension num_cats x embedding_dim with row index given by cats names (SPELLING SHOULD MATCH THE ONES IN cats_names ). This embedding should contain a vector for ALL the cats in cats_names
-  `M_E_embeddings_bh`: (Table) Chemical prior input similar to M_E_embeddings_neutral but adjusted with a contrastive learning algorithm using Buchwald-Hartwig reactions.
-  `M_E_embeddings_sm`: (Table) Chemical prior input similar to M_E_embeddings_neutral but adjusted with a contrastive learning algorithm using Suzuki-Miyaura reactions.
-  `rxn_category`: (String) Tells the algorithm which prior to use (if used), options are `bh`, `sm` and `neutral`.
- `lambdaa`: (Real) lambda parameter for the weight (positive) of the prior in the optimisation problem. A higher value of lambda means more weight is given to the embeddings (our prior knowledge), and therefore, the resulting good catalysts will tend to be from the same region of the embedding space. A value equal to zero is equivalent to the case without a chemical prior. It is suggested to start by experimenting with a small value (around 1) and then progressively increase it and rerun the analysis. No substantial difference is observed in the result empirically when 
- `n_hitters_to_return`: (Integer) Number of top hitters (in order of score) that we will be returned from the compressed sensing algorithm. If the number of hitters (catalysts where the reconstructed score is more than 0.1) is less than this number, then all the hitters will be returned.

  **Output:**
- `chem_cond_df`: (Table) Dataframe giving stats on the area% found in the different conditions. As a rule of thumb, having a condition with a small minimum and high maximum indicates a condition where the mixed experiments gave a good signal/noise ratio
- `result_final`: (Table) Dataframe containing a ranking of the top 12 catalysts found in this experiment that can be tested in the follow-up 24 vials plate
- `df_perc_of_sm`: (Table) Dataframe giving information about the percentage of vials in the late with some limiting starting material (for each sample point). As a rule of thumb, high ratios (when there are also some vials with larger area%) indicate that in many vials the catalysts were all exhausted, and the measure we record of area% is not a lower bound but the actual value
- `other_df`: (Table) Dataframe giving information about the amounts of other materials (not "Produc"t or " lim SM" or "other SM") at each sample point. High values are to be considered as an indicator of bad quality of the experiment, and indicate that many unknown compounds are present in the vial
- `df_reshaped_plate`: (Table) Dataframe giving the result of the experimental plate in the shape: N_conditions X N_ mixtures
- `df_complete_data`: (Table) Dataframe giving the result of the deconvolution (normalised so that maximum score is 100) in the shape of N_conditions X N_catalysts. Good catalysts are the ones that have high scores across all conditions.

## Other information

- Every time you have changed the parameters and you want to rerun the analysis, you can go to the "three connected squares visualisation" (bottom left). Select the data function you are interested in and click the button `Refresh` at the top
- You can change the Python script by choosing `Data > Data Functions Properties > Edit Script...`
- You can change the parameters values/origin by choosing `Data > Data Functions Properties > Edit Parameters...`
