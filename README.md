# sdm-disease-ecology-multi-scale
> Code for: _Species distribution modeling for disease ecology: a multi-scale case study for schistosomiasis host snails in Brazil._ 2024. Singleton, A. L., Glidden, C. K., Chamberlin, A. J., Tuan, R., Palasio, R. G. S., Pinter A., Caldeira, R. L., Mendonça C. L. F., Carvalho O. S., Monteiro, M. V., Athni, T. S., Sokolow S. H., Mordecai, E. A., De Leo, G. A. PLoS Global Public Health (under review). Code written by Alyson Singleton and Caroline Glidden.

__Data__:
* _full_dataset.rds_: All presence and background species data, including coordinates, source, species, subregion, and all covariates used in this analysis. Includes both expert-collected and GBIF data (designated by the "source" column). Made up of:
  1. Expert-collected snail data set, including coordinates, source, and species (_Biomphalaria_ and many others). Compiled from the Coleção de Malacologia Médica, Fundação Oswaldo Cruz (CMM-Fiocruz) and the Coordination for Disease Control of the State Health Secretariat of São Paulo (CCD-SP).
  2. GBIF snail data set, including coordinates, source, and species (_Biomphalaria_ and all others reported in the expert-collected dataset). Collection methods detailed in the manuscript.
* _preprocessed_datasets_: Subsets of the full dataset including presence and background data for each scale and species investigated in this analysis. These datasets have already been thinned and are ready as model input.
  * _mg_: Minas Gerais, _sp_: São Paulo, _limited_: presence points randomly reduced to match the number of total GBIF points for B. tenagophila in São Paulo (Fig E in S1 text).


__Code__:
* _00_functions.R_: Load pROC functions from the NicheToolBox package in R by Luis Osorio: https://github.com/luismurao/ntbox/blob/master/R/pROC.R.
* _01_preprocessing_species_data.R_: Script used to create _preprocessed_datasets_ above. Posted for reproducibility. Recommend to just use the _preprocessed_datasets_ in the _data_ folder, rather than recreate from scratch.
* _02_maxent_function.R_, _03_rf_function.R_, and _04_brt_function.R_: Main functions used to evaluate model performance (i.e., using spatial cross-validation).
* _05_prediction_bootstrapping.R_: Bootstrapping analysis to produce prediction maps (i.e., random 80/20 splits).
* _06_manuscript_figures.Rmd_: Rmd file including code for constructing all figures, tables, and supplementary information.
