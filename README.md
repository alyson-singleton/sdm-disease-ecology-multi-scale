# sdm-disease-ecology-multi-scale
> Code for: [_Species distribution modeling for disease ecology: a multi-scale case study for schistosomiasis host snails in Brazil._](https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0002224) 2024. Singleton, A. L., Glidden, C. K., Chamberlin, A. J., Tuan, R., Palasio, R. G. S., Pinter A., Caldeira, R. L., Mendonça C. L. F., Carvalho O. S., Monteiro, M. V., Athni, T. S., Sokolow S. H., Mordecai, E. A., De Leo, G. A. PLOS Global Public Health. Code written by Alyson Singleton and Caroline Glidden.
>
> Please feel free to email me with any specific, code-related questions at asinglet@stanford.edu.

__Data__:
* _full_dataset.rds_: All presence and background species data, including coordinates, source, species, subregion, and all covariates used in this analysis. Includes both expert-collected and GBIF data (designated by the "source" column):
  1. Expert-collected snail data set, including coordinates, source, and species (_Biomphalaria_ and many others). Compiled from the Coleção de Malacologia Médica, Fundação Oswaldo Cruz (CMM-Fiocruz) and the Coordination for Disease Control of the State Health Secretariat of São Paulo (CCD-SP).
  2. GBIF snail data set, including coordinates, source, and species (_Biomphalaria_, all species reported in the expert-collected dataset, and all freshwater animals found in South America, as defined by the International Union for Conservation of Nature). Collection methods detailed in the manuscript.
* _preprocessed_datasets_: Subsets of the full dataset including presence and background data for each scale and species investigated in this analysis. These datasets have already been thinned and are ready as model input.
  * _mg_: Minas Gerais, _sp_: São Paulo, _limited_: presence points randomly reduced to match the number of total GBIF points for _B. tenagophila_ in São Paulo (Fig E in S1 text).
* _mapping_datasets_: Datasets used to create Fig 1 maps.


__Code__:
* _00_functions.R_: Load pROC functions from the NicheToolBox package in R by Luis Osorio: https://github.com/luismurao/ntbox/blob/master/R/pROC.R.
* _01_preprocessing_species_data.R_: Script used to create _preprocessed_datasets_ above. Posted for reproducibility. Recommend to just use the _preprocessed_datasets_ in the _data_ folder, rather than recreate from scratch.
* _02_maxent_function.R_, _03_rf_function.R_, and _04_brt_function.R_: Main functions used to evaluate model performance (i.e., using spatial cross-validation), create partial dependence plots, and calculate variable importance.
* _05_prediction_bootstrapping.R_: Bootstrapping analysis to produce prediction maps (i.e., random 80/20 splits).
* _06_manuscript_figures_tables.Rmd_: Rmd file including code for constructing all figures and tables in the main manuscript.
* _07_supplementary_figures_tables.Rmd_: Rmd file including code for constructing all figures and tables in the supplementary information (S1_text).
