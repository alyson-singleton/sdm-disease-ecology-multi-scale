# sdm-disease-ecology-multi-scale
> Code for: _Species distribution modeling for disease ecology: a multi-scale case study for schistosomiasis host snails in Brazil._ 2024. Singleton, A. L., Glidden, C. K., Chamberlin, A. J., Tuan, R., Palasio, R. G. S., Pinter A., Caldeira, R. L., Mendonça C. L. F., Carvalho O. S., Monteiro, M. V., Athni, T. S., Sokolow S. H., Mordecai, E. A., De Leo, G. A. PLoS Global Public Health (under review). Code written by Alyson Singleton and Caroline Glidden.

__Data__:
* _full_dataset.rds_: All presence and background species data, including coordinates, source, species, subregion, and all covariates used in this analysis. Includes both expert-collected and GBIF data (designated by the "source" column). Made up of:
  1. Expert-collected snail data set, including coordinates, source, and species (_Biomphalaria_ and many others). Compiled from the Coleção de Malacologia Médica, Fundação Oswaldo Cruz (CMM-Fiocruz) and the Coordination for Disease Control of the State Health Secretariat of São Paulo (CCD-SP).
  2. GBIF snail data set, including coordinates, source, and species (_Biomphalaria_ and all others reported in the expert-collected dataset). Collection methods detailed in the manuscript.
* _preprocessed_datasets_: Subsets of the full dataset including presence and background data for each scale and species investigated in this analysis. These datasets have already been thinned and are ready as model input.
  * mg: Minas Gerais, sp: São Paulo, limited: presence points randomly reduced to match the number of total GBIF points for B. tenagophila in São Paulo (Fig E in S1 text).


__Code__:
* Coming soon.
