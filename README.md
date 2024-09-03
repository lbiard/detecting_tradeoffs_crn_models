# Detecting context dependent tradeoffs

This repository hosts data and R code for *Bliard L., *Martin JS., Paniw M., Blumstein, D, Martin JGA., Pemberton JM, Nussey DH, Childs DZ, Ozgul A. Detecting context dependence in the expression of life history trade-offs. https://doi.org/10.1111/1365-2656.14173 

*Shared first authorship

Preprint version: https://doi.org/10.32942/X2D31R 

The repository is archived on Zenodo at https://doi.org/10.5281/zenodo.12800618

## GENERAL INFORMATION

1. Title: Data and scripts from "Detecting context-dependence in the expression of life history tradeoffs".

2. Author Information:
	
        A.  Name: Louis Bliard
		Institution: University of Zurich
		Address: Winterthurerstrasse 190, 8057 Zurich, Switzerland
		Email: louis.bliard@uzh.ch
	
        B.  Name: Jordan S Martin
		Institution: University of Zurich
	
        C.  Name: Maria Paniw
		Institution: Estación Biológica de Doñana
		
        D.  Name: Dan Blumstein
		Institution: University of California Los Angeles
		
        E.  Name: Julien GA Martin
		Institution: University of Ottawa
		
        F.  Name: Josephine M Pemberton
		Institution: University of Edinburgh
    
        G.  Name: Dan H Nussey
		Institution: University of Edinburgh
		
        H.  Name: Dylan Z Childs
		Institution: University of Sheffield
		
        I.  Name: Arpat Ozgul
		Institution: University of Zurich
		Email: arpat.ozgul@uzh.ch
		
		
3. Date of data collection: Yellow-bellied marmots - 1979-2020 ; Soay sheep - 1985-2021

4. Geographic location of data collection: Yellow-bellied marmots - 38°57′N, 106°59′W ; Soay sheep - 57°48′N, 8°37′W


## SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC-BY 4.0

2. Links to publications that cite or use the data: https://doi.org/10.32942/X2D31R 

3. Links to other publicly accessible locations of the data: None

4. Was data derived from another source? NA

5. Recommended citation for this dataset: Bliard L., Martin JS., Paniw M., Blumstein, D, Martin JGA., Pemberton JM, Nussey DH, Childs DZ, Ozgul A. (XXXX). Detecting context-dependence in the expression of life history tradeoffs [Data set].



## DATA & FILE OVERVIEW

1. File List: 
- `figure_1_bias.R`
- `functions_postprocessing.R`

- `simulate_data_non_repeated_crn.R`
- `non_repeated_crn_analysis.R`
- `non_repeated_crn_model.stan`

- `simulate_data_hybrid_crn.R`
- `hybrid_crn_analysis.R`
- `hybrid_crn_model.stan`

- `marmot_data_for_analysis.txt`
- `marmot_analysis.R`
- `marmot_model.stan`

- `sheep_data_for_analysis.txt`
- `sheep_analysis.R`
- `sheep_model.stan`

2. Relationship between files, if important: 

The R script `figure_1_bias.R` highlights the amount of bias in estimating among-individual correlations from observation-level correlations. This will produce Figure 1.

The R script `functions_postprocessing.R` contains some functions for postprocessing of the CRN model output.

The script `simulate_data_non_repeated_crn.R` simulates the data needed to then run the Stan model `non_repeated_crn_model.stan` using the script `non_repeated_crn_analysis.R`. This will produce Figure 2.

The script `simulate_data_hybrid_crn.R` simulates the data needed to then run the Stan model `hybrid_crn_model.stan` using the script `hybrid_crn_analysis.R`. This will produce Figure 3.

The data `marmot_data_for_analysis.txt` is needed to run the Stan model `marmot_model.stan` using the script `marmot_analysis.R`. This will produce Figure 4.

The data `sheep_data_for_analysis.txt` is needed to run the Stan model `sheep_model.stan` using the script `sheep_analysis.R`. This will produce Figure 5.


## METHODOLOGICAL INFORMATION
 
1. Methods for processing the data: R. Only the formatted data for analysis is provided. Note that individuals were anonymised in the sheep and marmot datasets (they will not match the individual IDs from other marmot and sheep publications). Please contact the data custodians in charge of the marmot data (Dan Blumstein) or the sheep data (Dan Nussey) if you wish to access the whole marmot or sheep database for your own analyses.

2. Instrument- or software-specific information needed to interpret the data: 
- R v.4.3.2 https://www.r-project.org/
- CmdStanR v.0.6.1 https://mc-stan.org/cmdstanr/

3. People involved with sample collection: Blumstein, D, Martin JGA., Pemberton JM, Nussey DH.

4. People involved with data formatting and analysis: Bliard L., Martin JS., Paniw M., Childs DZ, Ozgul A.

5. for more general informations regarding the methods, see https://doi.org/10.32942/X2D89H and the related code https://github.com/Jordan-Scott-Martin/covariance-reaction-norms


### DATA-SPECIFIC INFORMATION FOR: `marmot_data_for_analysis.txt`

1. Number of variables: 16

2. Number of cases/rows: 2540

3. Variable List: 
- "yrborn" = year of birth of the given marmot.
- "litter_id" = unique ID of the litter the marmot was born in.
- "uid" = unique ID of the marmot.
- "sex" = sex of the marmot (M = male / F = female).
- "col_area" = colony area where the marmot was born.
- "dam" = Unique ID of the mother of the given marmot.
- "pup_emerjdate" = date of emergence of the litter the marmot is part of (in Julian days).
- "littersizeborn" = size of the litter the marmot is part of.
- "massjun" = imputed weight of the marmot at the day of emergence of the litter.
- "massaug" = imputed weight of the marmot on August 31st.
- "age" = age of the mother of the marmot at the time of birth.
- "massjun_mom" = imputed weight of the mother of the marmot on June 1st.
- "date_bare_ground" = first date of bare ground (in Julian days).
- "total_snow" = total amount of snow during the winter (in mm).
- "summer_max_tj" = average daily maximum temperature during the month of June (in degree C).

4. Missing data codes: NA

### DATA-SPECIFIC INFORMATION FOR: `sheep_data_for_analysis.txt`

1. Number of variables: 9

2. Number of cases/rows: 2497

3. Variable List: 
- "obsY" = year of the observation.
- "id" = unique identifier for the ewe (female sheep).
- "ageY" = age of the ewe at the time of the obsevration.
- "capWgt" = weight at the time of capture (summer t).
- "isRepro1" = whether the ewe reproduced (0 = no / 1 = yes).
- "lamdNum1" = number of lamb produced (0 / 1 / 2).
- "growth" = weight at the time of next capture (summer t+1).
- "nao_winter" = overwinter North Atlantic Oscillation.
- "abundance" = number of adult sheep in the population.

4. Missing data codes: NA
