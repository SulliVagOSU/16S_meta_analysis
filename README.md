# Cervicovaginal microbiome of people living with HIV in DRC
Code accompanying the manuscript "The cervicovaginal microbiome of pregnant people living with HIV on antiretroviral therapy in the Democratic Republic of Congo: A Pilot Study and Global Meta-analysis" ![image](https://github.com/user-attachments/assets/7dff1da4-b098-4e02-984d-4ab79a6260ac)


## Contributors
Kimberley S Ndlovu<br/> 
Ricardo R Pavan, PhD

## File descriptions:
### DATA folder
contains the .RDS or .RDATA phyloseq files 

ps_congo_filt_spp is phyloseq  of CVMB count data and taxa from the CQI-PMTCT (CONGO) samples that is filtered and agglomerated to species level<br/> 
ps_congo_filt_spp is phyloseq of CVMB count data and taxa from the CQI-PMTCT (CONGO) samples that is filtered, agglomerated to species level and transformed using robust centered log ratio<br/>  
ps_congo_picrust is a phyloseq of predicted metabolic pathways for the CQI-PMTCT (CONGO) samples<br/> 
ps_filt_spp_clr_zero is a phyloseq  of CVMB count data and taxa from the intergated datasets that is filtered, zero sum ASVs removed, and transformed using centered log ratio<br/> 
ps_filt_spp_clr_zero_adj is a phyloseq  of CVMB count data and taxa from the intergated datasets that is filtered, transformed using centered log ratio and batch adjusted<br/> 
ps_filt_spp_zero is a ps_filt_spp_clr_zero phyloseq  of CVMB count data and taxa from the intergated datasets that is filtered and transformed using centered log ratio<br/> 
ps_meta_picrust is a phyloseq of predicted metabolic pathways for integrated datasets<br/>  
ps_raw_all is a phyloseq of raw count data for the integrated datasets<br/> 
### CODE folder
DRC_CVMB_Figures2.Rmd contains all scripts used to make the figures in the paper. Also includes analyses used to make the inputfiles for the figures<br/> 
functions_drc_cvmb.R functions used in the DRC_CVMB_Figures2.Rmd<br/> 
dada2_congo_samples.R contains dada2 analyses for the CQI-PMTCT (CONGO) samples<br/> 
batch_adj_plots.Rmd contains batch adjustment scripts and figures for the integrated datasets<br/> 
### TABLES folder
contains all the input files for the code/DRC_CVMB_Figures2.Rmd. Scripts to obtian the files are included in the same Rmd file.
