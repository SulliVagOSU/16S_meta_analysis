#MODIFIED BY: KSNdlovu
#DATE LASTMODIFIED: March 12, 2025
#ORIGINAL AUTHOR: RRPavan

colours <- readRDS("~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Possible thesis projects/HIV project/Vaginal study w Ann&Ricardo/16S_data_vaginal_DNA/Attempt_2/scripts/colours.RDS")

########DATA PROCESSING######
#assign level 1,2,4 pathways to a picrust 2 output
categorize_by_function_l3 <- function(in_ko, kegg_brite_mapping) {
  # Function to categorize by KEGG Pathway levels, including Level 1, 2, and 3 groupings.
  # Input KO table is assumed to have rownames as KOs and sample names as columns.
  
  # Initialize the output data frame with columns for Level 1, Level 2, and Level 3
  out_pathway <- data.frame(matrix(NA, nrow=0, ncol=(ncol(in_ko) + 3)))
  colnames(out_pathway) <- c("Level1", "Level2", "Level3", colnames(in_ko))
  
  for(ko in rownames(in_ko)) {
    
    # Skip KO if not in KEGG BRITE mapping df
    if(! ko %in% rownames(kegg_brite_mapping)) {
      next
    }
    
    # Extract the list of pathways for the current KO
    pathway_list <- strsplit(kegg_brite_mapping[ko, "metadata_KEGG_Pathways"], "\\|")[[1]]
    
    for(pathway in pathway_list) {
      
      # Extract Level 1, Level 2, and Level 3 from the pathway string
      levels <- strsplit(pathway, ";")[[1]]
      level1 <- levels[1]  # Level 1 pathway
      level2 <- levels[2]  # Level 2 pathway
      level3 <- levels[3]  # Level 3 grouping
      
      # Create a new row with NA and KO-specific data
      new_row <- data.frame(matrix(c(NA, NA, NA, as.numeric(in_ko[ko,])), nrow=1, ncol=ncol(out_pathway)))
      colnames(new_row) <- colnames(out_pathway)
      
      # Assign Level 1, Level 2, and Level 3 to the new row
      new_row$Level1 <- level1
      new_row$Level2 <- level2
      new_row$Level3 <- level3
      
      # Append the new row to the output data frame
      out_pathway <- rbind(out_pathway, new_row)
    }
  }
  
  # Aggregate data by Level 1, Level 2, and Level 3
  out_pathway <- data.frame(aggregate(. ~ Level1 + Level2 + Level3, data = out_pathway, FUN=sum))
  
  # Set row names to a combination of Level 1, Level 2, and Level 3
  rownames(out_pathway) <- paste(out_pathway$Level1, out_pathway$Level2, out_pathway$Level3, sep = "_")
  
  # Remove the Level 1, Level 2, and Level 3 columns if not needed anymore
  out_pathway <- out_pathway[, -which(colnames(out_pathway) %in% c("Level1", "Level2", "Level3"))]
  
  # Remove rows with zero sums
  if(length(which(rowSums(out_pathway) == 0)) > 0) {
    out_pathway <- out_pathway[-which(rowSums(out_pathway) == 0), ]
  }
  
  return(out_pathway)
}
#agglomerate phyloseq
glom_ps <- function(physeq, level = level, NArm=TRUE) {
  physeq <- tax_glom(physeq, taxrank = level, NArm=NArm)
  return(physeq)
}
#remove outliers using the 1.5IQR methods
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}
#rename phyloseq object to have taxa names in the otu and remove tax table
rename_taxa_phyloseq <- function(ps) {
  # Extract OTU table from phyloseq object
  
  otu <- otu_table(ps) #taxa must be rows
  if (!taxa_are_rows(ps)) {
    otu <-t(otu)
  }
  
  otu <- as.data.frame(otu)
  # Get the taxonomic annotations from the phyloseq object
  taxa <- tax_table(ps) %>% as.data.frame()  %>% rownames_to_column("asvid") %>% mutate(Taxon = case_when(
    !is.na(Species) ~ paste(sub("^g__", "", Genus), sub("^s__", "", Species), sep = " "),
    !is.na(Genus) ~ (Genus),
    !is.na(Family) ~ Family,
    !is.na(Order) ~ Order,
    !is.na(Class) ~ Class,
    !is.na(Phylum) ~ Phylum,
    TRUE ~ "Unknown"
  )) %>%select(-Kingdom, -Species, -Genus, -Family, -Order, -Class, -Phylum)
  
  # Rename ASVs to actual taxa names
  renamed_otu_table <- otu %>%
    rownames_to_column(var = "asvid") %>%
    left_join(., as.data.frame(taxa), by = "asvid") %>%
    select(-asvid) %>%
    column_to_rownames(var = "Taxon") %>%
    otu_table(taxa_are_rows = T)
  
  # Update the OTU table in the phyloseq object
  ps2 <- phyloseq(otu_table(renamed_otu_table), sample_data(sample_data(ps)))
  
  # Return the updated phyloseq object
  return(ps2)
}
#faster psmel
psmelt.dplyr = function(physeq, glom_by = NULL) {
  sd = data.frame(sample_data(physeq)) %>% rownames_to_column("Sample")
  TT = data.frame(tax_table(physeq)) %>% rownames_to_column("OTU")
  otu.table = data.frame(otu_table(physeq), check.names = FALSE) %>% 
    rownames_to_column("OTU")
  if(!is.null(glom_by)){
    otu.table %>% 
      data.table() %>%
      setDT() %>%
      melt(id.vars = "OTU", variable.name = "Sample", value.name = "abundance") %>%
      lazy_dt() %>%
      left_join(TT) %>%
      group_by(.data[[glom_by]], Sample) %>%
      summarise(Abundance = sum(Abundance)) %>%
      left_join(sd)  %>%
      collect()
  }else{
    otu.table %>% 
      data.table() %>%
      setDT() %>%
      melt(id.vars = "OTU", variable.name = "Sample", value.name = "abundance") %>%
      lazy_dt() %>%
      left_join(sd) %>%
      left_join(TT) %>%
      collect()
  }
}

#' Subsample.Table from https://github.com/jbisanz/MicrobeR/blob/master/R/Subsample.Table.R 
#'
#' @description Takes a Feature/OTU/SV table of counts where samples are columns and feature names are row names and returns a table with even coverage on a per-sample basis. Now this function is really just an alias for rarefy_even_depth from phyloseq or rtk from rtk. NOTE: Sampling with replacement for sinlge rarefraction method!
#'
#' @param OTUTABLE Table of feature/OTU/SV counts where Samples are columns, and IDs are row names
#' @param DEPTH Count depth, defaults to min(colSums(OTUTABLE)) if not passed
#' @param SEED A randomization SEED, defaults to 182. This is ignored for multiple subsamples and instead the seeds 1:NSAMPS is used.
#' @param NSAMPS The number of samples that should be taken of the table for multiple subsamples. Use an odd number to avoid decimals. Default=101
#' @param THREADS Number of CPUs to use for multiple rarefraction, defaults to 2/3 of available.
#' @param VERBOSE Should progress and metrics be printed to screen via message()? Default=TRUE
#' @return Subsampled Table
#' @export

Subsample.Table<-function(FEATURES,DEPTH, SEED, VERBOSE){
  if(missing(VERBOSE)){VERBOSE=T}
  if(missing(DEPTH)){DEPTH=min(colSums(FEATURES))}
  if(missing(SEED)){SEED=182}

  if(VERBOSE==T){message(paste("Subsampling feature table to", DEPTH, ", currently has ", nrow(FEATURES), " taxa."))}
  subsampled.FEATURES<-as.data.frame(phyloseq::rarefy_even_depth(otu_table(FEATURES, taxa_are_rows = T), sample.size=DEPTH, rngseed=SEED, verbose = FALSE)) #expecting the transpose for otu table layout so transpose then transpose back
  if(VERBOSE==T){message(paste("...sampled to",DEPTH, "reads with", nrow(subsampled.FEATURES), "taxa"))}


  return(subsampled.FEATURES)
}

Multiple.Subsample.Table<-function(FEATURES,DEPTH, VERBOSE, NSAMPS, THREADS, AGGFUNCTION){
  if(missing(VERBOSE)){VERBOSE=T}
  if(missing(DEPTH)){DEPTH=min(colSums(FEATURES))}
  if(missing(NSAMPS)){NSAMPS=101}
  
  if(NSAMPS%%2==0){stop("NSAMPS must be an odd number to prevent fractions in read counts.")}
  
  if(missing(THREADS)){THREADS=round(parallel::detectCores()*2/3, 0)}
  if(missing(AGGFUNCTION)){AGGFUNCTION="median"}
  
  if(!AGGFUNCTION %in% c("median","mean")){stop("Aggregation function (AGGFUNCTION) must be either mean or median")}
  
  if(VERBOSE){message(paste("Subsampling feature table to", DEPTH, ", currently has ", nrow(FEATURES), " taxa."))}
  
  tables<-rtk::rtk(FEATURES, repeats=NSAMPS, depth=DEPTH, threads=THREADS, ReturnMatrix = NSAMPS, verbose=FALSE, margin=2)$raremat
  
  if(VERBOSE){message("Merging individual features tables:")}
  subsampled.FEATURES<-matrix(nrow=nrow(FEATURES), ncol=ncol(FEATURES))
  rownames(subsampled.FEATURES)<-rownames(FEATURES)
  colnames(subsampled.FEATURES)<-colnames(FEATURES)
  for(i in 1:ncol(FEATURES)){
    if(VERBOSE){message(i/ncol(FEATURES)*100, "%") }
    subsampled.FEATURES[,i]<-  
      lapply(tables, function(x) x[,i]) %>%
      do.call(cbind, .) %>%
      apply(., 1, AGGFUNCTION)
  }
  
  subsampled.FEATURES<-subsampled.FEATURES[rowSums(subsampled.FEATURES)>0,]
  
  if(VERBOSE){message(paste("...multiply sampled to",DEPTH, "with the median feature count reported. A total of",nrow(subsampled.FEATURES), "taxa have been returned."))}
  if(VERBOSE){print(summary(colSums(subsampled.FEATURES)))}
  
  return(subsampled.FEATURES)
}

#######TRANSFORMATIONS########
# Zero imputation and CLR transformation 
# references Sugnet Lubbe, Peter Filzmoser, Matthias Templ (2021)
# code based on https://github.com/thomazbastiaanssen/Tjazi/blob/master/R/clr_lite.R

#' Impute zeroes and perform a centered log-ratio (CLR) transformation
#' @description Microbiome data is compositional. When compositional data is examined using non-compositional methods, many problems arise.
#' Performing a centered log-ratio transformation is a reasonable way to address these problems reasonably well.
#' \cr \cr
#' A major problem with this approach is that microbiome data typically contains lots of zeroes and the logarithm of zero is undefined.
#' Here, we implemented a few methods discussed by Lubbe \emph{et al.} 2021 to replace zeroes with non-zero values in such a way that the structure of the data remains reasonably well preserved.
#' \cr \cr
#' Some of these methods (namely 'logunif' and 'runif') involve imputing small values between 0 and the lowest non-zero value in the dataset.
#' For these methods, we have implemented a resampling approach in order to stabilize the inter-run variability. See \code{method} for more information.
#' @param counts A compositional count table.
#' @param samples_are Either "cols" or "rows". Default is "cols". Denotes whether the columns or rows depict individual samples.
#' @param method The method for zero imputation. One of \code{"logunif"}, \code{"unif"} or \code{"const"}.
#' \code{'logunif'} samples small numbers from a log-uniform distribution, whereas \code{'unif'} samples from a uniform one. On the other hand, \code{"const"} simply replaces zeroes with \code{0.65 * [the lowest value]}.
#' @param replicates An integer. For the two random sampling methods, if this is larger than 1, every zero will be imputed that many times. The median of the CLR of all those replicates will be returned. If \code{method} is set to \code{"const"}, replicates will be automatically set to 1 as no random numbers are generated.
#' @return A CLR-transformed count table.
#' @references Sugnet Lubbe, Peter Filzmoser, Matthias Templ (2021)
#' \emph{Comparison of zero replacement strategies for compositional data with large numbers of zeros.}
#' doi:https://doi.org/10.1016/j.chemolab.2021.104248
#' @export
clr_lite = function(counts, samples_are = "cols", method = "logunif", replicates = 1000)
{
  temp_counts = counts
  
  if(! method %in% c("logunif", "unif", "const"))
  {stop("`method` must be exactly `logunif`, `unif` or `const`")}
  
  if(method == "const"){replicates = 1}
  
  if(samples_are == "rows"){
    temp_counts = data.frame(t(temp_counts))
  }
  
  temp_counts = apply(X          = temp_counts,
                      MARGIN     = 2,
                      FUN        = clr_imputed,
                      method     = method,
                      replicates = replicates)
  
  if(samples_are == "rows"){
    temp_counts = data.frame(t(temp_counts))
  }
  
  clr_counts = data.frame(temp_counts)
  rownames(clr_counts) = rownames(counts)
  colnames(clr_counts) = colnames(counts)
  return(clr_counts)
}

#' compute CLR using Aitchison's method
#' @description See \code{clr_lite}.
#' @seealso \code{\link{clr_lite}}
#' @param x A vector of compositional data without zeroes.
#' @return A vector of CLR-transformed data
#'
anansi_compute_clr = function(x){
  #compute CLR using Aitchison's method
  return(log(x/exp(mean(log(x)))))
}

#' Replace zeroes with non-zero values in order to perform a CLR-transformation
#' @description See \code{\link{clr_lite}}.
#' @seealso \code{\link{clr_lite}}
#' @param vec A vector of compositional data that may contain zeroes.
#' @param method The method for zero imputation. One of "logunif", "unif" or "const".
#' @return A vector with all the zeroes replaced with non-zero values.
#' @importFrom stats runif
#'
impute_zeroes = function(vec, method = "logunif"){
  if(! method %in% c("logunif", "unif", "const")){stop("`method` must be exactly `logunif`, `unif` or `const`")}
  
  #Find detection limit
  DL = min(vec[vec != 0])
  if(method == "logunif"){
    vec[vec == 0] = DL/(10^(runif(n = sum(vec == 0), min =  0, max = 1)))
  }
  else if(method == "unif"){
    vec[vec == 0] = runif(n = sum(vec == 0), min =  0.1*DL, max = DL)
  }
  else if(method == "const"){
    vec[vec == 0] = 0.65 * DL
    print(paste0("The DL is equal to ", DL))
  }
  return(vec)
}

#' Resample random values, perform CLR over each iteration and return the median result.
#' @description See \code{clr_lite}.
#' @seealso \code{\link{clr_lite}}
#' @param vec A vector of compositional data that may contain zeroes.
#' @param method The method for zero imputation. One of "logunif", "unif" or "const".
#' @param replicates A positive integer. Default is 1000. Controls how many replicates the median should be taken over.
#' @return a vector of CLR-transformed data
#' @importFrom stats median
#'
clr_imputed = function(vec, method = "logunif", replicates = 1000){
  if(! method %in% c("logunif", "unif", "const")){stop("`method` must be exactly `logunif`, `unif` or `const`")}
  return(apply(replicate(replicates, anansi_compute_clr(impute_zeroes(vec = vec, method = method))), 1, median))
}

#' @rdname clr_lite
#' @section Functions:
#' \code{clr_c:}
#'  A wrapper for \code{clr_lite(counts, method = "const", replicates = 1)}.
#' @export
#'
clr_cc <- function(counts, samples_are = "cols"){
  clr_lite(counts, samples_are = samples_are, method = "const", replicates = 1)
} 
clr_cr <-function(counts, samples_are = "rows"){
  clr_lite(counts, samples_are = samples_are, method = "const", replicates = 1)
} 

#' @rdname clr_lite
#' @section Functions:
#' \code{clr_unif:}
#'  A wrapper for \code{clr_lite(counts, method = "unif")}.
#' @export
#'
clr_unif <- function(counts, samples_are = "cols", replicates = 1000){
  clr_lite(counts, samples_are = samples_are, method = "unif", replicates = replicates)
}

#' @rdname clr_lite
#' @section Functions:
#' \code{clr_logunif:}
#'  A wrapper for \code{clr_lite(counts, method = "logunif")}.
#' @export
#'
clr_logunif <- function(counts, samples_are = "cols", replicates = 1000){
  clr_lite(counts, samples_are = samples_are, method = "logunif", replicates = replicates)
}

#hellinger transformation taxa_are_rows = T
hellinger_transform <- function(ps) {
  otu_table_ps <- otu_table(ps)
  otu_table_ps <- t(otu_table_ps)
  hellinger_transformed <- decostand(as.matrix(otu_table_ps), "hellinger")
  ps_new <- phyloseq(otu_table(hellinger_transformed, taxa_are_rows = T), sample_data(ps), phy_tree(ps), refseq(ps), tax_table(ps))
  return(ps_new)
}
#robust clr taxa_are_rows = F
robust_clr_transform <- function(ps) {
  otu_table_ps <- otu_table(ps)
  rclr_transformed <- decostand(as.matrix(otu_table_ps), "rclr")
  ps_new <- phyloseq(otu_table(rclr_transformed, taxa_are_rows = T), sample_data(ps), phy_tree(ps), refseq(ps), tax_table(ps))
  return(ps_new)
}
#same as above but omits the phy_tree
robust_clr_transform2 <- function(ps) {
  otu_table_ps <- otu_table(ps)
  rclr_transformed <- decostand(as.matrix(otu_table_ps), "rclr")
  ps_new <- phyloseq(otu_table(rclr_transformed, taxa_are_rows = T), sample_data(ps), refseq(ps), tax_table(ps))
  return(ps_new)
}

#log transformation taxa_are_rows=T
log_transform <- function(ps) {
  otu_table_ps <- as.data.frame(otu_table(ps))
  otu_table_ps <- t(otu_table_ps)
  log_transformed <- decostand(as.matrix(otu_table_ps), "log")
  ps_new <- phyloseq(otu_table(log_transformed, taxa_are_rows = T), sample_data(ps), phy_tree(ps), refseq(ps), tax_table(ps))
  return(ps_new)
}
#sqrt root taxa_are_rows = F
sqrt_transform <- function(ps) {
  otu_table_ps <- otu_table(ps)
  sqrt_transformed <- sqrt(as.matrix(otu_table_ps))
  ps_new <- phyloseq(otu_table(sqrt_transformed, taxa_are_rows = T), sample_data(ps),phy_tree(ps), refseq(ps), tax_table(ps))
  return(ps_new)
}
#z transformatoion x is an otu table
 z_transform <- function(x) {
   z <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)  # Z transformation
   return(z)
 }
#phylogenetic isometric log
 philr_transform <- function(ps) {
   otu_table_ps <- as.data.frame(otu_table(ps))
   ilr_transformed <- decostand(as.matrix(otu_table_ps), "ilr")
   ps_new <- phyloseq(
     otu_table(ilr_transformed, taxa_are_rows = TRUE),
     sample_data(ps),
     phy_tree(ps),
     tax_table(ps)
   )
   return(ps_new)
 }

#additive log taxa_are_rows = T
alr_transform <- function(ps) {
  otu_table_ps <- as.data.frame(otu_table(ps))
  otu_table_ps <- t(otu_table_ps)
  alr_transformed <- decostand(as.matrix(otu_table_ps), "alr")
  ps_new <- phyloseq(otu_table(alr_transformed, taxa_are_rows = T), sample_data(ps),phy_tree(ps), refseq(ps), tax_table(ps))
  return(ps_new)
}

#cubed root taxa_are_rows = F
cubed_root_transform <- function(ps) {
  otu_table_ps <- otu_table(ps)
  cubed_root_transformed <- (as.matrix(otu_table_ps))^(1/3)
  ps_new <- phyloseq(otu_table(cubed_root_transformed, taxa_are_rows = T), sample_data(ps),phy_tree(ps), refseq(ps), tax_table(ps))
  return(ps_new)
}

########ORDINATIONS#######
# Principal component analysis (PCA) with density plots per componen
Scatter_Density_batch <- function(data = data,
                                  data_meta = data_meta,
                                  batch = batch, 
                                  trt = trt, 
                                  scale = F, 
                                  batch.legend.title = batch.legend.title, 
                                  trt.legend.title = trt.legend.title, 
                                  title = title,
                                  batch_color = batch_color,
                                  trt_color = trt_color){
  
  # set color
  # friendly_cols <- dittoSeq::dittoColors()
  # pre-processing
  df <- data
  df_meta <- meta_data %>% mutate(across(where(is.character), as.factor))
  # n = n_distinct(df_meta[batch])
  
  # fit PCA
  pca_fit <- df %>% 
    prcomp(center = T, scale. = scale) 
  # tidy the results
  var <- pca_fit %>% 
    tidy(matrix = "eigenvalues")
  
  pcs_fit <- pca_fit %>% 
    tidy(matrix = "pcs")
  
  perc_expl_var <- pcs_fit %>%
    filter(PC < 6) %>%
    ggplot(aes(x = PC, y = percent)) +
    geom_bar(stat = "identity")
  
  pca_df_variances <- pca_fit %>%
    augment() %>%
    rename("sample_name" = 1) %>%
    left_join(df_meta) %>%
    rename_with(~str_replace(.,".fitted",""))
  
  # make the plots
  pMain <- pca_df_variances %>%
    ggplot(aes(PC1, PC2, color = .data[[batch]], shape = .data[[trt]])) + 
    geom_point(alpha = 0.7, size =3) + 
    #geom_text(check_overlap = TRUE, hjust = 'inward', family = 'IBM Plex Sans') +
    geom_hline(yintercept=0, linetype="dashed", alpha = 0.3) +
    geom_vline(xintercept=0, linetype="dashed", alpha = 0.3) +
    labs(x = paste("PC1", round(var$percent[1]*100,2), "%"),
         y = paste("PC2", round(var$percent[2]*100,2), "%")) +
    # scale_colour_viridis_d(option = "plasma") 
    scale_color_brewer(palette = "Paired") +
    scale_fill_brewer(palette = "Paired")+
    scale_shape_manual(values = seq_along(unique(pca_df_variances[[trt]]))) +
    theme_bw() +
    labs(colour = batch.legend.title, shape = trt.legend.title)
  
  pTop <- pca_df_variances %>%
    ggplot(aes(PC1, fill = .data[[batch]], linetype = .data[[trt]])) + 
    geom_density(size = 0.2, alpha = 0.5) + 
    ylab('Density - PC1') + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_text(size = rel(0.8)), 
          plot.title = element_text(hjust = 0.5, size = rel(1.5)), 
          axis.line = element_blank(), 
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') + 
    scale_color_brewer(palette = "Paired") +
    labs(title = title)
  
  pRight <- pca_df_variances %>%
    ggplot(aes(PC2, fill = .data[[batch]], linetype = .data[[trt]])) + 
    geom_density(size = 0.2, alpha = 0.5) + 
    coord_flip() +
    ylab('Density - PC2') + 
    theme(axis.title.x = element_text(size = rel(0.8)), 
          axis.title.y = element_blank(), 
          axis.line = element_blank(),
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') + 
    scale_color_brewer(palette = "Paired")
  
  # plot all together
  pTop + theme(legend.position = 'none') + guide_area() + pMain  + pRight + theme(legend.position = 'none')  + plot_layout(guides = 'collect') +
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4)) & theme(legend.box = "horizontal")
}

#normal PCA
run_PCA <- function(ps = ps,
                    batch = batch, 
                    scale = FALSE, 
                    title = title,
                    color_palette = color_palette){
  
  #df <- t(otu_table(ps)) #if taxa_are_rows =T
  df <- otu_table(ps)
  
  if (taxa_are_rows(ps)) {
    df <-t(df)
  }
  
  df_meta <- sample.data.frame(ps) 
  if (!("sampleID" %in% colnames(df_meta))) {
    df_meta <- df_meta %>%
      mutate(across(where(is.character), as.factor))  %>% 
      rownames_to_column("sampleID")
  } else {
    df_meta <- df_meta %>%
      mutate(across(where(is.character), as.factor))
  }
  
  # Fit PCA
  pca_fit <- prcomp(na.omit(df), center = TRUE, scale. =scale)
  
  # Tidy the results
  var <- pca_fit %>% 
    tidy(matrix = "eigenvalues")
  
  pcs_fit <- pca_fit %>% 
    tidy(matrix = "pcs")
  
  perc_expl_var <- pcs_fit %>%
    filter(PC < 6) %>%
    ggplot(aes(x = PC, y = percent)) +
    geom_bar(stat = "identity")
  
  pca_df_variances <- pca_fit %>%
    augment() %>%
    dplyr::rename("sampleID" = 1) %>%
    left_join(df_meta) %>%
    rename_with(~str_replace(.,".fitted",""))
  
  # Make the plots
  # num_colors <- length(unique(df_meta[[batch]]))
  # colours <- color_palette[1:num_colors]
  pca_df_variances <- pca_df_variances %>% filter(!is.na(.data[[batch]]))
  
  ord_plot <- pca_df_variances %>%
    ggplot(aes(PC1, PC2, color = as.factor(.data[[batch]]))) + 
    geom_point(alpha = 0.7, size = 3) + 
    #geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    #geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(x = paste("PC1", round(var$percent[1]*100, 2), "%"),
         y = paste("PC2", round(var$percent[2]*100, 2), "%")) +
    scale_color_manual(values = color_palette) + theme_minimal() + theme(panel.background = element_rect(fill = NA,colour = "white"),  panel.border = element_blank(), axis.line = element_line()) + ggtitle(title) + guides(color = guide_legend(ncol = 2))
  
  return(ord_plot)
}
# PCA biplot
run_PCA_biolot <- function(ps = ps,
                    batch = batch, 
                    scale = FALSE, 
                    title = title,
                    color_palette = color_palette) {
  
  # Extract the OTU table
  df <- otu_table(ps)
  
  # Transpose OTU table if taxa are rows
  if (taxa_are_rows(ps)) {
    df <- t(df)
  }
  
  # Extract metadata
  df_meta <- sample.data.frame(ps) 
  if (!("sampleID" %in% colnames(df_meta))) {
    df_meta <- df_meta %>%
      mutate(across(where(is.character), as.factor)) %>% 
      rownames_to_column("sampleID")
  } else {
    df_meta <- df_meta %>%
      mutate(across(where(is.character), as.factor))
  }
  
  # Fit PCA
  pca_fit <- prcomp(na.omit(df), center = TRUE, scale. = scale)
  
  # Extract eigenvalues (variance explained)
  var <- pca_fit %>%
    tidy(matrix = "eigenvalues")
  
  # Extract PC scores (samples)
  pcs_fit <- pca_fit %>%
    tidy(matrix = "pcs")
  
  # Extract loadings (variables/OTUs)
  loadings <- pca_fit %>%
    tidy(matrix = "rotation")
  
  # Prepare the data for variance explained bar plot
  perc_expl_var <- pcs_fit %>%
    filter(PC < 6) %>%
    ggplot(aes(x = PC, y = percent)) +
    geom_bar(stat = "identity")
  
  # Prepare the PCA data for plotting
  pca_df_variances <- pca_fit %>%
    augment() %>%
    dplyr::rename("sampleID" = 1) %>%
    left_join(df_meta, by = "sampleID") %>%
    rename_with(~str_replace(., ".fitted", ""))
  
  # Ensure no NA values in the batch
  pca_df_variances <- pca_df_variances %>% filter(!is.na(.data[[batch]]))
  
  # Make the PCA biplot (samples and loadings)
  ord_plot <- pca_df_variances %>%
    ggplot(aes(PC1, PC2, color = as.factor(.data[[batch]]))) +
    geom_point(alpha = 0.7, size = 3) +
    labs(x = paste("PC1", round(var$percent[1] * 100, 2), "%"),
         y = paste("PC2", round(var$percent[2] * 100, 2), "%")) +
    scale_color_manual(values = color_palette) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = NA, colour = "white"),
          panel.border = element_blank(),
          axis.line = element_line()) +
    ggtitle(title) +
    guides(color = guide_legend(ncol = 2))
  
  # Add the loadings (arrows) to the plot
  loadings_plot <- loadings %>%
    filter(PC %in% c(1, 2)) %>%  # Only PC1 and PC2
    mutate(PC1 = 1 * 5,  # Adjust the scale of loadings (optional)
           PC2 = 2 * 5) %>%
    ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "blue", alpha = 0.6) +
    geom_text(aes(x = PC1, y = PC2, label = column), vjust = 1.5, color = "blue")
  
  # Combine both the PCA plot and the loadings arrows
  final_plot <- ord_plot + loadings_plot
  
  return(final_plot)
}


#PCA for continuous values
run_PCA_cont <- function(ps = ps,
                     batch = batch, 
                     scale = FALSE, 
                     title = title){
  
  #df <- t(otu_table(ps)) #if taxa_are_rows =T
  df <- otu_table(ps)
  
  if (taxa_are_rows(ps)) {
    df <-t(df)
  }
  
  df_meta <- sample.data.frame(ps)
  
  if (!("sampleID" %in% colnames(df_meta))) {
    df_meta <- df_meta %>%
      mutate(across(where(is.character), as.factor))  %>% 
      rownames_to_column("sampleID")
  } else {
    df_meta <- df_meta %>%
      mutate(across(where(is.character), as.factor))
  }
  # Calculate the distance matrix (Euclidean distance) from the OTU table
  dist_mat <- dist(df, method = "euclidean")
  #dist_mat <- vegdist(df, method = "bray")
  # Fit PCA
  pca_fit <- prcomp(na.omit(dist_mat), center = TRUE, scale. =scale)
  
  # Tidy the results
  var <- pca_fit %>% 
    tidy(matrix = "eigenvalues")
  
  pcs_fit <- pca_fit %>% 
    tidy(matrix = "pcs")
  
  perc_expl_var <- pcs_fit %>%
    filter(PC < 6) %>%
    ggplot(aes(x = PC, y = percent)) +
    geom_bar(stat = "identity")
  
  pca_df_variances <- pca_fit %>%
    augment() %>%
    dplyr::rename("sampleID" = 1) %>%
    left_join(df_meta) %>%
    rename_with(~str_replace(.,".fitted",""))
  
  # Make the plots
  # num_colors <- length(unique(df_meta[[batch]]))
  # colours <- color_palette[1:num_colors]
  pca_df_variances <- pca_df_variances %>% filter(!is.na(.data[[batch]]))
  
  ord_plot <- pca_df_variances %>%
    ggplot(aes(PC1, PC2, color = .data[[batch]])) + 
    geom_point(alpha = 0.7, size = 3) + 
    #geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    #geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(x = paste("PC1", round(var$percent[1]*100, 2), "%"),
         y = paste("PC2", round(var$percent[2]*100, 2), "%")) +
    scale_color_gradient(low = "#2A81FA", high = "#FC006E") + theme_minimal() + theme(panel.background = element_rect(fill = NA,colour = "white"),  panel.border = element_blank(), axis.line = element_line()) + ggtitle(title)
  
  return(ord_plot)
}

#ordination plot - PCoA, NMDS etc
run_ordplot <-
  function(ps = ps,
           distance = distance, method = method, 
           color_palette = color_palette,
           batch = batch,
           title = title) {
    # Perform ordination
    ordination <- ordinate(ps, method = method, distance = distance)
    ord_df <- as.data.frame(ordination$vectors)
    rownames(ord_df) <- rownames(sample_data(ps))
    ord_df <- cbind(sample_data(ps), ord_df)
    
    # Calculate percentages of variance explained by the first and second axes
    var_explained <-round(ordination$values$Eigenvalues / sum(ordination$values$Eigenvalues) * 100,
                          2)
    var_explained_first_axis <- var_explained[1]
    var_explained_second_axis <- var_explained[2]
    
    # Create plot
    ord_df <- ord_df %>% filter(!is.na(.data[[batch]]))
    #num_colors <- length(unique(ord_df)[[batch]])
    #colss <- color_palette[1:num_colors]
    #
    ord_plot <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = .data[[batch]])) +
      geom_point(alpha = 0.7, size =3) +
      #geom_hline(yintercept=0, linetype="dashed", alpha = 0.3) +
      #geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
       labs(
        x = paste(method, "1 (", var_explained_first_axis, "%)"),
        y = paste(method, "2 (", var_explained_second_axis, "%)"),
        title = title
      )+ scale_color_manual(values = color_palette) +
      theme_minimal() +theme(panel.background = element_rect(fill = NA,colour = "white"),  panel.border = element_blank(), axis.line = element_line())
    
    return(ord_plot)
  }
#ordination for continuous values
run_ordplot_cont <-
  function(ps = ps,
           distance = distance, method = method, 
           color_palette = color_palette,
           batch = batch,
           title = title) {
    # Perform ordination
    ordination <- ordinate(ps, method = method, distance = distance)
    ord_df <- as.data.frame(ordination$vectors)
    rownames(ord_df) <- rownames(sample_data(ps))
    ord_df <- cbind(sample_data(ps), ord_df)
    
    # Calculate percentages of variance explained by the first and second axes
    var_explained <-round(ordination$values$Eigenvalues / sum(ordination$values$Eigenvalues) * 100,
                          2)
    var_explained_first_axis <- var_explained[1]
    var_explained_second_axis <- var_explained[2]
    
    # Create plot
    ord_df <- ord_df %>% filter(!is.na(.data[[batch]]))
    #num_colors <- length(unique(ord_df)[[batch]])
    #colss <- color_palette[1:num_colors]
    #
    ord_plot <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = .data[[batch]])) +
      geom_point(alpha = 0.7, size =3) +
      #geom_hline(yintercept=0, linetype="dashed", alpha = 0.3) +
      #geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(
        x = paste(method, "1 (", var_explained_first_axis, "%)"),
        y = paste(method, "2 (", var_explained_second_axis, "%)"),
        title = title
      )+ scale_color_gradient(low = "#2A81FA", high = "#FC006E") + theme_minimal() + theme(panel.background = element_rect(fill = NA,colour = "white"),  panel.border = element_blank(), axis.line = element_line())
    
    return(ord_plot)
  }

#PCA with shape variables
run_PCA_shape <- function(ps = ps,
                          batch = batch, shape = shape,
                          scale = FALSE, 
                          title = title,
                          color_palette = color_palette,
                          shape_palette = shape_palette){
  
  df <- otu_table(ps)
  
  if (taxa_are_rows(ps)) {
    df <- t(df)
  }
  
  df_meta <- sample.data.frame(ps) 
  if (!("sampleID" %in% colnames(df_meta))) {
    df_meta <- df_meta %>%
      mutate(across(where(is.character), as.factor))  %>% 
      rownames_to_column("sampleID")
  } else {
    df_meta <- df_meta %>%
      mutate(across(where(is.character), as.factor))
  }
  
  pca_fit <- prcomp(na.omit(df), center = TRUE, scale. = scale)
  
  var <- pca_fit %>% 
    tidy(matrix = "eigenvalues")
  
  pcs_fit <- pca_fit %>% 
    tidy(matrix = "pcs")
  
  perc_expl_var <- pcs_fit %>%
    filter(PC < 6) %>%
    ggplot(aes(x = PC, y = percent)) +
    geom_bar(stat = "identity")
  
  pca_df_variances <- pca_fit %>%
    augment() %>%
    dplyr::rename("sampleID" = 1) %>%
    left_join(df_meta) %>%
    rename_with(~str_replace(.,".fitted",""))
  
  pca_df_variances <- pca_df_variances %>% filter(!is.na(.data[[batch]]))
  
  ord_plot <- pca_df_variances %>%
    ggplot(aes(PC1, PC2, fill = as.factor(.data[[batch]]), color = as.factor(.data[[batch]]), shape = as.factor(.data[[shape]]))) + 
    geom_point(alpha = 0.7, size = 3) + 
    labs(x = paste("PC1", round(var$percent[1] * 100, 2), "%"),
         y = paste("PC2", round(var$percent[2] * 100, 2), "%")) +
    scale_fill_manual(values = color_palette) + scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    theme_minimal() + 
    theme(panel.background = element_rect(fill = NA, colour = "white"),  
          panel.border = element_blank(), 
          axis.line = element_line()) + 
    ggtitle(title) + 
    guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))
  
  return(ord_plot)
}



########Diversity Analyses #######

# The function removes self comparisons, and colors points in the same group red. You'll get an S4 error, but that is OK.
# p = your phyloseq object
# m = distance metric
# s = phyloseq sample data column with your sample names
# d = sample data for groups
#t = ptot title


plotDistances = function(p = p, m = m, s = s, d = d, t =t) {
  
  # calc distances
  wu = phyloseq::distance(p, m)
  wu.m = melt(as.matrix(wu))
  
  # remove self-comparisons
  wu.m = wu.m %>%
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor,as.character)
  
  # get sample data (S4 error OK and expected)
  sd = sample.data.frame(p) %>% rownames_to_column("sampleID") %>%
  select(s, d) %>%
    mutate_if(is.factor, as.character)
  
  # combined distances with sample data
  colnames(sd) = c("Var1", "Type1")
  wu.sd = left_join(wu.m, sd, by = "Var1")
  
  colnames(sd) = c("Var2", "Type2")
  wu.sd = left_join(wu.sd, sd, by = "Var2")
  
  # plot
  ggplot(wu.sd, aes(x = Type2, y = value)) +
    theme_bw() +
    geom_point() +
    geom_boxplot(aes(color = ifelse(Type1 == Type2, "#FC006E", "black"))) +
    scale_color_identity() + theme_bw() +
    facet_wrap(~ Type1, scales = "free_x") + labs(x = d, y = "Value", title = paste0("Distance Metric = ", t)) +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) 
}


# Asymptotic alpha diversity From this paper: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13682 
#' Functional wrapper to get the output of the wonderful iNext library in the format I use to pipe into ggplot2.
#' @export
get_asymptotic_alpha = function(species, verbose = TRUE){
  
  if(any(colSums(species) < 1)){
    if(verbose){print("reads seem to be relative abundance. This is not ideal. Applying transformation to address this.")}
    species <- apply(species, 2, un_cpm)
  }
  #create output df, use progress bar if available
  
  alpha_diversity <- if(requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pbapply(X = species,
                     MARGIN =  2,
                     FUN = div_applier) } else {
                       apply(X = species,
                             MARGIN =  2,
                             FUN = div_applier)
                     }
  alpha_diversity = data.frame(t(alpha_diversity))
  
  row.names(alpha_diversity)  <- colnames(species)
  
  colnames(alpha_diversity) = c("Chao1", "Simpson Index", "Shannon Entropy")
  
  return(alpha_diversity)
}

#' Functional wrapper to get the output of the wonderful iNext library in the format I use to pipe into ggplot2.
#'
div_applier <- function(x){
  out_div <- vector(mode = "numeric", length = 3)
  
  out_div[1] <- iNEXT::ChaoRichness(x)$Estimator
  out_div[2] <- iNEXT::ChaoSimpson( x)$Estimator
  out_div[3] <- iNEXT::ChaoShannon( x)$Estimator
  
  return(out_div)
}


#' Undo CPM transformation
#'
un_cpm <- function(x) 1E6 * x/sum(x)

#########STATS################
# Linda Function -----
# https://github.com/zhouhj1994/LinDA/blob/master/R/linda.R
#' @param otu.tab data frame or matrix representing observed OTU table. Row: taxa; column: samples.
#' @param meta data frame of covariates. The rows of \code{meta} correspond to the columns of \code{otu.tab}.
linda <- function(otu.tab, meta, formula, type = 'count',
                  adaptive = TRUE, imputation = FALSE, pseudo.cnt = 0.5, corr.cut = 0.1,
                  p.adj.method = 'BH', alpha = 0.05,
                  prev.cut = 0, lib.cut = 1, winsor.quan = NULL, n.cores = 1) {
  if(any(is.na(otu.tab))) {
    stop('The OTU table contains NAs! Please remove!\n')
  }
  allvars <- all.vars(as.formula(formula))
  Z <- as.data.frame(meta[, allvars])
  
  ## preprocessing
  keep.sam <- which(colSums(otu.tab) >= lib.cut & rowSums(is.na(Z)) == 0)
  Y <- otu.tab[, keep.sam]
  Z <- as.data.frame(Z[keep.sam, ])
  names(Z) <- allvars
  
  n <- ncol(Y)
  keep.tax <- which(rowSums(Y > 0) / n >= prev.cut)
  Y <- Y[keep.tax, ]
  m <- nrow(Y)
  
  ## some samples may have zero total counts after screening taxa
  if(any(colSums(Y) == 0)) {
    ind <- which(colSums(Y) > 0)
    Y <- Y[, ind]
    Z <- as.data.frame(Z[ind, ])
    names(Z) <- allvars
    keep.sam <- keep.sam[ind]
    n <- ncol(Y)
  }
  
  ## scaling numerical variables
  ind <- sapply(1 : ncol(Z), function(i) is.numeric(Z[, i]))
  Z[, ind] <- scale(Z[, ind])
  
  ## winsorization
  if(!is.null(winsor.quan)) {
    Y <- winsor.fun(Y, winsor.quan)
  }
  
  ##
  if(grepl('\\(', formula)) {
    random.effect <- TRUE
  } else {
    random.effect <- FALSE
  }
  
  if(is.null(rownames(otu.tab))) {
    taxa.name <- (1 : nrow(otu.tab))[keep.tax]
  } else {
    taxa.name <- rownames(otu.tab)[keep.tax]
  }
  if(is.null(rownames(meta))) {
    samp.name <- (1 : nrow(meta))[keep.sam]
  } else {
    samp.name <- rownames(meta)[keep.sam]
  }
  
  ## handling zeros
  if(type == 'count') {
    if(any(Y == 0)) {
      N <- colSums(Y)
      if(adaptive) {
        logN <- log(N)
        if(random.effect) {
          tmp <- lmer(as.formula(paste0('logN', formula)), Z)
        } else {
          tmp <- lm(as.formula(paste0('logN', formula)), Z)
        }
        corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
        if(any(corr.pval <= corr.cut)) {
          cat('Imputation approach is used.\n')
          imputation <- TRUE
        } else {
          cat('Pseudo-count approach is used.\n')
          imputation <- FALSE
        }
      }
      if(imputation) {
        N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
        N.mat[Y > 0] <- 0
        tmp <- N[max.col(N.mat)]
        Y <- Y + N.mat / tmp
      } else {
        N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
        N.mat[Y > 0] <- 0
        tmp <- N[max.col(N.mat)]
        Y <- Y + N.mat / tmp
      }
    }
  }
  
  if(type == 'proportion') {
    if(any(Y == 0)) {
      ## Half minimum approach
      Y <- t(apply(Y, 1, function (x) {
        x[x == 0] <- 0.5 * min(x[x != 0])
        return(x)
      }))
    }
  }
  
  ## CLR transformation
  logY <- log2(Y)
  W <- t(logY) - colMeans(logY)
  
  ## linear regression
  oldw <- getOption('warn')
  options(warn = -1)
  if(!random.effect) {
    suppressMessages(fit <- lm(as.formula(paste0('W', formula)), Z))
    res <- do.call(rbind, coef(summary(fit)))
    d <- ncol(model.matrix(fit))
    df <- rep(n - d, m)
    tmp <- vcov(fit)
    res.cov <- foreach(i = 1 : m) %do% {tmp[((i-1)*d+1) : (i*d), ((i-1)*d+1) : (i*d)]}
    res.cov <- do.call(rbind, res.cov)
    rownames(res.cov) <- rownames(res)
    colnames(res.cov) <- rownames(res)[1 : d]
  } else {
    fun <- function(i) {
      w <- W[, i]
      fit <- lmer(as.formula(paste0('w', formula)), Z)
      list(coef(summary(fit)), vcov(fit))
    }
    if(n.cores > 1) {
      tmp <- mclapply(c(1 : m), function(i) fun(i), mc.cores = n.cores)
    } else {
      suppressMessages(tmp <- foreach(i = 1 : m) %do% fun(i))
    }
    res <- do.call(rbind, lapply(tmp, `[[`, 1))
    res.cov <- do.call(rbind, lapply(tmp, `[[`, 2))
  }
  options(warn = oldw)
  
  res.intc <- res[which(rownames(res) == '(Intercept)'), ]
  rownames(res.intc) <- NULL
  options(warn = -1)
  suppressMessages(tmp <- mlv(sqrt(n) * res.intc[, 1],
                              method = 'meanshift', kernel = 'gaussian') / sqrt(n))
  options(warn = oldw)
  baseMean <- 2 ^ (res.intc[, 1] - tmp)
  baseMean <- baseMean / sum(baseMean) * 1e6
  
  output.fun <- function(x) {
    res.voi <- res[which(rownames(res) == x), ]
    rownames(res.voi) <- NULL
    
    if(random.effect) {
      df <- res.voi[, 3]
    }
    
    log2FoldChange <- res.voi[, 1]
    lfcSE <- res.voi[, 2]
    oldw <- getOption('warn')
    options(warn = -1)
    suppressMessages(bias <- mlv(sqrt(n) * log2FoldChange,
                                 method = 'meanshift', kernel = 'gaussian') / sqrt(n))
    options(warn = oldw)
    log2FoldChange <- log2FoldChange - bias
    stat <- log2FoldChange / lfcSE
    
    pvalue <- 2 * pt(-abs(stat), df)
    padj <- p.adjust(pvalue, method = p.adj.method)
    reject <- padj <= alpha
    output <- cbind.data.frame(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, reject, df)
    rownames(output) <- taxa.name
    return(list(bias = bias, output = output))
  }
  
  cov.fun <- function(x) {
    tmp <- (1 : ncol(res.cov))[-c(1, which(colnames(res.cov) == x))]
    covariance <- as.data.frame(as.matrix(res.cov[which(rownames(res.cov) == x), tmp]))
    rownames(covariance) <- taxa.name
    colnames(covariance) <- colnames(res.cov)[tmp]
    return(covariance)
  }
  
  variables <- unique(rownames(res))[-1]
  variables.n <- length(variables)
  bias <- rep(NA, variables.n)
  output <- list()
  if(variables.n == 1) {
    covariance <- NULL
  } else {
    covariance <- list()
  }
  for(i in 1 : variables.n) {
    tmp <- output.fun(variables[i])
    output[[i]] <- tmp[[2]]
    bias[i] <- tmp[[1]]
    if(variables.n > 1) {
      covariance[[i]] <- cov.fun(variables[i])
    }
  }
  names(output) <- variables
  if(variables.n > 1) {
    names(covariance) <- variables
  }
  
  rownames(Y) <- taxa.name
  colnames(Y) <- samp.name
  rownames(Z) <- samp.name
  return(list(variables = variables, bias = bias, output = output, covariance = covariance, otu.tab.use = Y, meta.use = Z, clr.table = W))
}

winsor.fun <- function(Y, quan) {
  N <- colSums(Y)
  P <- t(t(Y) / N)
  cut <- apply(P, 1, quantile, quan)
  Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
  ind <- P > Cut
  P[ind] <- Cut[ind]
  Y <- round(t(t(P) * N))
  return(Y)
}

perform_anosim <- function(physeq, group_var, distance_method) {
  # Extract the distance matrix
  dist_matrix <- phyloseq::distance(physeq, method = distance_method)
  
  # Extract the grouping variable
  grouping_var <- sample.data.frame(physeq)[[group_var]]
  
  # Perform ANOSIM
  anosim_result <- anosim(dist_matrix, grouping_var)
  
  # Print the results
  cat("ANOSIM Results:\n")
  print(anosim_result)
  
  # Summarize the results
  cat("\nSummary of ANOSIM Results:\n")
  summary(anosim_result)
  
  # Plot the ANOSIM result
  plot(anosim_result, main = paste("ANOSIM Plot for", group_var))
} #group_vae and distance_method in ""
perform_pairwise_anosim <- function(physeq, group_var,distance_method, padj, perm=999){
  dist_matrix <- phyloseq::distance(physeq, method = distance_method)
  # Extract the grouping variable
  grouping_var <- sample.data.frame(physeq)[[group_var]]
  #calculate pairwise anosim
  pairwise_results <- anosim.pairwise(dist_matrix, grouping_var,sim.method = distance_method,p.adjust.m = padj, perm = perm)
  cat("Pairwise Anosim Results:\n")
  print(pairwise_results)
}

perform_pairwise_adonis <- function(physeq, group_var, distance_method, padj){
  dist_matrix <- phyloseq::distance(physeq, method = distance_method)
  
  # Extract the grouping variable
  grouping_var <- sample.data.frame(physeq)[[group_var]]
  #calculate pairwise adonis
  pairwise_results <- pairwise.adonis(dist_matrix, grouping_var, p.adjust.m = padj)
  cat("Pairwise Permanova Results:\n")
  print(pairwise_results)
}
  
perform_pairwise_adonis2 <- function(physeq, group_var, distance_method = distance_method, strata = NULL, nperm = 999) {
  # Extract the distance matrix
  dist_matrix <- phyloseq::distance(physeq, method = distance_method)
  
  # Extract the grouping variable
  meta_df <- sample.data.frame(physeq)

  # Perform pairwise PERMANOVA using pairwise.adonis2
  pairwise_results <- pairwise.adonis2(
    x = as.formula(paste("dist_matrix ~", group_var)),  # Model formula
    data = meta_df,           # Data frame with the grouping variable
    strata = strata,          # Strata (optional)
    nperm = nperm             # Number of permutations
  )
  
  # Print the results
  cat("Pairwise PERMANOVA Results:\n")
  print(pairwise_results)
  
  # Return the results for further use
  return(pairwise_results)
}


runPermanovaAnalysis <- function(pseq, num_permutations = 999, seed, method) {
  set.seed(seed)
  all <- NULL  # Initialize all as NULL
  
  asv_table <- otu_table(pseq)
  metadata_table <- sample.data.frame(pseq)
  
  metadata_cols <- colnames(metadata_table)
  
  for (i in 1:ncol(metadata_table)) {
    meta_id <- metadata_cols[i]
    metadatai <- as.data.frame(cbind(rownames(metadata_table), metadata_table[, i]))
    metadatai <- metadatai[complete.cases(metadatai), ]
    colnames(metadatai) <- c("V1", "V2")
    samples_w_data <- metadatai$V1
    metadata_each <- metadatai$V2
    
    # Check if the data has at least two dimensions
    if (length(unique(metadata_each)) > 1 && nrow(asv_table) > 1) {
      OTU_subset_1 <- subset(asv_table, rownames(asv_table) %in% samples_w_data)
      OTU_subset_1 <- as.matrix(OTU_subset_1)
      adon <- adonis2(formula = OTU_subset_1 ~ metadata_each, permutations = num_permutations, method = method)
      pval <- adon[5][[1]][1]
      Fa <- adon[4][[1]][1]
      r2 <- adon[3][[1]][1]
      df <- adon[1][[1]][1]
      obs <- adon[1][[1]][3]
      adjr2 <- RsquareAdj(r2, obs, df)
      all <- rbind(all, c(meta_id, Fa, r2, adjr2, pval))
    }
  }
  
  # Check if all has data before setting column names
  if (!is.null(all) && nrow(all) > 0) {
    colnames(all) <- c("metadata", "F", "r2", "adjr2", "pvalue")
    qval <- p.adjust(all[, "pvalue"], method = "BH")
    all2 <- cbind(all, qval)
    
    # Return the results as a data frame
    return(all2)
  } else {
    # Return NULL if there's no data
    return(NULL)
  }
}

kruskal_bh_test <- function(x, y) {
  kruskal.test(x ~ y)$p.value
}

permanova.pseq <- function(pseq, formula, seed, permutations = 999){
  set.seed(seed)
  perm.tb <- c()
  # extract abundance table
  asv.tb <- pseq %>% otu_tibble(column.id = "asv") %>% column_to_rownames("asv")
  # extract meta data
  metadata <- pseq %>% sample_tibble(column.id = "sample_id") %>% column_to_rownames("sample_id")
  
  #Compute euclidean distance over CLR-transformed values (i.e. Aitchison distance).
  dis_ait = dist(t(asv.tb), method = "euclidean")
  
  #ADONIS test (tests homogeneity of dispersion among groups)
  W <- dis_ait
  perm.df <- vegan::adonis2(formula = as.formula(paste0('W', formula)),
                            data = metadata,
                            permutations = permutations)
  pval <- perm.df[5][[1]][1]
  Fa <- perm.df[4][[1]][1]
  r2 <- perm.df[3][[1]][1]
  df <- perm.df[1][[1]][1]
  obs <- perm.df[1][[1]][3]
  adjr2 <- RsquareAdj(r2, obs, df)
  perm.tb <- rbind(perm.tb, c(Fa, r2, adjr2, pval))
 
  return(perm.tb)
}


