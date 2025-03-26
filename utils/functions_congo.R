# title: Functions for CONGO pilot project 
# Author: Ricardo RPS
# Date: 2025

pcks <- c('patchwork', 'broom', 'dplyr', 'ggplot2', 'dtplyr', 'data.table', 'vegan', 'patchwork')

for(package in pcks){
  if(package %in% (.packages())){
      print(paste("the package", package ,"is loaded"))
     }else{
       pacman::p_load(package, character.only= TRUE) 
      }
}


# faster psmelt----
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



# Zero imputation and CLR transformation ----
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
clr_c <- function(counts, samples_are = "cols"){
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


# Asymptotic alpha diversity ----
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



#---------------------------------------------------------------------
# Principal component analysis (PCA) with density plots per component
#---------------------------------------------------------------------
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
  df_meta <- data_meta %>% mutate(across(where(is.character), as.factor))
  # n = n_distinct(df_meta[batch])
  
  # fit PCA
  pca_fit <- df %>% 
    prcomp(center = T, scale = scale) 
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
    rename("sample_id" = 1) %>%
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
    scale_color_manual(values = batch_color) +
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
    scale_fill_manual(values = batch_color) +
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
    scale_fill_manual(values = batch_color)
  
  # plot all together
  pTop + theme(legend.position = 'none') + guide_area() + pMain  + pRight + theme(legend.position = 'none')  + plot_layout(guides = 'collect') +
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4)) & theme(legend.box = "horizontal")
}

#------------------------
# Box and density plots
#------------------------

box_plot_batch <- function(data, asv_name, batch, title = NULL){
  
  box_plot <- data %>% 
    ggplot(aes(x = .data[[batch]], y = .data[[asv_name]], fill = .data[[batch]])) + 
    stat_boxplot(geom = "errorbar", width = 0.4, alpha = 0.4) + 
    geom_boxplot(alpha = 0.5) + 
    scale_fill_manual(values = batch_color) + 
    theme_bw() + 
    theme(# axis.text.x = element_text(angle = 0, hjust = 0.5), 
          panel.grid = element_blank(),
          axis.title.x = element_blank(), 
          # axis.text = element_text(size = 10),
          # axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") + 
    labs(fill = paste0("Batch - ", batch), y = 'Value', title = title) 
  
  dens_plot <- data %>% 
    ggplot(aes(x = .data[[asv_name]], fill = .data[[batch]])) +  
    geom_density(alpha = 0.5) + 
    scale_fill_manual(values = batch_color) + 
    labs(title = asv_name, y = "Density", x = 'Value', fill = paste0("Batch - ", batch)) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          panel.grid = element_blank())
  
  box_plot + dens_plot + plot_layout(guides = 'collect')
  
}

#-------------
# RLE plot
#-------------

RleMicroRna2 <- function (object, maintitle = NULL, batch = batch, xlab = NA,
                          legend = TRUE, cex.lab = 1.2, cex.xaxis = 1, 
                          cex.yaxis = 1, abline.lwd=0.5, legend.cex = 0.8,
                          xaxis.dist.ratio = 0.1, outcex = 1, title.cex = 1.3) 
{
  colorfill =  mixOmics::color.mixo(batch)
  nARR = dim(object)[2]
  nGEN = dim(object)[1]
  y = apply(object, 1, median)
  mva = matrix(nrow = nGEN, ncol = nARR)
  for (i in 1:nARR) {
    
    x = object[, i]
    mva[ ,i] = (x - y)
  }
  med = apply(mva, 2, median)
  MIN = min(mva, na.rm = TRUE)
  MAX = max(mva, na.rm = TRUE)
  par(las = 3)
  plot(med, xlim = c(0, nARR + 1), ylim = c(MIN, MAX), axes = FALSE, 
       xlab = xlab, ylab = 'Deviations',cex.lab = cex.lab)
  colnames(mva) = colnames(object)
  res = boxplot(data.frame(mva), outline = TRUE, add = TRUE, col = colorfill,
                xaxt = 'n', outcex = outcex, cex.axis = cex.yaxis) #outcex for outlier
  axis(1, cex.axis = cex.xaxis, at = 1:ncol(object), labels = NA)
  points(med, type = 'p', col = 'blue', cex = outcex)
  lines(med, type = 'l', col = 'blue', lty = 'dotted')
  title(main = maintitle, cex.main = title.cex)
  abline(0, 0, col = 'red', lwd = abline.lwd)
  par(las = 0)
  end_point = 0.5 + ncol(object)  # add degrees to the x axis
  box.max = max(max(res$stats), max(res$out))
  box.min = min(min(res$stats), min(res$out))
  box.range = box.max - box.min
  text(seq(1.2, end_point, by = 1), par("usr")[3] - xaxis.dist.ratio*box.range, 
       srt = 60, adj = 1, xpd = TRUE,
       labels = paste(colnames(object)), cex = cex.xaxis)
  if(legend == TRUE){
    legend('topright', legend = unique(batch), pch=15, col = unique(colorfill), cex = legend.cex)
  }
}



