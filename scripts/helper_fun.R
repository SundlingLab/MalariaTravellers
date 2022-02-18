# script for helper functions
library(tidyverse)
library(tidyr) # data reformatting
library(stringi)
library(gridExtra)
library(viridis)
#BiocManager::install(version='devel')
#BiocManager::install("MOFA2",version='devel')
#devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))

library(MOFA2)

#library(reticulate)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(readxl)
library("correlation")
#BiocManager::install(("ComplexHeatmap"))
library(ComplexHeatmap)
#BiocManager::install(("circlize"))
library(circlize)

OMICS_colors <- setNames(brewer.pal(3,"Dark2")[c(1,3)], c("Plasma proteins","Immune cells"))
OMICS_colors3 <- setNames(brewer.pal(3,"Dark2"), c("Plasma proteins","IgG class response","Immune cells"))
#ENDEMIC_colors <- setNames(brewer.pal(3,"Paired")[c(2,3)], c("previously_exposed","primary_infected"))
ENDEMIC_colors <- setNames(brewer.pal(3,"PuOr")[c(1,3)], c("previously_exposed","primary_infected"))
TIME_colors <- setNames(brewer.pal(6,"PiYG"), c("Acute","D10","M1","M3","M6","Y1"))

theme12 <- theme(axis.title.x = element_text(size = 12),
                 axis.text.x = element_text(size = 12),
                 axis.title.y = element_text(size = 12),
                 axis.text.y = element_text(size = 12),
                 text = element_text(size=12))
theme_figpanel <- theme(axis.title.x = element_text(size = 5),
                 axis.text.x = element_text(size = 5),
                 axis.title.y = element_text(size = 6),
                 axis.text.y = element_text(size = 5),
                 text = element_text(size=5),
                 legend.title = element_text(size=5), #change legend title font size
                 legend.text = element_text(size=5))

## ==== Lollipop plot ====
lollipop_plot <- function(data){
  plot <- data %>% 
    ggplot() +
    geom_segment(aes(x= fct_reorder(feature, value), xend=feature, y=0, yend=value), size=0.5, 
                 color=ifelse(data$view == "Plasma proteins","#1B9E77","#7570B3")) +
    geom_point(aes(x=feature, y=value), 
               
               color=ifelse(data$view == "Plasma proteins","#1B9E77","#7570B3"),
               size=0.8) +
    geom_hline(yintercept=00, size=0.5) +
    theme_classic(base_size = 6) +
    #scale_x_discrete(position = "top") +
    labs(x= NULL,
         y= "scaled weight") +
    coord_flip() + 
    theme(legend.position="none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.title.x = element_text(size = 6),
          axis.text.x = element_text(size = 5),
          axis.title.y = element_text(size = 6),
          axis.text.y = element_text(size = 5),
          text = element_text(size=6))
  return(plot)
}



#function that outputs mean, lower limit and upper limit of 95% CI
data_summary <- function(x) {
  m <- mean(x, na.rm=TRUE)
  sem <-sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))
  ymin<-m-1.96*sem
  ymax<-m+1.96*sem
  return(c(y=m,ymin=ymin,ymax=ymax))
}



## rank norm function 
rank_norm_fun <- function (x, FUN = qnorm, ties.method = "average", na.action) 
{
  if (missing(na.action)) {
    na.action <- get(getOption("na.action"))
  }
  if (!is.function(na.action)) {
    stop("'na.action' must be a function")
  }
  
  value <- na.action(x)
  
  FUN(rank(x, ties.method = ties.method)/(length(x) + 1))
}

## range01 function
range01 <- function(x){(x-min(x))/(max(x)-min(x))}



## ==== Linear mixed effect moddel code with mefisto model input data ====
model.input.lm.fun <- function(data, feat){
  
  plot.input.data <- data %>% 
    mutate(Time = as.factor(gsub(".*\\|","",sample)),
           ID = gsub("\\|.*","",sample)) %>% 
    filter(feature == feat) %>% 
    mutate(group = factor(group, levels = c("primary_infected","previously_exposed")))
  
  #plot.input.data %>% ggplot(aes(x=Time,y=value,color=group, group=group)) +  geom_smooth()
  
  fit2 <- lmer(value ~ Time * group + (1|ID), REML = F, data = plot.input.data)
  summary(fit2)
  anova(fit2)
  #Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
  stat.sum <- as.data.frame(coef(summary(fit2))) %>% 
    rownames_to_column() %>% 
    
    filter(grepl("grouppreviously_exposed",rowname)) %>% 
    mutate(Time = ifelse(grepl("TimeD10",rowname),"D10",
                         ifelse(grepl("TimeM1",rowname),"M1",
                                ifelse(grepl("TimeM3",rowname),"M3",
                                       ifelse(grepl("TimeM6",rowname),"M6",
                                              ifelse(grepl("TimeY1",rowname),"Y1","Acute")))))) %>% 
    mutate(pval = `Pr(>|t|)`) %>% 
    mutate(Significant = ifelse(pval<=0.001,"***",
                                ifelse(pval<=0.01,"**",
                                       ifelse(pval<=0.05,"*","ns")))) %>% 
    #mutate(Time = gsub(":.*","",rowname),
    #       Time = gsub("Time","",Time)) %>% 
    column_to_rownames("rowname")
  
  plot <- cat_plot(fit2, 
                   pred = Time, 
                   modx = group, 
                   main.title = unique(plot.input.data$feature),
                   geom = "line") +
    labs(y="fitted value",
         x=NULL) +
    geom_text(data = stat.sum,
              mapping = aes(x = Time, y = 10, label = Significant),
              size=8,
              inherit.aes = FALSE) +
    scale_color_manual(values = ENDEMIC_colors) +
    theme(axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          text = element_text(size=12),
          legend.position = "bottom")
  
  
  return(plot)
  
}

## ===== Single time point plot ====
single_timepoint_plot <- function(plot.data) {
  plot.data <- plot.data %>% 
    mutate(group = as.factor(factor(group, levels = c("primary_infected", "previously_exposed"))))
  
  stat.test <- plot.data %>% 
    group_by(feature) %>% 
    wilcox_test(value ~ group) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    add_xy_position(x = "group", #fun = "median_iqr", 
                    dodge = 0.8) 
  
  plot <- plot.data %>%
    ggplot(aes(x=group, y= value,fill=group))+#fill=as.factor(Time))) +
    geom_boxplot(outlier.size = 0.1,size=0.2) +
    geom_point(pch = 21,size=1, position = position_jitterdodge()) + 
    #facet_wrap(~feature, scales = "free",ncol=3 ) +
    #geom_boxplot(alpha=0.5) +
    #geom_jitter(width = 0.25) +
    scale_fill_manual(values= ENDEMIC_colors) +
    stat_pvalue_manual(stat.test,
                       size = 2,
                       bracket.size = .2,
                       label = "p.adj.signif",
                       #tip.length = 0.01,
                       hide.ns = T,
                       inherit.aes = FALSE) +
    jtools::theme_nice() +
    theme(legend.position='top',
          #axis.title.x = element_text(size = 12),
          axis.text.x = element_blank(),#element_text(size = 12),
          #axis.title.y = element_text(size = 12),
          #axis.text.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          #text = element_text(size=12)
          ) +
    labs(x = "Acute time point",
         title = "",
         caption = "Wilcox test with fdr correction")
  return(plot)
}

## ===== longitudinal time point plot ====
longitudinal_timepoint_plot <- function(plot.data) {
  
  plot <- plot.data %>% 
    mutate(group = as.factor(factor(group, levels = c("primary_infected", "previously_exposed")))) %>% 
    ggplot(aes(x = Time, y = value, fill = group ,color = group, group = group)) +
    #facet_wrap(~feature, ncol = 5) +
    
    jtools::theme_nice() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
    geom_smooth(method="loess", alpha=0.4,size=0.2) +
    scale_fill_manual(values=ENDEMIC_colors) +
    scale_color_manual(values=ENDEMIC_colors) +
    labs(x = "") +
    jtools::theme_nice() +
    theme(legend.position='top',
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          text = element_text(size=12))
  return(plot)
}



## === MOFA2 functions to modify heatmap input ====


.check_and_get_factors <- function(object, factors) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(factors)))
  if (is.numeric(factors)) {
    stopifnot(all(factors <= object@dimensions$K))
    factors_names(object)[factors] 
  } else {
    if (paste0(factors, collapse = "") == "all") { 
      factors_names(object)
    } else {
      stopifnot(all(factors %in% factors_names(object)))
      factors
    }
  }
  
}

.check_and_get_covariates <- function(object, covariates) {
  if (!.hasSlot(object, "covariates") || is.null(object@covariates))
    stop("No covariates found in object.")
  stopifnot(!any(duplicated(covariates)))
  if (is.numeric(covariates)) {
    stopifnot(all(covariates <= object@dimensions$C))
    covariates_names(object)[covariates] 
  } else {
    if (paste0(covariates, collapse = "") == "all") { 
      covariates_names(object)
    } else {
      stopifnot(all(covariates %in% covariates_names(object)))
      covariates
    }
  }
}

.check_and_get_views <- function(object, views, non_gaussian=TRUE) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(views)))
  if (is.numeric(views)) {
    stopifnot(all(views <= object@dimensions$M))
    views <- views_names(object)[views]
  } else {
    if (paste0(views, sep = "", collapse = "") == "all") { 
      views <- views_names(object)
    } else {
      stopifnot(all(views %in% views_names(object)))
    }
  }
  
  # Ignore non-gaussian views  
  if (isFALSE(non_gaussian)) {
    non_gaussian_views <- names(which(object@model_options$likelihoods!="gaussian"))
    views <- views[!views%in%non_gaussian_views]
  }
  
  return(views)
}


.check_and_get_groups <- function(object, groups) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(groups)))
  if (is.numeric(groups)) {
    stopifnot(all(groups <= object@dimensions$G))
    groups_names(object)[groups] 
  } else {
    if (paste0(groups, collapse = "") == "all") { 
      groups_names(object)
    } else {
      stopifnot(all(groups %in% groups_names(object)))
      groups
    }
  }
}


.check_and_get_samples <- function(object, samples) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(samples)))
  if (is.numeric(samples)) {
    stopifnot(all(samples <= sum(object@dimensions$N)))
    unlist(samples_names(object))[samples] 
  } else {
    if (paste0(samples, collapse = "") == "all") { 
      unlist(samples_names(object))
    } else {
      stopifnot(all(samples %in% unlist(samples_names(object))))
      samples
    }
  }
}



get_imputed_data <- function(object, views = "all", groups = "all", features = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (length(object@imputed_data)==0) stop("imputed data not found, did you run: 'object <- impute(object)'?")
  
  # Get views and groups
  views <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  
  # Get features
  if (is(features, "list")) {
    stopifnot(all(sapply(seq_len(length(features)), function(i) all(features[[i]] %in% features_names(object)[[views[i]]]))))
    stopifnot(length(features)==length(views))
    if (is.null(names(features))) names(features) <- views
  } else {
    if (paste0(features, collapse="") == "all") { 
      features <- features_names(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }
  
  # Fetch mean
  imputed_data <- lapply(object@imputed_data[views], function(x) x[groups] )
  imputed_data <- lapply(seq_len(length(imputed_data)), function(m) lapply(seq_len(length(imputed_data[[1]])), function(p) imputed_data[[m]][[p]][as.character(features[[m]]),,drop=FALSE]))
  imputed_data <- .name_views_and_groups(imputed_data, views, groups)
  
  # Add feature intercepts
  # tryCatch( {
  #
  #   if (add_intercept & length(object@intercepts[[1]])>0) {
  #     intercepts <- lapply(object@intercepts[views], function(x) x[groups])
  #     intercepts <- .name_views_and_groups(intercepts, views, groups)
  #
  #     for (m in names(imputed_data)) {
  #       for (g in names(imputed_data[[m]])) {
  #         imputed_data[[m]][[g]] <- imputed_data[[m]][[g]] + intercepts[[m]][[g]][as.character(features[[m]])]
  #       }
  #     }
  #   } }, error = function(e) { NULL })
  
  # Convert to long data frame
  if (isTRUE(as.data.frame)) {
    
    imputed_data <- lapply(views, function(m) { 
      lapply(groups, function(g) { 
        tmp <- reshape2::melt(imputed_data[[m]][[g]])
        colnames(tmp) <- c("feature", "sample", "value")
        tmp <- cbind(view = m, group = g, tmp)
        return(tmp) 
      })
    })
    imputed_data <- do.call(rbind, do.call(rbind, imputed_data))
    
    
    factor.cols <- c("view","group","feature","sample")
    imputed_data[factor.cols] <- lapply(imputed_data[factor.cols], factor)
  }
  return(imputed_data)
}


#object <- model
#factor="Factor1"
#view="Immune cells"
#groups = "all" 
#features = 20
#annotation_features = NULL
#annotation_samples = NULL
#transpose = FALSE
#imputed = TRUE
#denoise = TRUE 
#max.value = NULL
#min.value = NULL

plot_data_heatmap_ordered <- function(object, factor, view = 1, groups = "all", features = 50, 
                                      annotation_features = NULL, annotation_samples = NULL, transpose = FALSE, 
                                      imputed = FALSE, denoise = FALSE, max.value = NULL, min.value = NULL, ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(factor)==1)
  stopifnot(length(view)==1)
  
  # Define views, factors and groups
  groups <- .check_and_get_groups(object, groups)
  factor <- .check_and_get_factors(object, factor)
  view <- .check_and_get_views(object, view)
  
  # Get weights
  W <- do.call(rbind, get_weights(object, views=view, factors=factor, as.data.frame = FALSE))
  
  # NOTE: By default concatenate all the groups
  Z <- lapply(get_factors(object)[groups], function(z) as.matrix(z[,factor]))
  Z <- do.call(rbind, Z)[,1]
  Z <- Z[!is.na(Z)]
  
  
  # Get data
  if (isTRUE(denoise)) {
    data <- predict(object, views=view, groups=groups)[[1]]
  } else {
    if (isTRUE(imputed)) {
      data <- get_imputed_data(object, view, groups)[[1]]
    } else {
      data <- get_data(object, views=view, groups=groups)[[1]]
    }
  }
  
  # Concatenate groups
  if (is(data, "list")) {
    data <- do.call(cbind, data)
  }
  
  # Subset features
  if (is(features, "numeric")) {
    if (length(features)==1) {
      features <- rownames(W)[tail(order(abs(W)), n=features)]
    } else {
      features <- rownames(W)[order(-abs(W))[features]]
    }
    # Sort features according to the weights
    features <- names(W[features,])[order(W[features,])]
  } else if (is(features, "character")) {
    stopifnot(all(features %in% features_names(object)[[view]]))
  } else {
    stop("Features need to be either a numeric or character vector")
  }
  data <- data[features,]
  
  
  # Select respective samples
  data <- data[,names(Z)]
  
  # Ignore samples with full missing views
  data <- data[, apply(data, 2, function(x) !all(is.na(x)))]
  
  # By default, sort samples according to the factor values
  order_samples <- names(sort(Z, decreasing = TRUE))
  order_samples <- order_samples[order_samples %in% colnames(data)]
  data <- data[,order_samples]
  
  # Add sample annotations
  if (!is.null(annotation_samples)) {
    
    # Predefined data.frame
    if (is.data.frame(annotation_samples)) {
      message("'annotation_samples' provided as a data.frame, please make sure that the rownames match the sample names")
      if (any(!colnames(data)%in%rownames(annotation_samples))) {
        stop("There are rownames in annotation_samples that do not correspond to sample names in the model")
      }
      annotation_samples <- annotation_samples[colnames(data), , drop = FALSE]
      
      # Extract metadata from the sample metadata  
    } else if (is.character(annotation_samples)) {
      stopifnot(annotation_samples%in%colnames(object@samples_metadata))
      # tmp <- tibble::column_to_rownames(object@samples_metadata,"sample")[order_samples,,drop=F]
      tmp <- object@samples_metadata
      rownames(tmp) <- tmp$sample
      tmp$sample <- NULL
      tmp <- tmp[order_samples,,drop=FALSE]
      annotation_samples <- tmp[,annotation_samples, drop=FALSE]
      rownames(annotation_samples) <- rownames(tmp)
    } else {
      stop("Input format for 'annotation_samples' not recognised ")
    }
    
    # Convert character columns to factors
    foo <- sapply(annotation_samples, function(x) is.logical(x) || is.character(x))
    if (any(foo)) annotation_samples[,which(foo)] <- lapply(annotation_samples[,which(foo),drop=FALSE], as.factor)
  }
  
  
  # Add feature annotations
  if (!is.null(annotation_features)) {
    stop("'annotation_features' is currently not implemented")
  }
  
  # Transpose the data
  if (transpose) {
    data <- t(data)
    if (!is.null(annotation_samples)) {
      annotation_features <- annotation_samples
      annotation_samples <- NULL
    }
    if (!is.null(annotation_features)) {
      annotation_samples <- annotation_features
      annotation_features <- NULL
    }
  }
  
  # Cap values
  if (!is.null(max.value)) data[data>=max.value] <- max.value
  if (!is.null(min.value)) data[data<=min.value] <- min.value
  return(data)
  dim(data)
  # Plot heatmap
  #pheatmap(data)#, 
  #annotation_row = annotation_features, 
  #annotation_col = annotation_samples, 
  #...
  #)
  
}


## for cow plot title
## https://stackoverflow.com/questions/50973713/ggplot2-creating-themed-title-subtitle-with-cowplot
draw_label_theme <- function(label, theme = NULL, element = "text", ...) {
  if (is.null(theme)) {
    theme <- ggplot2::theme_get()
  }
  if (!element %in% names(theme)) {
    stop("Element must be a valid ggplot theme element name")
  }
  
  elements <- ggplot2::calc_element(element, theme)
  
  cowplot::draw_label(label, 
                      fontfamily = elements$family,
                      fontface = elements$face,
                      colour = elements$color,
                      size = elements$size,
                      ...
  )
}

