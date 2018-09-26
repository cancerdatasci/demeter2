library(plyr);
library(magrittr);
library(reshape2)
library(tibble)
library(dplyr)
library(taigr);
library(jsonlite)
library(ggplot2)
library(dependr)
library(readr)
library(stringr)
library(pheatmap)
library(jmmBase)

load_demeter2_model <- function(model_dir, use_bayes = FALSE) {
    hp_df <- read.csv(file.path(model_dir, 'hp_data.csv'), stringsAsFactors = F, check.names = F) %>%
        column_to_rownames(var = 'hp')
    CL_df <- read.csv(file.path(model_dir, 'CL_data.csv'), row.names=1, stringsAsFactors = F, check.names = F)
    other_data <- read_json(file.path(model_dir, 'other_info.json'))
    results <- list(hairpin_data = hp_df,
                    CL_data = CL_df,
                    other_dat = other_data)
    if (file.exists(file.path(model_dir, 'Geff.csv'))) {
        results[['Geff']] <- read.csv(file.path(model_dir, 'Geff.csv'), row.names=NULL, stringsAsFactors = F)
        # results[['Seff']] <- read.csv(file.path(model_dir, 'Seff.csv'), row.names=NULL, stringsAsFactors = F)
    }
    if (use_bayes) {
        gene_means <- jmmBase::load_matrix_fast(file.path(model_dir, 'gene_means.csv'))
        gene_SDs <- jmmBase::load_matrix_fast(file.path(model_dir, 'gene_SDs.csv'))
        results %<>% append(list(
            gene_scores = t(gene_means),
            gene_SDs = t(gene_SDs)
        ))
    } else {
        gene_df <- jmmBase::load_matrix_fast(file.path(model_dir, 'gene_data.csv'))
        results %<>% append(list(
            gene_scores = t(gene_df)
        ))
    }
    return(results)
}


get_feature_dep_corrs <- function(df, feature_name, gene_controls = NULL, use_weights = TRUE) {
    res <- ddply(df, .(Gene), function(df) {
        df <- df[!is.na(df[[feature_name]]),]
        if (length(unique(df$D2)) < 2 || length(unique(df[[feature_name]])) < 2) {
            D2_cor <- NA
            D1_cor <- NA
            ATARiS_cor <- NA
            GA_cor <- NA
            # lm_cor <- NA
            CRISPR_cor <- NA
        } else {
            if (use_weights) {
                D2_cor <- wtd.cors(df$D2, df[[feature_name]], df$w)[1]
            } else {
                D2_cor <- cor(df$D2, df[[feature_name]], use = 'pairwise.complete.obs')
            }
            # lm_cor <- cor(df$lm, df[[feature_name]], use = 'pairwise.complete.obs')
            ATARiS_cor <- cor(df$ATARiS, df[[feature_name]], use = 'pairwise.complete.obs')
            GA_cor <- cor(df$GA, df[[feature_name]], use = 'pairwise.complete.obs')
            CRISPR_cor <- cor(df$CRISPR, df[[feature_name]], use = 'pairwise.complete.obs')
            D1_cor <- cor(df$D1, df[[feature_name]], use = 'pairwise.complete.obs')
        }
        return(data.frame(D1 = D1_cor, D2 = D2_cor, GA = GA_cor, ATARiS = ATARiS_cor, CRISPR = CRISPR_cor))
    })
    if (!is.null(gene_controls)) {
        res %<>% mutate(gene_type = 'none',
                        gene_type = ifelse(Gene %in% gene_controls[['pos']], 'ess', gene_type),
                        gene_type = ifelse(Gene %in% gene_controls[['neg']], 'non_ess', gene_type))
    }
    return(res)
}

get_feature_dep_corrs_sig <- function(df, feature_name, gene_controls = NULL) {
    get_cor_test <- function(df, dep_name, feature_name, min_samps = 5, weight_name = NULL) {
        if (sum(!is.na(df[[dep_name]])) >= min_samps) {
            if (!is.null(weight_name)) {
                if (length(unique(df[[dep_name]])) < 2 | length(unique(df[[feature_name]])) < 2) {
                    return(list(estimate = NA, p.value=NA))
                }
                wc <- wtd.cor(df[[dep_name]], df[[feature_name]], df[[weight_name]])
                c <- list(estimate = wc[1], p.value = wc[4])
            } else {
                c <- cor.test(df[[dep_name]], df[[feature_name]], use = 'pairwise.complete.obs')
            }
        } else {
            c <- list(estimate=NA, p.value=NA)
        }
        return(c)
    }
    df <- df[!is.na(df[[feature_name]]),]
    res <- ddply(df, .(Gene), function(df) {
        D2_cor <- get_cor_test(df, 'D2', feature_name, weight_name = 'w')
        D1_cor <- get_cor_test(df, 'D1', feature_name)
        # lm_cor <- get_cor_test(df, 'lm', feature_name)
        ATARiS_cor <- get_cor_test(df, 'ATARiS', feature_name)
        GA_cor <- get_cor_test(df, 'GA', feature_name)
        CRISPR_cor <- get_cor_test(df, 'CRISPR', feature_name)
        return(data.frame(D1 = D1_cor$estimate,
                          D2 = D2_cor$estimate,
                          # lm = lm_cor$estimate,
                          ATARiS = ATARiS_cor$estimate,
                          GA = GA_cor$estimate,
                          CRISPR = CRISPR_cor$estimate,
                          D1_p = D1_cor$p.value,
                          D2_p = D2_cor$p.value,
                          # lm_p = lm_cor$p.value,
                          ATARiS_p = ATARiS_cor$p.value,
                          GA_p = GA_cor$p.value,
                          CRISPR_p = CRISPR_cor$p.value))
    })
    if (!is.null(gene_controls)) {
        res %<>% mutate(gene_type = 'none',
                        gene_type = ifelse(Gene %in% gene_controls[['pos']], 'ess', gene_type),
                        gene_type = ifelse(Gene %in% gene_controls[['neg']], 'non_ess', gene_type))
    }
    return(res)
}

get_matrix_vector_corrs <- function(dep_data, target_dsets, vector, use_bayes = TRUE) {
    #Get correlations between each gene and a target vector for a set of gene score dsets
    ldply(target_dsets, function(dset) {
        use_CLs <- intersect(rownames(dep_data[[dset]]$gene_scores), names(vector))
        if (dep_data[[dset]]$mod_type == 'D2' & use_bayes) {
            cmat <- ldply(colnames(dep_data[[dset]]$gene_scores), function(gene) {
                r <- wtd.cors(dep_data[[dset]]$gene_scores[use_CLs, gene], vector[use_CLs],
                              weight = 1/dep_data[[dset]]$gene_score_SDs[use_CLs, gene]^2)
                data.frame(cor = r[1],
                           Gene = gene,
                           model = dset)
            })
        } else {
            cmat <- data.frame(cor = cor(dep_data[[dset]]$gene_scores[use_CLs,], vector[use_CLs], use = 'pairwise.complete.obs')[,1],
                               Gene = colnames(dep_data[[dset]]$gene_scores),
                               model = dset)
        }
    })
}


#pulled from https://github.com/mlist/R/blob/master/SSMD.R and modified
ssmdrobust <- function(pop1, pop2)
{
    beta <- median(pop1, na.rm=T) - median(pop2, na.rm=T)
    beta <- -beta / sqrt((mad(pop1, na.rm=T))^2 + (mad(pop2, na.rm=T))^2)
    return (beta)
}
ssmdMM <- function(pop1, pop2)
{
    beta <- mean(pop1, na.rm=T) - mean(pop2, na.rm=T)
    beta <- beta / sqrt((sd(pop1, na.rm=T))^2 + (sd(pop2, na.rm=T))^2)
    return (beta)
}
est_ssmd <- function(scores, calls, weights = NULL)
{
    if (is.null(weights)) {
        m1 <- mean(scores[calls], na.rm=T)
        m2 <- mean(scores[!calls], na.rm=T)
        v1 <- var(scores[calls], na.rm=T)
        v2 <- var(scores[!calls], na.rm=T)
    } else {
        m1 <- wtd.mean(scores[calls], weights[calls], na.rm=T)
        m2 <- wtd.mean(scores[!calls], weights[!calls], na.rm=T)
        v1 <- wtd.var(scores[calls], weights[calls], na.rm=T)
        v2 <- wtd.var(scores[!calls], weights[!calls], na.rm=T)
    }
    return((m2 - m1) / sqrt(v1 + v2))
}

get_pred_auc <- function(pred, res) {
    if (all(is.na(pred)) | all(is.na(res))) {
        return(list(performance_df = list(auc_roc = NA, auc_pr = NA)))
    } else {
        return(cdsr::evaluate_prediction_auc(pred, res))
    }
}

# get_GS_SSMD <- function(gene_scores, gene_type) {
#     if (all(is.na(gene_scores))) {
#         return(data.frame(rob_SSMD_non_ess = NA, rob_SSMD_unexp = NA, SSMD_non_ess = NA, SSMD_unexp = NA))
#     }
#     rob_SSMD_non_ess <- ssmdrobust(gene_scores[gene_type == 'non_ess'],
#                                    gene_scores[gene_type == 'ess'])
#     SSMD_non_ess <- ssmdMM(gene_scores[gene_type == 'non_ess'],
#                            gene_scores[gene_type == 'ess'])
#     if ('unexp' %in% unique(gene_type)) {
#         rob_SSMD_unexp <- ssmdrobust(gene_scores[gene_type == 'unexp'],
#                                        gene_scores[gene_type == 'ess'])
#         SSMD_unexp <- ssmdMM(gene_scores[gene_type == 'unexp'],
#                                      gene_scores[gene_type == 'ess'])
#     } else {
#         rob_SSMD_unexp <- NA
#         SSMD_unexp <- NA
#     }
#     return(data.frame(rob_SSMD_non_ess = rob_SSMD_non_ess, rob_SSMD_unexp = rob_SSMD_unexp,
#                 SSMD_non_ess = SSMD_non_ess, SSMD_unexp = SSMD_unexp))
# }


get_dep_feat_cor_ <- function(dep_dset, targ_gene, feature, use_bayes) {
  u_cls <- intersect(rownames(dep_data[[dep_dset]]$gene_scores), names(feature))
  if (dep_data[[dep_dset]]$mod_type == 'D2' & use_bayes) {
    wtd.cors(dep_data[[dep_dset]]$gene_scores[u_cls, targ_gene],
             feature[u_cls], weight = 1/dep_data[[dep_dset]]$gene_score_SDs[u_cls, targ_gene]^2)[,1]
  } else {
    cor(dep_data[[dep_dset]]$gene_scores[u_cls, targ_gene], feature[u_cls],
        use = 'pairwise.complete.obs')
  }
}

get_dep_feat_cors_ <- function(dep_dset, targ_gene, features, feat_gene,use_bayes) {
  u_cls <- intersect(rownames(dep_data[[dep_dset]]$gene_scores), rownames(features))
  if (dep_data[[dep_dset]]$mod_type == 'D2' & use_bayes) {
    feat_cors <- wtd.cors(dep_data[[dep_dset]]$gene_scores[u_cls, targ_gene],
             features[u_cls,], weight = 1/dep_data[[dep_dset]]$gene_score_SDs[u_cls, targ_gene]^2)[1,]
  } else {
    feat_cors <- cor(dep_data[[dep_dset]]$gene_scores[u_cls, targ_gene], features[u_cls,],
        use = 'pairwise.complete.obs')[1,]
  }
  feat_cors_z <- (feat_cors - mean(feat_cors, na.rm=T))/sd(feat_cors, na.rm=T)
  return(feat_cors_z[feat_gene])
}

#get stronger correlation (magnitude) between dependency and either the gene's own CN or GE
get_best_assoc <- function(dep_dset, targ_gene, feat_gene, z_score = FALSE, use_abs = TRUE) {
  if (!(targ_gene %in% colnames(dep_data[[dep_dset]]$gene_scores))) {
    return(NA)
  }
  if (feat_gene %in% colnames(feature_data$GE)) {
    if (z_score) {
      GE_cor <- get_dep_feat_cors_(dep_dset, targ_gene,
                                  feature_data$GE,
                                  feat_gene,
                                  use_bayes = use_bayes)
    } else {
      GE_cor <- get_dep_feat_cor_(dep_dset, targ_gene,
                                  feature_data$GE[, feat_gene],
                                  use_bayes = use_bayes)
    }
   } else {
    GE_cor <- NA
  }
  if (feat_gene %in% colnames(feature_data$CN)) {
    if (z_score) {
      CN_cor <- get_dep_feat_cors_(dep_dset, targ_gene,
                                  feature_data$CN,
                                  feat_gene,
                                  use_bayes = use_bayes)
    } else {
      CN_cor <- get_dep_feat_cor_(dep_dset, targ_gene,
                                  feature_data$CN[, feat_gene],
                                  use_bayes = use_bayes)
    }
  } else {
    CN_cor <- NA
  }
  if (use_abs) {
    return(pmax(abs(GE_cor), abs(CN_cor), na.rm=T))
  } else {
    if (is.na(GE_cor) & is.na(CN_cor)) {return(NA)}
    if (is.na(GE_cor)) {GE_cor <- 0}
    if (is.na(CN_cor)) {CN_cor <- 0}
    if (abs(GE_cor) > abs(CN_cor)) {
      return(GE_cor)
    } else {
      return(CN_cor)
    }
  }
}


#get correlation between dependency and feature
get_spec_assoc <- function(dep_dset, targ_gene, feat_gene, feat_type, use_bayes = TRUE) {
  if (!(targ_gene %in% colnames(dep_data[[dep_dset]]$gene_scores))) {
    return(NA)
  }
  if (!(feat_gene %in% colnames(feature_data[[feat_type]]))) {
    return(NA)
  }
  return(get_dep_feat_cor_(dep_dset, targ_gene,
                               feature_data[[feat_type]][, feat_gene],
                               use_bayes = use_bayes)
  )
}



renorm_gene_scores <- function(mod) {
    guide_eff <- mod$guide_Geff %>% set_names(mod$data_names$hps)
    gene_max_Geff <- rep(NA, length(mod$mapGH))
    gene_avg_Geff <- rep(NA, length(mod$mapGH))
    gene_avg_Seff <- rep(NA, length(mod$mapGH))
    gene_avg_hpoff <- rep(NA, length(mod$mapGH))
    gene_nGuides <- rep(NA, length(mod$mapGH))

    for (ii in seq(length(mod$mapGH))) {
        gene_max_Geff[ii] <- max(guide_eff[mod$mapGH[[ii]]])
        gene_avg_Geff[ii] <- mean(guide_eff[mod$mapGH[[ii]]])
        gene_avg_Seff[ii] <- mean(mod$guide_Seff[mod$mapGH[[ii]]])
        gene_nGuides[ii] <- length(mod$mapGH[[ii]])
        if ((gene_max_Geff[ii] == 0) | (is.null(mod$mapGH[[ii]]))) {
            mod$gene_score[, ii] <- 0
        } else {
            mod$gene_score[, ii] <- mod$gene_score[, ii]  * gene_max_Geff[ii]
            mod$q_gene_score_mean[, ii] <- mod$q_gene_score_mean[, ii] * gene_max_Geff[ii]
            # mod$guide_Geff[mod$mapGH[[ii]]] <- mod$guide_Geff[mod$mapGH[[ii]]] / gene_max_Geff[ii]
            gene_avg_hpoff[ii] <- mean(mod$hairpin_offset[mod$mapGH[[ii]]])
        }
    }
    mod$guide_Geff[mod$guide_Geff > 1] <-  1
    return(list(mod = mod,
                gene_stats = data.frame(max_Geff = gene_max_Geff,
                                        avg_Geff = gene_avg_Geff,
                                        avg_Seff = gene_avg_Seff,
                                        avg_hpoff = gene_avg_hpoff,
                                        nGuides = gene_nGuides)))
}

pheatmap_with_NA <- function(data, ...) {
    d_row <- dist(data)
    d_row[is.na(d_row)] <- mean(d_row, na.rm=T)
    row_clust <- hclust(d_row)
    d_col <- dist(t(data))
    d_col[is.na(d_col)] <- mean(d_col, na.rm=T)
    col_clust <- hclust(d_col)
    pheatmap(data, cluster_rows = row_clust, cluster_cols = col_clust, ...)
}

rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

make_scatter_compare <- function(x, y, xname, yname, gene_type = NULL, use_abs = FALSE) {
    c <- cor.test(x,y)
    df <- data.frame(x = x, y = y)
    if (!is.null(gene_type)) {df %<>% mutate(gene_type = gene_type)}
    g1 <- ggplot(df, aes(x, y)) +
        geom_point(alpha = 0.25) +
        # geom_density_2d() +
        geom_abline(slope = 1) +
        xlab(xname) + ylab(yname) +
        ggtitle(sprintf('Pearson Corr: %.3f, p: %.3g', c$estimate, c$p.value))
    if (!is.null(gene_type)) {
        g1 <- g1 + geom_density_2d(data = filter(df, gene_type %in% c('ess', 'non_ess')), aes(color = gene_type))
    }
    if (use_abs) {
        x = abs(x)
        y = abs(y)
    }
    t_stats <- t.test(y - x)
    g2 <- ggplot(df, aes(y - x)) + geom_vline(xintercept = 0, linetype = 'dashed')
    if (!is.null(gene_type)) {
        g2 <- g2 + geom_density(aes(color = gene_type))
    } else {
        g2 <- g2 + geom_density()
    }
    if (use_abs) {
        g2 <- g2 + xlab(paste0('|', yname, '| - |', xname, '|')) +
            ggtitle(sprintf('Avg magnitude diff: %.3f, p: %.3g', t_stats$estimate, t_stats$p.value))
    } else {
        g2 <- g2 + xlab(paste0(yname, ' - ', xname)) +
            ggtitle(sprintf('Avg diff: %.3f, p: %.3g', t_stats$estimate, t_stats$p.value))
    }
    g3 <- ggplot(df, aes(x)) + ggtitle(xname) + geom_vline(xintercept = 0, linetype = 'dashed')
    g4 <- ggplot(df, aes(y)) + ggtitle(yname) + geom_vline(xintercept = 0, linetype = 'dashed')
    if (!is.null(gene_type)) {
        g3 <- g3 + geom_density(aes(color = gene_type))
        g4 <- g4 + geom_density(aes(color = gene_type))
    } else {
        g3 <- g3 + geom_density()
        g4 <- g4 + geom_density()
    }
    xl <- layer_scales(g3)$x$range$range
    g4 <- g4 + xlim(xl)
    plot_grid(g1, g3, g2, g4, ncol = 2)
}

weight.cov <- function(x1, x2, w = NULL) {
    if (is.null(w)) {w <- rep(1, length(x1))}
    w <- w / sum(w, na.rm=T)
    x1 <- x1 - weighted.mean(x1, w, na.rm=T)
    x2 <- x2 - weighted.mean(x2, w, na.rm=T)
    N <- sum(!is.na(x1) & !is.na(x2))
    return(sum(w * x1 * x2, na.rm=T) * N/(N-1))
}

# weight.cor <- function(x1, x2, w = NULL) {
#     if (is.null(w)) {w <- rep(1, length(x1))}
#     w <- w / sum(w, na.rm=T)
#     C <- weight.cov(x1, x2, w)
#     C1 <- weight.cov(x1, x1, w)
#     C2 <- weight.cov(x2, x2, w)
#     return(C / sqrt(C1 * C2))
# }


plot_colorByDensity = function(df, var1, var2, size = 1) {
    #modified from http://knowledge-forlife.com/r-color-scatterplot-points-density/
    x <- densCols(df[[var1]],df[[var2]], colramp=colorRampPalette(c("black", "white")))
    df$dens <- col2rgb(x)[1,] + 1L
    ggplot(df, aes_string(var1, var2, 'color' = 'dens')) +
        geom_point(size = size) +
        guides(color = guide_colorbar(title = 'Density', ticks = FALSE, label = FALSE)) +
        scale_color_distiller(palette = "Spectral")
}

unpack_gene_families <- function(gene_scores) {
    #unpack data from gene families
    gene_families <- which(grepl('&', colnames(gene_scores)))
    for (id in gene_families) {
        genes <- str_split(colnames(gene_scores)[id], '&')[[1]]
        a <- matrix(gene_scores[, id], ncol = length(genes), nrow = nrow(gene_scores), byrow=F,
                    dimnames = list(rownames(gene_scores), genes))
        gene_scores %<>% cbind(a)
    }
    return(gene_scores[, !grepl('&', colnames(gene_scores))])
}

read_ATARiS_data <- function(GS_file) {
    GS <-read.gct(GS_file)
    GS <- GS[!grepl('NO_CURRENT', rownames(GS)), ] #get rid of scores for non-genes
    GS <- GS[grepl('[[:alnum:]+]_1_', rownames(GS)), ] #take first gene score if multiple
    rownames(GS) <- str_match(rownames(GS), '([:alnum:]+)_')[,2]
    return(t(GS))
}

read_gene_avgs <- function(rds_file) {
    GS <- read_rds(rds_file)
    GS <- GS[, !grepl('NO_CURRENT', colnames(GS))]
    GS <- GS[, colnames(GS) != 'NA']
    return(GS)
}


get_combined_data <- function(dep_data, feature_data, use_bayes) {
    dlist <- list()
    for (dname in names(dep_data)) {
        dlist[[dname]] <- dep_data[[dname]]$gene_scores
        if (!is.null(dep_data[[dname]]$gene_score_SDs)) {
            dlist[[paste0(dname, '_sd')]] <- dep_data[[dname]]$gene_score_SDs
        }
    }
    for (fname in names(feature_data)) {
        dlist[[fname]] <- feature_data[[fname]]
    }
    comb_GS <- make.a.tidy.dataset(dlist, na.rows = TRUE, na.cols = TRUE) %>%
        dplyr::rename(CCLE_ID = Gene, Gene = Sample) %>%
        filter(!grepl('no_current', Gene, ignore.case = TRUE),
               !grepl('control', Gene, ignore.case = TRUE))
    return(comb_GS)
}

get_pr <- function(df, targ_var, neg_con = 'non_ess') {
    if (neg_con == 'non_ess') {
        df %<>% filter(ov_gene_type %in% c('ess', 'non_ess'))
        if (sum(!is.na(df[[targ_var]])) < 5) {
            return(list(performance_df = list(auc_roc = NA, auc_pr = NA)))
        }
        return(cdsr::evaluate_prediction_auc(-df[[targ_var]], df[["ov_gene_type"]] == 'ess'))
    } else {
        df %<>% filter(gene_type %in% c('ess', 'unexp'))
        if (sum(!is.na(df[[targ_var]])) < 5) {
            return(list(performance_df = list(auc_roc = NA, auc_pr = NA)))
        }
        return(cdsr::evaluate_prediction_auc(-df[[targ_var]], df[["gene_type"]] == 'ess'))

    }
}

# make_CL_vs_avg_plot <- function(dname) {
#     test_CLs <- per_CL_SSMD[!is.na(per_CL_SSMD[[paste0(dname, '_SSMD_unexp')]]),] %>%
#         arrange(desc(ref)) %>% head(3) %>% .[['CCLE_ID']]
#     test_CLs %<>% append(per_CL_SSMD[!is.na(per_CL_SSMD[[paste0(dname, '_SSMD_unexp')]]),] %>% arrange((ref)) %>% head(3) %>% .[['CCLE_ID']])
#     used_data <- comb_GS %>%
#         filter(CCLE_ID %in% test_CLs) %>%
#         mutate(CCLE_ID = factor(CCLE_ID, levels = test_CLs))
#     ggplot(used_data, aes_string(paste0(dname, '_avg'), dname)) +
#         geom_point(alpha = 0.5, size = 1) +
#         geom_density_2d(data = filter(used_data, ov_gene_type %in% c('ess', 'non_ess')),
#                         aes(color = ov_gene_type)) +
#         geom_abline() +
#         geom_smooth() +
#         facet_wrap(~CCLE_ID)
# }
make_CL_vs_avg_plot <- function(dep_name, test_CLs, gene_avgs) {
    df <- dep_data[[dep_name]]$gene_scores[test_CLs,] %>%
        as.data.frame() %>%
        rownames_to_column(var = 'CCLE_ID') %>%
        melt(id.vars = 'CCLE_ID') %>%
        set_colnames(c('CCLE_ID', 'Gene', 'score')) %>%
        left_join(gene_avgs %>%
                      filter(dset == dep_name) %>%
                      .[, c('Gene', 'avg_score', 'gene_type')], by = 'Gene') %>%
        # mutate(CCLE_ID = factor(CCLE_ID, levels = test_CLs))
        mutate(CL_name = str_match(CCLE_ID, '([:alnum:]+)_')[,2]) %>% 
        mutate(gene_type = plyr::revalue(gene_type, c(non_essential = 'non-essential')))
    ggplot(df, aes(avg_score, score)) +
        geom_point(size = 0.5, color = 'grey') +
        geom_density_2d(data = filter(df, gene_type %in% c('essential', 'non-essential')),
                        aes(color = gene_type)) +
        geom_abline() +
        # geom_smooth(color = 'orange') +
        xlab('Cell line avg. dependency score') + ylab('Cell line dependency score') +
        guides(color = guide_legend(title = 'Gene type', override.aes = list(lwd = 3))) +
        ggtitle(dep_name) +
        facet_wrap(~CL_name)
}

make_gs_dist_plots <- function(dataset, gene_avgs) {
    #Make plots showing dists of per-gene averages
    if (dataset == 'Achilles') {
        cur_dset <- c('D2_Ach', 'GA_Ach', 'RSA_Ach')
    } else {
        cur_dset <- c('D2_DRIVE', 'GA_DRIVE', 'RSA_DRIVE')
    }
    r <- gene_avgs %>% filter(dset %in% cur_dset) %>% .[['avg_score']] %>% quantile(c(0.0001, 0.9999), na.rm=T)

    cur_gene_avgs <- filter(gene_avgs, dset %in% cur_dset) %>%
        mutate(model = str_match(dset, '([:alnum:]+)_.+')[,2])
    p <- ggplot(cur_gene_avgs, aes(avg_score, fill = gene_type)) +
        geom_density(alpha = 0.5) +
        xlim(r) +
        xlab('Cell line avg. dependency score') +
        ylab('Density') +
        guides(fill = guide_legend(title = 'Gene Type')) +
        # scale_fill_Publication() +
        scale_fill_manual(values = gene_type_pal) +
        facet_wrap(~model)
    return(p)
    # title <- ggdraw() + draw_label(dataset, fontface='bold', size = 20)
    # plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
}

get_feature_cor <- function(df, feature_name, dname) {
    if (dep_data[[dname]]$mod_type == 'D2' & use_bayes) {
        return(wtd.cors(df[[feature_name]], df[[dname]], 1/df[[paste0(dname, '_sd')]]^2))
    } else {
        return(cor(df[[feature_name]], df[[dname]], use = 'pairwise.complete.obs'))
    }
}
get_pairs <- function(vector) {
    pairs <- data.frame()
    for (ii in seq(length(vector) - 1)) {
        for (jj in seq(from = ii+1, to = length(vector))) {
            if (!is.na(vector[ii]) & !is.na(vector[jj])) {
                pairs %<>% rbind(data.frame(a = vector[ii], b = vector[jj]))
            }
        }
    }
    pairs %<>% mutate(a = as.vector(a), b = as.vector(b))
    return(pairs)
}

plot_with_cor <- function(df, var1, var2) {
    c <- cor.test(df[[var1]], df[[var2]], use = 'pairwise.complete.obs', method = 'spearman')
    ggplot(df, aes_string(var1, var2)) +
        geom_point(alpha = 0.5) +
        geom_smooth() +
        ggtitle(sprintf('Spearman r: %.3f, p: %.3f', c$estimate, c$p.value))
}

get_gene_type <- function(Gene, ess_genes, non_ess_genes) {
    gene_type <- 'other'
    gene_type <- ifelse(Gene %in% ess_genes, 'essential', gene_type)
    gene_type <- ifelse(Gene %in% non_ess_genes, 'non_essential', gene_type)
    return(gene_type)
}

get_mat_dep_cors <- function(dep_names, target_mat, use_bayes) {
    cors_df <- ldply(dep_names, function(cur_dset) {
        # print(cur_dset)
        cur_gene_scores <- dep_data[[cur_dset]]$gene_scores
        if (dep_data[[cur_dset]]$mod_type == 'D2' & use_bayes) {
            cur_gene_score_SDs <- dep_data[[cur_dset]]$gene_score_SDs
        }
        common_CLs <- intersect(rownames(cur_gene_scores), rownames(target_mat))
        common_genes <- intersect(colnames(cur_gene_scores), colnames(target_mat))

        ldply(common_genes %>% set_names(common_genes), function(gene) {
            if (dep_data[[cur_dset]]$mod_type == 'D2' & use_bayes) {
                return(wtd.cors(cur_gene_scores[common_CLs, gene],
                                target_mat[common_CLs, gene],
                                1/cur_gene_score_SDs[common_CLs, gene]^2))
            } else {
                return(cor(cur_gene_scores[common_CLs, gene],
                           target_mat[common_CLs, gene],
                           use = 'pairwise.complete.obs'))
            }
        }) %>% set_colnames(c('Gene', 'cor')) %>% mutate(dset = cur_dset)
    }) %>% mutate(gene_type = get_gene_type(Gene, hart_ess, hart_non_ess))
}


get_SSMD <- function(gene_scores, gene_calls, weights = NULL, robust = FALSE) {
    if (all(is.na(gene_scores))) {
        return(NA)
    }
    if (robust) {
        stopifnot(is.null(weights))
        return(ssmdrobust(gene_scores[gene_calls], gene_scores[!gene_calls]))
    } else {
        return(est_ssmd(gene_scores, gene_calls, weights))
    }
}

get_gene_calls_NE <- function(gene_type) {
    gene_calls <- rep(NA, length(gene_type))
    gene_calls[gene_type == 'essential'] <- TRUE
    gene_calls[gene_type == 'non_essential'] <- FALSE
    return(gene_calls)
}
get_gene_calls_UE <- function(gene_type, gene_exp, GE_thresh) {
    gene_calls <- rep(NA, length(gene_type))
    gene_calls[gene_type == 'essential'] <- TRUE
    gene_calls[gene_exp < GE_thresh] <- FALSE
    return(gene_calls)
}

clean_CL_names_with_missing <- function(CL_names) {
    new_names <- CleanCellLineName(CL_names)
    new_names[is.na(new_names)] <- CL_names[is.na(new_names)]
    return(new_names)
}

#load specified set of dependency dataset, and apply some initial preprocessing
load_all_dep_data <- function(dep_dnames, 
                              dep_datasets, 
                              include_gene_families = FALSE, 
                              use_bayes=TRUE) {
    dep_data <- lapply(dep_dnames, function(dep_name) {
        print(sprintf('Loading %s', dep_name))
        cur_dep_info <- dep_datasets[[dep_name]]
        if (cur_dep_info$mod_type == 'D2') {
            cur_dep_info$mod <- load_demeter2_model(cur_dep_info$path, use_bayes = use_bayes)
            cur_dep_info$gene_scores <- cur_dep_info$mod$gene_scores
            if (use_bayes) {
                cur_dep_info$gene_score_SDs <- cur_dep_info$mod$gene_SDs
            }
            if (include_gene_families) {
                cur_dep_info$gene_scores <- unpack_gene_families(cur_dep_info$gene_scores)
                if (use_bayes) {
                    cur_dep_info$gene_score_SDs <- unpack_gene_families(cur_dep_info$gene_score_SDs)
                }
            } else {
              #throw out gene families, which include gene names catted with &
                cur_dep_info$gene_scores <- cur_dep_info$gene_scores[, !grepl('&', colnames(cur_dep_info$gene_scores))]
                if (use_bayes) {
                    cur_dep_info$gene_score_SDs <- cur_dep_info$gene_score_SDs[, !grepl('&', colnames(cur_dep_info$gene_score_SDs))]
                }
            }
        }
        if (cur_dep_info$mod_type %in% c('D1', 'GA', 'RSA', 'mageck', 'ATARiS')) {
            #load this from taiga
            cur_dep_info$gene_scores <- load.from.taiga(
                data.name=cur_dep_info$taiga_name,
                data.version=cur_dep_info$taiga_version,
                data.file = cur_dep_info$taiga_file,
                transpose = T)
            if (grepl('\\(', colnames(cur_dep_info$gene_scores)[1])) {
                colnames(cur_dep_info$gene_scores) <- str_match(colnames(cur_dep_info$gene_scores), '\\((.+)\\)')[,2]
            }
            if (sum(duplicated(colnames(cur_dep_info$gene_scores))) > 0) {
              print('removing duplicated genes')
              cur_dep_info$gene_scores <- cur_dep_info$gene_scores[, !duplicated(colnames(cur_dep_info$gene_scores))]
            }
        }
        if (cur_dep_info$mod_type == 'CERES') {
            #load this from taiga
            cur_dep_info$gene_scores <- load.from.taiga(
                data.name=cur_dep_info$taiga_name,
                data.version=cur_dep_info$taiga_version,
                data.file = cur_dep_info$taiga_file,
                transpose = F)
            colnames(cur_dep_info$gene_scores) <-
                str_match(colnames(cur_dep_info$gene_scores), ' \\((.+)\\)')[,2]
        }
        # if (cur_dep_info$dataset == 'DRIVE' | cur_dep_info$dataset == 'Marcotte') {
        #     #map CL names to CCLE_ID, and only take those in CCLE
        #     rownames(cur_dep_info$gene_scores) <- clean_CL_names_with_missing(rownames(cur_dep_info$gene_scores))
        #     cur_dep_info$gene_scores <- cur_dep_info$gene_scores[!is.na(rownames(cur_dep_info$gene_scores)),]
        #     if (cur_dep_info$mod_type == 'D2' & use_bayes) {
        #         rownames(cur_dep_info$gene_score_SDs) <- clean_CL_names_with_missing(rownames(cur_dep_info$gene_score_SDs))
        #         cur_dep_info$gene_score_SDs <- cur_dep_info$gene_score_SDs[!is.na(rownames(cur_dep_info$gene_score_SDs)),]
        #     }
        # }
        
        #ignore genes which are NA for all cell lines
        used_genes <- colSums(!is.na(cur_dep_info$gene_scores)) > 0
        cur_dep_info$gene_scores <- cur_dep_info$gene_scores[, used_genes]
        if (cur_dep_info$mod_type == 'D2' & use_bayes) {
            cur_dep_info$gene_score_SDs <- cur_dep_info$gene_score_SDs[, used_genes]
        }
        # #ignore XLOC genes
        # cur_dep_info$gene_scores <- cur_dep_info$gene_scores[,!grepl('XLOC', colnames(cur_dep_info$gene_scores))]
        # if (cur_dep_info$mod_type == 'D2' & use_bayes) {
        #     cur_dep_info$gene_score_SDs <- cur_dep_info$gene_score_SDs[,!grepl('XLOC', colnames(cur_dep_info$gene_score_SDs))]
        # }
        # #ignore non Entrez genes
        # cur_dep_info$gene_scores <- cur_dep_info$gene_scores[,grepl('[0-9]', colnames(cur_dep_info$gene_scores))]
        # if (cur_dep_info$mod_type == 'D2' & use_bayes) {
        #   cur_dep_info$gene_score_SDs <- cur_dep_info$gene_score_SDs[,grepl('[0-9]', colnames(cur_dep_info$gene_score_SDs))]
        # }
        cur_dep_info$mod <- NULL
        return(cur_dep_info)
    })
    return(dep_data)
}

load_all_feature_data <- function(feature_dnames, feature_datasets, all_dep_genes = NULL, all_dep_CLs = NULL) {
    feature_data <- lapply(feature_dnames, function(feature_name) {
        print(sprintf('Loading %s', feature_name))
        cur_feature_info <- feature_datasets[[feature_name]]
        if (is.null(cur_feature_info$transpose)) {
            cur_feature_info$transpose <- FALSE
        }
        #load this from taiga
        cur_data <- load.from.taiga(
            data.name=cur_feature_info$taiga_name,
            data.version=cur_feature_info$taiga_version,
            data.file = cur_feature_info$data.file,
            transpose = cur_feature_info$transpose)
        if (feature_name == 'GE') {
            colnames(cur_data) <- str_match(colnames(cur_data), ' \\((.+)\\)')[, 2]
            cur_data <- cur_data[, !duplicated(colnames(cur_data))] #get rid of any duplicated columns
        }
        if (feature_name == 'CN') {
            colnames(cur_data) <- str_match(colnames(cur_data), ' \\((.+)\\)')[, 2]
        }
        if (feature_name %in% c('MUT_MIS', 'MUT_DAM', 'MUT_HOT')) {
            colnames(cur_data) <- str_match(colnames(cur_data), ' \\((.+)\\)')[, 2]
        }
        return(cur_data)
    })

    if ('MUT' %in% names(feature_datasets)) {
        #parse hot-spot misense and damaging mutation datasets
        feature_data$MUT_HOT <- feature_data$MUT %>%
            mutate(cell_line = Tumor_Sample_Barcode) %>%
            filter(
                Variant_Classification == 'Missense_Mutation',
                (isTCGAhotspot) | (isCOSMIChotspot)) %>%
            select(cell_line, Entrez_Gene_Id) %>%
            mutate(mutation = 1) %>%
            unique() %>%
            reshape2::acast(cell_line ~ Entrez_Gene_Id, value.var = 'mutation', fill = 0)
        feature_data$MUT_DAM <- feature_data$MUT %>%
            mutate(cell_line = Tumor_Sample_Barcode) %>%
            filter(isDeleterious) %>%
            select(cell_line, Entrez_Gene_Id) %>%
            mutate(mutation = 1) %>%
            unique() %>%
            reshape2::acast(cell_line ~ Entrez_Gene_Id, value.var = 'mutation', fill = 0)
        feature_data <- feature_data[setdiff(names(feature_data), 'MUT')]
    }

    #subset feature datasets to use only specified lists of genes and CLs
    for (feature_dname in names(feature_data)) {
        if (!is.null(all_dep_genes)) {
            feature_data[[feature_dname]] <- feature_data[[feature_dname]][, colnames(feature_data[[feature_dname]]) %in% all_dep_genes]
        }
        if (!is.null(all_dep_CLs)) {
            feature_data[[feature_dname]] <- feature_data[[feature_dname]][rownames(feature_data[[feature_dname]]) %in% all_dep_CLs, ]
        }
    }

    #align rownames across mutation datasets
    if (all(c('MUT_DAM', 'MUT_HOT') %in% names(feature_datasets))) {
        common_MUT_CLs <- rownames(feature_data$MUT_DAM)
        # feature_data$MUT_MIS <- feature_data$MUT_MIS[common_MUT_CLs,]
        feature_data$MUT_DAM <- feature_data$MUT_DAM[common_MUT_CLs,]
        feature_data$MUT_HOT <- feature_data$MUT_HOT[common_MUT_CLs,]
    }
    return(feature_data)
}


scatter_compare_GA <- function(dset1, dset2, overlapping_genes = FALSE) {
    if (overlapping_genes) {
        use_genes <- gene_avgs_wide %>%
            na.omit()
    } else {
        use_genes <- gene_avgs_wide
    }
    c <- cor.test(use_genes[[dset1]], use_genes[[dset2]], use = 'pairwise.complete.obs')
    ggplot(use_genes, aes_string(dset1, dset2)) +
        geom_point(data = use_genes %>% filter(gene_type == 'other'), aes(color = 'other'), alpha = 0.25, cex = 1, color = 'black') +
        geom_point(data = use_genes %>% filter(gene_type == 'non_essential'), aes(color = 'non_essential'), alpha = 0.5, color = 'blue') +
        geom_point(data = use_genes %>% filter(gene_type == 'essential'), aes(color = 'essential'), alpha = 0.5, color = 'red') +
        geom_abline(slope = 1) +
        ggtitle(sprintf('Pearson cor: %.3f', c$estimate))
}

# get_cors_internal <- function(ref, dep, SD, genes, cur_use_bayes, cur_dset) {
#     res <- ldply(genes, function(gene) {
#         if (cur_use_bayes) {
#             cur_cor <- wtd.cors(ref,
#                                 dep[, gene],
#                                 1/SD[, gene]^2)
#         } else {
#             cur_cor <- cor(ref,
#                            dep[, gene],
#                            use = 'pairwise.complete.obs')
#         }
#         cur_cor[is.infinite(cur_cor)] <- NA
#         best_ind <- which.max(abs(cur_cor))
#         cor_df <- data.frame(cor = cur_cor[, best_ind],
#                              feat_gene = colnames(ref)[best_ind],
#                              targ_gene = gene)
#         return(cor_df)
#     }, .parallel=TRUE, .paropts = list(.noexport = c('dep_data', 'feature_data'))) %>% mutate(dset = cur_dset)
#     return(res)
# }

# get_cors_internal <- function(cache_file) {
#     inputs <- read_rds(cache_file)
#     res <- ldply(inputs$genes, function(gene) {
#         if (inputs$cur_use_bayes) {
#             cur_cor <- wtd.cors(inputs$ref,
#                                 inputs$dep[, gene],
#                                 1/inputs$SD[, gene]^2)
#         } else {
#             cur_cor <- cor(inputs$ref,
#                            inputs$dep[, gene],
#                            use = 'pairwise.complete.obs')
#         }
#         cur_cor[is.infinite(cur_cor)] <- NA
#         best_ind <- which.max(abs(cur_cor))
#         cor_df <- data.frame(cor = cur_cor[, best_ind],
#                              feat_gene = colnames(ref)[best_ind],
#                              targ_gene = gene)
#         return(cor_df)
#     }, .parallel=TRUE) %>% mutate(dset = inputs$cur_dset)
#     return(res)
# }

get_top_biomarker_cors <- function(targ_dsets, dep_data, feature_mat, IDs_to_symbols, related_genes, use_bayes, min_feature_var = 0.01, cache_file = NULL) {

    #find top correlated biomarker for each gene dep in a set of dep datasets
    # full_res <- ldply(targ_dsets, function(cur_dset) {
    full_res <- data.frame()
    for (cur_dset in targ_dsets) {
        #loop over datasets and get stats
        print(cur_dset)
        dep_mat <- dep_data[[cur_dset]]$gene_scores
        if (dep_data[[cur_dset]]$mod_type == 'D2' & use_bayes) {
            SD_mat <- dep_data[[cur_dset]]$gene_score_SDs
        } else {
            SD_mat <- NULL
        }

        common_CLs <- intersect(rownames(feature_mat), rownames(dep_mat))

        cur_use_bayes <- dep_data[[cur_dset]]$mod_type == 'D2' & use_bayes

        dep_mat <- dep_mat[common_CLs,]
        if (!is.null(SD_mat)) {
            SD_mat <- SD_mat[common_CLs,]
        }

        #get feature mat
        cur_ref_mat <- feature_mat[common_CLs,]
        cur_ref_mat <- cur_ref_mat[, apply(cur_ref_mat, 2, var, na.rm=T) > min_feature_var] #require minimum feature variance
        test_genes <- colnames(dep_mat)

        input_cache_file <- '~/CPDS/demeter2/results/temp_input_cache.rds'
        output_res_file <- '~/CPDS/demeter2/results/temp_output_cache.rds'
        inputs <- list(ref = cur_ref_mat, dep = dep_mat, SD = SD_mat, genes = test_genes, cur_use_bayes = cur_use_bayes, cur_dset = cur_dset)

        #isolate this part in separate script to avoid memory issues with foreach parallelization
        write_rds(inputs, input_cache_file)
        source('~/CPDS/packages/demeter2_pub/scripts/get_cors_internal.R')
        cur_res <- read_rds(output_res_file)

        full_res %<>% rbind(cur_res)
        # full_res %<>% rbind(get_cors_internal(cur_ref_mat, dep_mat, SD_mat, test_genes, cur_use_bayes, cur_dset))
        if(!is.null(cache_file)) {
            write_rds(full_res, cache_file)
        }
    }

    full_res %<>% filter(!is.na(feat_gene)) %>%
        mutate(feat_gene = as.character(feat_gene),
               targ_gene = as.character(targ_gene))

    full_res %<>% left_join(related_genes, by = c('feat_gene' = 'partner_ID', 'targ_gene' = 'target_ID'))

    full_res %<>% mutate(source = ifelse(feat_gene == targ_gene, 'same', source))
    return(full_res)
    #
    # #determine whether each feature was a related gene or not
    # rel_an <- adply(full_res, 1, function(row) {
    #     cur_rel <- related_genes %>%
    #         filter(target == row$targ_sym) %>%
    #         .[['partner']]
    #     return(row %>% mutate(is_related = row$feat_sym %in% cur_rel))
    # }, .parallel = FALSE) %>% mutate(is_same = feat_gene == targ_gene)
    #
    # rel_an %<>% mutate(is_related = ifelse(is_same, TRUE, is_related))
    # return(rel_an)
}


get_gene_avgs <- function(dep_data, target_dsets, pos_cons, neg_cons, use_bayes = TRUE) {
    ## Calcultate per-gene averages for set of models
    ldply(target_dsets, function(cur_dset) {
        if (dep_data[[cur_dset]]$mod_type == 'D2' & use_bayes) {
            cur_gene_avgs <- sapply(seq(ncol(dep_data[[cur_dset]]$gene_scores)), function(colind) {
                wtd.mean(dep_data[[cur_dset]]$gene_scores[, colind],
                         1/dep_data[[cur_dset]]$gene_score_SDs[, colind]^2)
            })
        } else {
            cur_gene_avgs <- colMeans(dep_data[[cur_dset]]$gene_scores, na.rm=T)
        }
        return(data.frame(Gene = colnames(dep_data[[cur_dset]]$gene_scores),
                          avg_score = cur_gene_avgs,
                          dset = cur_dset))
    }) %>% mutate(gene_type = get_gene_type(Gene, pos_cons, neg_cons))
}

get_gene_SDs <- function(dep_data, target_dsets, pos_cons, neg_cons, use_bayes = TRUE) {
    gene_SDs <- ldply(target_dsets, function(cur_dset) {
        #Calculate per-gene SDs for set of models
        if (dep_data[[cur_dset]]$mod_type == 'D2' & use_bayes) {
            cur_gene_sds <- sapply(seq(ncol(dep_data[[cur_dset]]$gene_scores)), function(colind) {
                if (sum(!is.na(dep_data[[cur_dset]]$gene_score_SDs[, colind])) < 3) {
                    return(NA)
                } else {
                    wtd.var(dep_data[[cur_dset]]$gene_scores[, colind],
                            1/dep_data[[cur_dset]]$gene_score_SDs[, colind]^2)
                }
            }) %>% sqrt()
        } else {
            cur_gene_sds <- apply(dep_data[[cur_dset]]$gene_scores, 2, sd, na.rm=T)
        }
        return(data.frame(Gene = colnames(dep_data[[cur_dset]]$gene_scores),
                          SD_score = cur_gene_sds,
                          dset = cur_dset))
    }) %>% mutate(gene_type = get_gene_type(Gene, pos_cons, neg_cons))
}

data_available <- function(feature_genes, feature_types, feature_data) {
    #helper function to find out whether feature data is available for gene genes
    is_avail <- rep(NA, length.out = length(feature_genes))
    is_avail[feature_types == 'CN'] <- feature_genes[feature_types == 'CN'] %in% colnames(feature_data$CN)
    is_avail[feature_types == 'GE'] <- feature_genes[feature_types == 'GE'] %in% colnames(feature_data$GE)
    is_avail[feature_types == 'MUT_DAM'] <- feature_genes[feature_types == 'MUT_DAM'] %in% colnames(feature_data$MUT_DAM)
    is_avail[feature_types == 'MUT_HOT'] <- feature_genes[feature_types == 'MUT_HOT'] %in% colnames(feature_data$MUT_HOT)
    is_avail[feature_types == 'MUT_MIS'] <- feature_genes[feature_types == 'MUT_MIS'] %in% colnames(feature_data$MUT_MIS)
    return(is_avail)
}

get_target_biomarker_dep_cors <- function(target_dsets, benchmark_set, feature_data, gene_set = NULL, CL_set = NULL, use_bayes = TRUE) {
    #evaluate correlations for a set of benchmark dep-biomarker relationships
    if (!is.null(gene_set)) {
        benchmark_set %<>% filter(Dep_Gene_ID %in% gene_set)
    }
    bench_res <- adply(benchmark_set, 1, function(bench) {
        if (bench$feat_type == 'CN') {
            feature <- feature_data$CN[, bench$Feature_gene_ID]
        } else if (bench$feat_type == 'GE') {
            feature <- feature_data$GE[, bench$Feature_gene_ID]
        } else if (bench$feat_type == 'MUT_HOT') {
            feature <- feature_data$MUT_HOT[, bench$Feature_gene_ID]
        } else if (bench$feat_type == 'MUT_MIS') {
            feature <- feature_data$MUT_MIS[, bench$Feature_gene_ID]
        } else if (bench$feat_type == 'MUT_DAM') {
            feature <- feature_data$MUT_DAM[, bench$Feature_gene_ID]
        }

        if (sum(!is.na(feature) < 5)) {
            return(data.frame())
        }
        all_cors <- ldply(target_dsets %>% set_names(target_dsets), function(dname) {
           if (!(bench$Dep_Gene_ID %in% colnames(dep_data[[dname]]$gene_scores))) {
                return(data.frame(cor = NA))
            }
            dep <- dep_data[[dname]]$gene_scores[, bench$Dep_Gene_ID]
            cur_CLs <- intersect(names(dep), names(feature))
            cur_CLs <- cur_CLs[!is.na(dep[cur_CLs])]
            if (!is.null(CL_set)) {
                cur_CLs %<>% intersect(CL_set)
            }
            if (dep_data[[dname]]$mod_type == 'D2' & use_bayes) {
                dep_SD <- dep_data[[dname]]$gene_score_SDs[, bench$Dep_Gene_ID]
                if (var(feature[cur_CLs], na.rm=T) == 0) {
                    res <- data.frame(cor = NA, p.value = NA)
                } else {
                  c <- wtd.cor(dep[cur_CLs], feature[cur_CLs],
                               weight = 1/dep_SD[cur_CLs]^2)
                  pvalue <- c[1, 'p.value']
                  ptype <- 'cor'
                  if (bench$feat_type %in% c('MUT_HOT', 'MUT_MIS', 'MUT_DAM')) {
                    inG <- cur_CLs[feature[cur_CLs] == 1]
                    outG <- cur_CLs[feature[cur_CLs] == 0]
                    if (length(inG) >= 2 & length(outG) >= 2) {
                      tres <- wtd.t.test(dep[inG], dep[outG], 
                                    weight = 1/dep_SD[inG]^2, weighty = 1/dep_SD[outG]^2)
                    pvalue <- tres$coefficients[['p.value']]
                    } else {
                      pvalue <- NA
                    }
                    ptype <- 'ttest'
                  }
                  res <- data.frame(cor = c[1, 'correlation'], p.value = pvalue, ptype = ptype)
                }
            } else {
                c <- cor.test(dep[cur_CLs], feature[cur_CLs], use = 'pairwise.complete.obs')
                pvalue <- c$p.value
                ptype <- 'cor'
                if (bench$feat_type %in% c('MUT_HOT', 'MUT_MIS', 'MUT_DAM')) {
                  inG <- cur_CLs[feature[cur_CLs] == 1]
                  outG <- cur_CLs[feature[cur_CLs] == 0] 
                  if (length(inG) >= 2 & length(outG) >= 2) {
                    tres <- t.test(dep[inG], dep[outG])
                    pvalue <- tres$p.value
                  } else {
                    pvalue <- NA
                  }
                  ptype <- 'ttest'
                }
                res <- data.frame(cor = c$estimate, p.value = pvalue, ptype = ptype)
            }
            return(res)
        }, .id = 'model')
        return(all_cors)
    })
    return(bench_res)
}



get_gene_avg_stats <- function(dep_data, gene_avgs, target_dsets, only_overlapping_genes, use_bayes = TRUE) {
    #Calculate pos-neg sep stats for per-gene avg estimates
    ldply(target_dsets, function(cur_dset) {
        if (dep_data[[cur_dset]]$type == 'abs') {
            cur_gene_avgs <- gene_avgs %>% filter(dset == cur_dset)
            if (only_overlapping_genes) {
                cur_gene_avgs %<>% filter(Gene %in% common_dep_genes)
            }
            cur_NE_calls <- get_gene_calls_NE(cur_gene_avgs$gene_type)
            used <- which(!is.na(cur_NE_calls) & !is.na(cur_gene_avgs$avg_score))
            if (dep_data[[cur_dset]]$mod_type == 'D2' & use_bayes) {
                gene_avg_var <- colMeans(dep_data[[cur_dset]]$gene_score_SDs^2, na.rm=T)
                cur_gene_avgs %<>% left_join(data.frame(Gene = names(gene_avg_var), gene_var = gene_avg_var), by = 'Gene')
                ssmd <- get_SSMD(cur_gene_avgs$avg_score, cur_NE_calls, 1/cur_gene_avgs$gene_var)
                pr <- pr.curve(-cur_gene_avgs$avg_score[which(cur_NE_calls)],
                               -cur_gene_avgs$avg_score[which(!cur_NE_calls)],
                               weights.class0=1/cur_gene_avgs$gene_var[which(cur_NE_calls)],
                               weights.class1=1/cur_gene_avgs$gene_var[which(!cur_NE_calls)])
            } else {
                ssmd <- get_SSMD(cur_gene_avgs$avg_score, cur_NE_calls)
                pr <- pr.curve(-cur_gene_avgs$avg_score[which(cur_NE_calls)],
                               -cur_gene_avgs$avg_score[which(!cur_NE_calls)])
            }
            rob_ssmd <- get_SSMD(cur_gene_avgs$avg_score, cur_NE_calls, robust = TRUE)
            return(data.frame(ssmd = ssmd, rob_ssmd = rob_ssmd, pr = pr$auc.integral, dset = cur_dset))
        }
    })
}

#calculate gene dep stats for each cell line
get_per_CL_stats <- function(dep_data, target_dsets, only_overlapping_genes, neg_control_type,
                             GE_mat, pos_cons, neg_cons, use_bayes = TRUE) {
    #Genes present in all target datasets
    common_dep_genes <- Reduce(intersect, llply(target_dsets, function(dname) {
        colnames(dep_data[[dname]]$gene_scores)
        }))

    ldply(target_dsets, function(cur_dset) {
        cur_CLs <- rownames(dep_data[[cur_dset]]$gene_scores)
        if (only_overlapping_genes) {
            cur_genes <- common_dep_genes
        } else {
            cur_genes <- colnames(dep_data[[cur_dset]]$gene_scores)
        }
        cur_gene_scores <- dep_data[[cur_dset]]$gene_scores[, cur_genes]
        if (dep_data[[cur_dset]]$mod_type == 'D2' & use_bayes) {
            cur_gene_score_SDs <- dep_data[[cur_dset]]$gene_score_SDs[, cur_genes]
        }
        cur_gene_types <- get_gene_type(cur_genes, pos_cons, neg_cons)
        cur_GE <- GE_mat[match(cur_CLs, rownames(GE_mat)), match(cur_genes, colnames(GE_mat))]
        rownames(cur_GE) <- cur_CLs; colnames(cur_GE) <- cur_genes
        cur_per_CL_stats <- ldply(cur_CLs, function(cur_CL) {
            has_GE <- cur_CL %in% rownames(GE_mat)
            pos_med <- median(cur_gene_scores[cur_CL, cur_gene_types == 'essential'], na.rm=T)
            pos_mean <- mean(cur_gene_scores[cur_CL, cur_gene_types == 'essential'], na.rm=T)
            if (neg_control_type == 'unexpressed') {
                if (!has_GE) {
                    return(data.frame(pr = NA, roc = NA, ssmd = NA, rob_ssmd = NA, pos_med = pos_med, pos_mean = pos_mean, CCLE_ID = cur_CL))
                }
                cur_calls <- get_gene_calls_UE(cur_gene_types, cur_GE[cur_CL,], GE_thresh)
            } else if (neg_control_type == 'non_essential') {
                cur_calls <- get_gene_calls_NE(cur_gene_types)
            }
            used <- which(!is.na(cur_calls) & !is.na(cur_gene_scores[cur_CL,]))
            if (dep_data[[cur_dset]]$mod_type == 'D2' & use_bayes) {
                pr <- pr.curve(-cur_gene_scores[cur_CL,used[cur_calls[used]]],
                               -cur_gene_scores[cur_CL,used[!cur_calls[used]]],
                               weights.class0=1/cur_gene_score_SDs[cur_CL,used[cur_calls[used]]]^2,
                               weights.class1=1/cur_gene_score_SDs[cur_CL,used[!cur_calls[used]]]^2)
                roc <- roc.curve(-cur_gene_scores[cur_CL,used[cur_calls[used]]],
                                 -cur_gene_scores[cur_CL,used[!cur_calls[used]]],
                                 weights.class0=1/cur_gene_score_SDs[cur_CL,used[cur_calls[used]]]^2,
                                 weights.class1=1/cur_gene_score_SDs[cur_CL,used[!cur_calls[used]]]^2)
                ssmd <- get_SSMD(cur_gene_scores[cur_CL,], cur_calls,
                                 1/cur_gene_score_SDs[cur_CL,]^2)
            } else {
                pr <- pr.curve(-cur_gene_scores[cur_CL,used], weights.class0=cur_calls[used])
                roc <- roc.curve(-cur_gene_scores[cur_CL,used], weights.class0=cur_calls[used])
                ssmd <- get_SSMD(cur_gene_scores[cur_CL,], cur_calls)
            }
            rob_ssmd <- get_SSMD(cur_gene_scores[cur_CL,], cur_calls, robust = TRUE)
            return(data.frame(pr = pr$auc.integral,
                              roc = roc$auc,
                              ssmd = ssmd,
                              rob_ssmd = rob_ssmd,
                              pos_med = pos_med,
                              pos_mean = pos_mean,
                              CCLE_ID = cur_CL))
        }) %>% mutate(dset = cur_dset)
        return(cur_per_CL_stats)
    })
}

get_dep_corrs <- function(dep_data, target_dsets, target_genes, gene_avgs, use_bayes = TRUE) {
    dep_cor_res <- ldply(target_dsets, function(cur_dset) {
    print(cur_dset) 
    if (use_bayes & !is.null(dep_data[[cur_dset]]$gene_score_SDs)) {
      #use average gene-level variance for precision weighting
      CL_avg_var <- rowMeans(dep_data[[cur_dset]]$gene_score_SDs^2, na.rm=T)
      cur_weight <- 1/CL_avg_var
      c_mat <- wtd.cors(dep_data[[cur_dset]]$gene_scores[, target_genes], weight = cur_weight)
    } else {
      c_mat <- cor(dep_data[[cur_dset]]$gene_scores[, target_genes], use = 'pairwise.complete.obs')
    }
    #calculate average essentiality of each gene pair in the matrix
    if (grepl('Ach', cur_dset)) {
      avg_ess <- gene_avgs %>% 
        filter(dset == 'D2_Ach')
    } else if (grepl('DRIVE', cur_dset)) {
      avg_ess <- gene_avgs %>% 
        filter(dset == 'D2_DRIVE')
    }
    avg_ess <- avg_ess[match(target_genes, avg_ess$Gene), 'avg_score']
    a <- outer(avg_ess, avg_ess, '+') / 2
    a <- a[lower.tri(a)]
    c_mat <- c_mat[lower.tri(c_mat)]
    
    data.frame(cor = c_mat, a = a, dset = cur_dset)
  })
  return(dep_cor_res)
}


normalize_genescores_global_z <- function(cur_data) {
  if (!is.null(cur_data$gene_score_SDs)) {
    cur_data$gene_score_SDs <- cur_data$gene_score_SDs / sd(cur_data$gene_scores, na.rm=T)
  }            
  cur_data$gene_scores <- cur_data$gene_scores / sd(cur_data$gene_scores, na.rm=T)
  return(cur_data)
}

normalize_genescores_posneg <- function(cur_data, hart_ess, hart_non_ess, cur_GA, normalize_D2_per_CL = FALSE) {
  if (normalize_D2_per_CL & cur_data$mod_type == 'D2') {
    cur_pos_set <- as.character(hart_ess[hart_ess %in% colnames(cur_data$gene_scores)])
    cur_neg_set <- as.character(hart_non_ess[hart_non_ess %in% colnames(cur_data$gene_scores)])
    cur_pos_medians <- apply(cur_data$gene_scores[, cur_pos_set], 1, median, na.rm=T)
    cur_neg_medians <- apply(cur_data$gene_scores[, cur_neg_set], 1, median, na.rm=T)
    cur_scale_facs <- pmax(cur_neg_medians - cur_pos_medians, 0.001)
    if (!is.null(cur_data$gene_score_SDs)) {
      cur_data$gene_score_SDs %<>% t() %>% scale(center = FALSE, scale = cur_scale_facs) %>% t()
    }
    cur_data$gene_scores %<>% t() %>% scale(center = cur_neg_medians, scale = cur_scale_facs) %>% t()
  } else {
    # cur_GA <- cur_gene_avgs %>% filter(dset == cur_dset)
    pos_median <- cur_GA %>% filter(gene_type == 'essential') %>% .[['avg_score']] %>% median(na.rm=T)
    neg_median <- cur_GA %>% filter(gene_type == 'non_essential') %>% .[['avg_score']] %>% median(na.rm=T)
    if (!is.null(cur_data$gene_score_SDs)) {
      cur_data$gene_score_SDs <- cur_data$gene_score_SDs / (neg_median - pos_median)
    }
    cur_data$gene_scores <- (cur_data$gene_scores - neg_median) / (neg_median - pos_median)
  }
  return(cur_data)
}


print_paired_two_group_stats <- function(var1, var2) {
  print(sprintf(
    'Wilcox p: %.3g',
    wilcox.test(var1, var2, paired = TRUE)$p.value
  ))
  
  print(sprintf(
    'Num used: %d',
    sum(!is.na(var1) & !is.na(var2))
  ))
  
  print(sprintf(
    'Median improvement: %.3f',
    median((var1 - var2)/var1, na.rm=TRUE)
  ))
}


###----------  PLOTTING HELPERS #########

make_net_plot <- function(cur_df, C_mat) {
  max_nodes <- 25
  min_size <- 10
  max_size <- 30
  cor_z_thresh <- 4
  
  #get max_nodes top correlated dependencies
  node_df <- cur_df %>% 
    arrange(desc(abs(cor))) %>% 
    head(max_nodes) %>%
    mutate(id = Gene,
           label = Gene_symbol)
  
  node_df$color <- 'black'
  node_df %<>% mutate(
    color = ifelse(is_related, 'pink', 'lightgray'),
    color = ifelse(Gene == q_gene_ID, 'skyblue', color)
  )
  
  edge_df <- 
    C_mat[node_df$Gene, node_df$Gene] %>% 
    melt() %>% 
    set_colnames(c('from', 'to', 'cor')) %>% 
    mutate(cor_rank = percent_rank(-abs(cor)),
           cor = abs(cor)) %>% 
    # filter(cor_rank < cor_rank_thresh, from != to)
    filter(cor > cor_z_thresh, from != to) 
  
  edge_widths <- edge_df$cor/sd(edge_df$cor, na.rm=T)
  
  #build initial graph
  graph <- igraph::graph_from_data_frame(d = edge_df,
                                         vertices = node_df,
                                         directed = FALSE)
  
  #remove any components that are not connected
  is_connected <- sapply(node_df$id, function(target_ID) {
    if (target_ID == q_gene_ID) {
      return(TRUE)
    } else {
      return(igraph::edge_connectivity(graph, source = q_gene_ID, target = target_ID) > 0)
    }
  })
  node_df %<>% 
    mutate(is_connected = is_connected) %>% 
    filter(is_connected)
  edge_df %<>% filter(to %in% node_df$id, from %in% node_df$id)
  
  #rebuild graph
  graph <- igraph::graph_from_data_frame(d = edge_df,
                                         vertices = node_df,
                                         directed = FALSE)
  
  #node size represents how common-essential each gene is
  vertex_size <- -node_df[["avg_score"]] 
  vertex_size <- vertex_size*(max_size - min_size) + min_size
  vertex_size[is.na(vertex_size)] <- mean(vertex_size, na.rm=T)
  
  igraph::plot.igraph(graph,
                      vertex.label.dist = 0,
                      vertex.label = node_df$label,
                      vertex.label.color = 'black',
                      vertex.frame.color = 'black',
                      vertex.color = node_df$color,
                      vertex.size = vertex_size*0.75,
                      vertex.alpha = 0.75,
                      vertex.shape = 'circle',
                      edge.arrow.size = 0,
                      edge.curved = 0.1,
                      layout = layout_nicely,
                      vertex.label.dist=1,
                      vertex.label.font=2,
                      vertex.label.cex = 0.7,
                      edge.width=edge_df$cor*0.2)
}

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

#https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=12, base_family="Helvetica") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.2), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(),
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                legend.text = element_text(size = rel(1.2)),
                legend.direction = "horizontal",
                legend.key.size= unit(0.3, "cm"),
                legend.margin = unit(0, "cm"),
                legend.title = element_text(face="italic"),
                plot.margin=unit(c(5,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))

}

scale_fill_Publication <- function(...){
    library(scales)
    discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_Publication <- function(...){
    library(scales)
    discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

gene_type_pal <- c(
    essential = 'red',
    `non-essential` = 'blue',
    other = 'grey'
)

# data_type_pal <- c(
#     D2_Achilles = 'forestgreen',
#     D2_DRIVE= 'orange',
#     D2_combined = 'firebrick',
#     GA_combined = 'purple',
#     CERES_CRISPR = 'gray',
#     Achilles = 'forestgreen',
#     DRIVE= 'orange',
#     CERES = 'gray'
# )

data_type_pal <- c(
    `D2 Ach` = 'forestgreen',
    `D2 DRIVE`= 'orange',
    `D2 comb` = 'firebrick',
    `GA comb` = 'purple',
    `D1 comb` = 'lightblue',
    CRISPR = 'gray',
    Achilles = 'forestgreen',
    DRIVE= 'orange',
    CERES = 'gray'
)

mod_type_pal <- c(
    D2 = "#386cb0",
    D2noUnc = 'lightblue',
    GA = "#fdb462",
    RSA = '#662506',
    ATARiS = "#7fc97f",
    D1 = "#ef3b2c",
    CERES = "gray",
    MAGeCK = 'darkmagenta',
    `D1-PC` = 'lightgreen',
    `GA-indnorm` = 'magenta'
)

