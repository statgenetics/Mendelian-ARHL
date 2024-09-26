# Functions to generate data
process_trait_data <- function(file_path, trait_name) {
    df <- read.delim(file_path, sep = " ", comment.char = "#", header = TRUE)
    df <- df[df$CHROM != "CHROM", ]
    df$CHROM <- as.numeric(df$CHROM)
    df$LOG10P <- as.numeric(df$LOG10P)
    df$GENPOS <- as.numeric(df$GENPOS)
    df$ID_gene <- sapply(strsplit(df$ID, "\\."), `[`, 1)

    filtered_data_list <- list()
    for (i in unique(df$ALLELE1)) {
        for (j in unique(df$TEST)) {
            tmp <- df %>% filter(ALLELE1 == i & TEST == j)
            filtered_data_list[[paste0(i, "_", j)]] <- tmp
        }
    }

    Mask1_skato <- filtered_data_list[["Mask1.0.005_ADD-SKATO"]] %>%
        mutate(LOG10P = as.numeric(LOG10P), P = 10^(-LOG10P))
    Mask2_skato <- filtered_data_list[["Mask2.0.005_ADD-SKATO"]] %>%
        mutate(LOG10P = as.numeric(LOG10P), P = 10^(-LOG10P))
    Mask1_burden <- filtered_data_list[["Mask1.0.005_ADD"]] %>%
        mutate(LOG10P = as.numeric(LOG10P), P = 10^(-LOG10P))
    Mask2_burden <- filtered_data_list[["Mask2.0.005_ADD"]] %>%
        mutate(LOG10P = as.numeric(LOG10P), P = 10^(-LOG10P))

    joined_df_skato <- merge(Mask1_skato, Mask2_skato, by = c("CHROM", "GENPOS", "ALLELE0", "ID_gene"), all.x = TRUE, all.y = TRUE)
    joined_df_skato$ID_gene <- factor(joined_df_skato$ID_gene)

    joined_df_burden <- merge(Mask1_burden, Mask2_burden, by = c("CHROM", "GENPOS", "ALLELE0", "ID_gene"), all.x = TRUE, all.y = TRUE)
    joined_df_burden$ID_gene <- factor(joined_df_burden$ID_gene)
    
    final_df_pLoF_skato <- joined_df_skato %>%
        filter(ID_gene != "TBCEL-TECTA") %>%
        select(ID_gene, CHROM, GENPOS, P.x) %>%
        rename(Chromosome = CHROM, SNP = ID_gene, Position = GENPOS, !!paste0(trait_name, "_pLoF") := P.x) %>%
        select(SNP, Chromosome, Position, !!paste0(trait_name, "_pLoF"))

    final_df_pLoF_missense_splicing_skato <- joined_df_skato %>%
        filter(ID_gene != "TBCEL-TECTA") %>%
        select(ID_gene, CHROM, GENPOS, P.y) %>%
        rename(Chromosome = CHROM, SNP = ID_gene, Position = GENPOS, !!paste0(trait_name, "_pLoF_missense_splice") := P.y) %>%
        select(SNP, Chromosome, Position, !!paste0(trait_name, "_pLoF_missense_splice"))

    final_df_pLoF_burden <- joined_df_burden %>%
        filter(ID_gene != "TBCEL-TECTA") %>%
        select(ID_gene, CHROM, GENPOS, P.x) %>%
        rename(Chromosome = CHROM, SNP = ID_gene, Position = GENPOS, !!paste0(trait_name, "_pLoF") := P.x) %>%
        select(SNP, Chromosome, Position, !!paste0(trait_name, "_pLoF"))

    final_df_pLoF_missense_splicing_burden <- joined_df_burden %>%
        filter(ID_gene != "TBCEL-TECTA") %>%
        select(ID_gene, CHROM, GENPOS, P.y) %>%
        rename(Chromosome = CHROM, SNP = ID_gene, Position = GENPOS, !!paste0(trait_name, "_pLoF_missense_splice") := P.y) %>%
        select(SNP, Chromosome, Position, !!paste0(trait_name, "_pLoF_missense_splice"))

    list(final_df_pLoF_skato = final_df_pLoF_skato, final_df_pLoF_missense_splicing_skato = final_df_pLoF_missense_splicing_skato, final_df_pLoF_burden = final_df_pLoF_burden, final_df_pLoF_missense_splicing_burden = final_df_pLoF_missense_splicing_burden)
}


format_variant_df <- function(filename) {
    df <- fread(filename)
    df <- df %>%
        filter(CHROM != "CHROM") %>% # Remove header rows if any
        mutate(
            A1FREQ = as.numeric(A1FREQ),
            LOG10P_numeric = as.numeric(LOG10P),
            P = 10^(-LOG10P_numeric),
            MAF = if_else(A1FREQ <= 0.5, A1FREQ, 1 - A1FREQ)
        ) %>%
        rename(SNP = ID, Chromosome = CHROM, Position = GENPOS) %>%
        select(SNP, Chromosome, Position, P, A1FREQ, MAF)
    return(df)
}


import_variant_data <- function(pheno, data_field) {
    df <- format_variant_df(paste0(
        "/home/dmc2245/UKBiobank/RAP/results/autosomal/univariate/",
        pheno, "_070824/ref_last/", data_field, ".regenie"
    ))
    df_all <- df %>%
        select(!A1FREQ) %>%
        rename(!!pheno := P)
    df_rare <- df %>% filter(MAF <= 0.005)
    df_common <- df %>% filter(MAF > 0.005)
    return(list(all = df_all, rare = df_rare, common = df_common))
}

import_variant_data_all_columns <- function(pheno, data_field) {
    df <- fread(paste0(
        "/home/dmc2245/UKBiobank/RAP/results/autosomal/univariate/",
        pheno, "_070824/ref_last/", data_field, ".regenie"))
    df <- df %>%
        filter(CHROM != "CHROM") %>% # Remove header rows if any
        mutate(
            A1FREQ = as.numeric(A1FREQ),
            LOG10P_numeric = as.numeric(LOG10P),
            P = 10^(-LOG10P_numeric),
            MAF = if_else(A1FREQ <= 0.5, A1FREQ, 1 - A1FREQ)
        ) %>%
        rename(SNP = ID, Chromosome = CHROM, Position = GENPOS)
    df_all <- df
    df_rare <- df %>% filter(MAF <= 0.005)
    df_common <- df %>% filter(MAF > 0.005)
    return(list(all = df_all, rare = df_rare, common = df_common))
}


annotate_rsID <- function(df, df_rsID, remove_NA = TRUE) {
  df <- df %>%
    separate(SNP, into = c("Chromosome", "Position", "REF", "ALT"), sep = ":", remove = FALSE)
  df$Chromosome <- as.integer(df$Chromosome)
  df$Position <- as.integer(df$Position)
  df$BETA <- as.numeric(df$BETA)
  df$LOG10P_numeric <- as.numeric(df$LOG10P)
  merged_df <- merge(df, df_rsID,
                     by.x = c("Chromosome", "Position", "REF", "ALT"),
                     by.y = c("chr", "start_pos", "ref", "alt"),
                     all.x = TRUE)
  merged_df <- merged_df %>% select(!SNP) %>% rename(SNP = rsID)
    if (remove_NA){
        merged_df <- merged_df %>%
            filter(!is.na(SNP))
        }
  return(merged_df)
}


# Functions to generate figures
generate_manhattan_plot_gene <- function(df, plot_file_name, display = TRUE, format = "pdf") {
    # if display = TRUE, then only display in the notebook. Otherwise print to output file
    # format can be either png, jpg or pdf
    last_col_name <- colnames(df)[ncol(df)]
    SNPs_to_highlight <- df[df[[last_col_name]] < 2.5e-06, "SNP"]
    SNPs_to_highlight <- na.omit(SNPs_to_highlight)
    CMplot(df,
        plot.type = "m", LOG10 = TRUE, col = c("grey30", "grey60"),
        threshold = c(3.1e-7, 2.5e-6), threshold.lty = c(2, 1),
        threshold.lwd = c(1, 1), threshold.col = c("grey", "black"),
        highlight = SNPs_to_highlight, highlight.col = "red", highlight.cex = 1.5, highlight.text.cex = 1.8,
        highlight.text = SNPs_to_highlight, highlight.text.col = rep("red", length(SNPs_to_highlight)),
        chr.labels.angle = 0.01,
        amplify = FALSE, file = format, file.name = plot_file_name, dpi = 300,
        file.output = !display, verbose = TRUE, width = 14, height = 6
    )
}


generate_qq_plot_gene <- function(df, plot_file_name, display = TRUE, format = "pdf") {
    CMplot(df,
        plot.type = "q", box = FALSE, file = format, file.name = plot_file_name, dpi = 300,
        conf.int = TRUE, conf.int.col = NULL, threshold.col = "red", threshold.lty = 2,
        file.output = !display, verbose = TRUE, width = 5, height = 5, main = ""
    )
}

generate_qq_plot_variant <- function(df, plot_file_name, display = TRUE, format = "pdf") {
    df <- df %>% select(!MAF)
    CMplot(df,
        plot.type = "q", box = FALSE, file = format, file.name = plot_file_name, dpi = 300,
        conf.int = TRUE, conf.int.col = NULL, threshold.col = "red", threshold.lty = 2,
        file.output = !display, verbose = TRUE, width = 5, height = 5, main = ""
    )
}

generate_manhattan_plot_variants_annot_gene <- function(pheno, df, plot_file_name, display = TRUE, format = "pdf") {
    # print(Sys.time())
    # time_1 = Sys.time()
    gene_file <- paste0(
        "/home/dmc2245/UKBiobank/RAP/results/autosomal/univariate/",
        pheno, "_070824/ref_last/", pheno, "_pval5e-08_rap1.hg38.annotated.csv"
    )
    gene_annotation <- fread(gene_file)
    gene_annotation <- gene_annotation %>%
        select(ID, Gene.refGene) %>%
        rename(SNP = ID)
    variant_gene_info <- left_join(df, gene_annotation, by = "SNP") %>%
        mutate(gene = ifelse(is.na(Gene.refGene), NA, Gene.refGene)) %>%
        rename(
            ID = SNP,
            SNP = gene
        ) %>%
        select(SNP, Chromosome, Position, P)
    
    variant_gene_info_min_p <- variant_gene_info %>%
        filter(P < 5e-08) %>%
        group_by(SNP) %>%
        filter(P == min(P)) %>%
        ungroup()

    variant_gene_info <- variant_gene_info %>%
        left_join(variant_gene_info_min_p %>% select(SNP, MinP = P), by = "SNP") %>%
        mutate(SNP = ifelse(P != MinP, NA, SNP)) %>%
        select(-MinP)

    genes_to_highlight <- variant_gene_info %>%
        filter(P < 5e-08) %>%
        filter(!is.na(SNP)) %>%
        pull(SNP)
    genes_to_highlight <- na.omit(genes_to_highlight)
    print(paste0("===== Genes to highlight for ", pheno, " ====="))
    print(genes_to_highlight)
    colnames(variant_gene_info)[ncol(variant_gene_info)] <- pheno

    CMplot(variant_gene_info,
        plot.type = "m", LOG10 = TRUE, col = c("grey30", "grey60"),
        axis.cex = 1.2, lab.cex = 2.5, axis.lwd = 2.5,
        threshold = c(5e-08), threshold.lty = c(2),
        threshold.lwd = c(1), threshold.col = c("black"),
        highlight = genes_to_highlight, highlight.col = "red", highlight.cex = 1.5, highlight.text.cex = 1.8,
        highlight.text = genes_to_highlight, highlight.text.col = rep("red", length(genes_to_highlight)),
        chr.labels.angle = 0.01,
        amplify = FALSE, file = format, file.name = plot_file_name, dpi = 300,
        file.output = !display, verbose = TRUE, width = 14, height = 6
    )
}
