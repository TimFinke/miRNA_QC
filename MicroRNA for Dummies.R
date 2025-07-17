# ==================================================================================================
# Libraries & Functions
# ==================================================================================================
library("tidyverse")
library("haven")
library("edgeR")
library("openxlsx")
library("flextable")
library("corrplot")
library("ggcorrplot")
library("cowplot")


"%!in%" <- function(x, y) {!("%in%"(x, y))}

# Winsorization
# Origin: https://github.com/inDEPTHlab/epigenetics/blob/main/epigenetics/age0/epic/dnam_age0_epic_prep.R
"winsorize" <- function(x, 
                        probs = NULL, 
                        cutpoints = NULL , 
                        replace = c(cutpoints[1], cutpoints[2]), 
                        verbose = FALSE){
  dummy = is.integer(x)
  if (!is.null(probs)){
    stopifnot(is.null(cutpoints))
    stopifnot(length(probs)==2)
    cutpoints <- quantile(x, probs, type = 1, na.rm = TRUE)
  } else if (is.null(cutpoints)){
    l <- quantile(x, c(0.25, 0.50, 0.75), type = 1, na.rm = TRUE) 
    cutpoints <- c(l[1]-3*(l[3]-l[1]), l[3]+3*(l[3]-l[1]))  ### Default was Median+-3*IQR but has been changed to +-1*IQR+-3*IQR
  } else{
    stopifnot(length(cutpoints)==2)
  }
  if (is.integer(x)) cutpoints <- round(cutpoints)
  bottom <- x < cutpoints[1]
  top <- x > cutpoints[2]
  if (verbose){
    length <- length(x)
    message(paste(100*sum(bottom, na.rm = TRUE)/length,"% observations replaced at the bottom"))
    message(paste(100*sum(top, na.rm = TRUE)/length,"% observations replaced at the top"))
  }
  x[bottom] <- replace[1]
  x[top] <- replace[2]
  if (dummy){
    x <- as.integer(x)
  }
  x
}
# Make correlation matrix
"make_cor_matrix" <- function(dataset1, dataset2, 
                              title = title){
  # Bind into 1 dataframe
  xx <- cbind.data.frame(dataset1, dataset2)
  
  # correlations & P values
  corr_covs_PCs <- cor(xx, use = "pairwise.complete.obs")
  corr_covs_PCs_Pvals <- cor.mtest(xx, conf.level = .95)$p
  
  # Construct figure
  ggcorrplot(corr_covs_PCs, hc.order = FALSE,
             title = title,
             p.mat = corr_covs_PCs_Pvals,
             sig.level = .05,
             insig = "blank",
             outline.color = "white",
             type = "upper",
             ggtheme = ggplot2::theme_classic,
             colors = c("#E46726", "white", "#6D9EC1"),
             lab = TRUE
  ) +
    theme(plot.title = element_text(hjust = 0.5, size = 25))
}

# ==================================================================================================
# Loading raw data
# ==================================================================================================
# Path to folder
path <- "./data"
# path to archived data
if(!dir.exists(file.path(path, "data_clean/"))){
  dir.create(file.path(path, "data_clean/"))
}
if(!dir.exists(file.path(path, "data_clean/archived/"))){
  dir.create(file.path(path, "data_clean/archived/"))
}
# path to DataWiki
if(!dir.exists(file.path(path, "data_clean/DataWiki/"))){
  dir.create(file.path(path, "data_clean/DataWiki/"))
}
# path to archived data
if(!dir.exists(file.path(path, "data_clean/Midas/"))){
  dir.create(file.path(path, "data_clean/Midas/"))
}

# IDCs
df_QCres <- read_sav(file.path(path, "QCresults_miRNA_allPlates_IDC.sav"))
# general data
df_gen <- read_sav(file.path(path, "CHILD-ALLGENERALDATA_05122024.sav"))
# ethnicity data
df_ethv3 <- read_sav(file.path(path, "ethnicity/PCA_Selection GWAv3_revised def_October2022.sav"))[, c(1, 3, 5, 6)]
colnames(df_ethv3)[c(3,4)] <- c("GWAv3C1", "GWAv3C2")
# cell type methylation data
MethylEPIC1 <- read_sav(file.path(path, "Cell_types/Selection_GENR_MethylEPIC_release1_birth_20230717.sav"))
load(file.path(path, "Cell_types/GENR_EPICv1METH_birth_CellTypes_combined.RData"))
Methyl450k3 <- read_sav(file.path(path, "Cell_types/Selection_GENR_450kmeth_release3_birth_20230608.sav"))
Methyl450k3_Salas <- read.csv(file.path(path, "Cell_types/GENR_450kmeth_release3_birth_Salas.csv"))

Methyl450k3$SampleID <- Methyl450k3$Sample_ID
Methyl450k3_Salas$SampleID <- Methyl450k3_Salas$Rownames
Methyl450k3$Sample_ID <- NULL
Methyl450k3_Salas$Rownames <- NULL

Methyl_450k_merged <- left_join(Methyl450k3, Methyl450k3_Salas, by = "SampleID")
Methyl_EPIC_merged <- left_join(MethylEPIC1, GENR_EPICv1METH_birth_CellTypes_combined.data, by = "SampleID")
Methyl_450k_merged[, c(2:5)] <- NULL
Methyl_EPIC_merged[, c(2, 4)] <- NULL
Methyl_merged <- rbind.data.frame(Methyl_450k_merged, Methyl_EPIC_merged)

# Read in data per plate
dataframes <- list()

n_plates <- 17 # first 17 plates, plate 18 newest version is under the name plate 20!!

for (i in c(1:n_plates, 20)) {
  ### Read in haemolysis
  ### Define parameters first, i.e., file name etc
  haemolysis_file <- paste0("Plate", i, "_Haemolysis.csv")
  name_haemolysis <- paste0("Plate", i, "_Haemolysis")
  ### Read in haemolysis matrix (8x12), translate as a matrix
  haemolysis <- as.matrix(read.csv(file = file.path(path, haemolysis_file), header = T, row.names = 1))
  ### Assign a well name to each well; QC step
  names(haemolysis) <- paste0(LETTERS[rep((1:8), 12)], rep(1:12, each = 8))
  ### Translate (8x12) matrix into a 96-element vector
  haemolysis <- c(haemolysis)
  assign(x = name_haemolysis, value = haemolysis)
  
  ### Read in concentrations
  ### Define parameters first, i.e., file name etc
  concentration_file <- paste0("Plate", i, "_Concentrations.csv")
  name_concentration <- paste0("Plate", i, "_Concentration")
  concentration <- read.csv(file = file.path(path, concentration_file), header = F, row.names = 1)
  colnames(concentration) <- c("Concentration")
  concentration$Well <- as.character(paste0(LETTERS[rep((1:8), 12)], rep(1:12, each = 8)))
  ### Translate into a 96-element vector
  concentration <- c(concentration)
  assign(x = name_concentration, value = concentration)
  
  ### Read in counts file
  ### Define parameters first, i.e., file name etc
  name_file <- paste0("Plate", i, ".csv")
  name <- paste0("Plate", i)
  table <- read.csv(file = file.path(path, name_file), header = T)
  
  ### Add haemolysis information into the last column, named #Haemolysis
  table$Haemolysis <- haemolysis
  table$Concentration <- concentration$Concentration
  
  assign(x = name, value = table)
  dataframes[[i]] <- get(name)
}

# All plates in 1 dataframe
df <- do.call("rbind", dataframes)
# Refactor Wells into each well being its own factor
df$Well <- factor(df$Well, levels = paste0(LETTERS[rep((1:8), 12)], rep(1:12, each = 8)))

### merge df with df_QCres
# Matching colnames for df & df_QCres
df_QCres$Sample.Name <- df_QCres$Sample
df_QCres$Sample <- NULL
# Join on Sample.Name and Well
df <- df %>% left_join(., df_QCres, by = c("Sample.Name", "Well"))

# save original data
original_data <- df[, c(2108, 1:2103)]
write.xlsx(original_data, file.path(path, paste0("data_clean/archived/miRNA_GenR_original_", format(Sys.Date(), "%d_%m_%Y"), ".xlsx")), 
           colNames = TRUE, rowNames = FALSE,
           overwrite = FALSE)


# save dataset covariates
covariate_cols <- c("IDC", "Sample.Name", "Plaat", "Well", 
                    "Concentration", "Haemolysis", "QCStatus", "POS", 
                    "RawTotalCounts", "RSD", "indexplaat")
covariates <- df[, covariate_cols]

write.xlsx(covariates, file.path(path, paste0("data_clean/Midas/miRNA_GenR_covariates_", format(Sys.Date(), "%d_%m_%Y"), ".xlsx")),
           sheetName = c("covariates"),
           colNames = TRUE, rowNames = FALSE,
           overwrite = FALSE)

# ==================================================================================================
# Applying Quality Control
# ==================================================================================================
### Implement filters here for HK- & control genes, controls, QC status failures etc
# dropping HK and control genes
df <- df[, -c(2:20)]
# dropping controls
controls <- df %>% filter(., grepl("pl", Sample.Name) | 
                            grepl("Control", Sample.Name) | 
                            grepl("control", Sample.Name) | 
                            grepl("ctrl", Sample.Name)
)
df <- df %>% filter(., Sample.Name %!in% controls$Sample.Name)
# dropping all samples that did not pass QCStatus1 or QCStatus2
df_QC_Fail <- df %>% filter(., QCStatus != "PASS") # 18, of which 3 controls and 15 participants, hence N = 1695
df <- df %>% filter(., QCStatus == "PASS")

# Make df technical covariates (columns 2085 - 2094)
cols_tech_covs <- which(names(df) %in% c(covariate_cols[-c(2)])) # c(2085:2094)
df_tech_covs <- df[, c(1, cols_tech_covs)]
df_tech_covs <- left_join(df_tech_covs, Methyl_merged, by = "IDC") # methylation data
df_tech_covs <- left_join(df_tech_covs, df_ethv3, by = "IDC") # ethnicity data
df_gen <- df_gen[, c("IDC", "GESTBIR", "GENDER", "WEIGHT")] # add sex, GA, and weight to tech covariate df
df_tech_covs <- left_join(df_tech_covs, df_gen, by = "IDC")

# construct standardized variables for 
df_tech_covs$C1_scaled <- scale(df_tech_covs$GWAv3C1, center = TRUE, scale = TRUE)
df_tech_covs$C2_scaled <- scale(df_tech_covs$GWAv3C2, center = TRUE, scale = TRUE)

# set plate and gender as factors
#df_tech_covs$Plaat <- df_tech_covs$Plaat %>% as.factor()
#df_tech_covs$GENDER <- df_tech_covs$GENDER %>% as.factor()

# drop covariate_cols from df
extracted_counts <- df[, -c(cols_tech_covs)]
#[, -covariate_cols]

# set Sample.Name as rowname and remove column from extracted_counts
extracted_counts <- extracted_counts %>% t()
colnames_extracted_counts <- extracted_counts[1, ]
colnames(extracted_counts) <- extracted_counts[1, ]
extracted_counts <- extracted_counts[-1, ]
counts <- data.frame(apply(extracted_counts, 2, function(x) as.numeric(as.character(x))))
rownames(counts) <- rownames(extracted_counts)
colnames(counts) <- colnames_extracted_counts # code line 127 does something weird with the colnames (puts X. in there), fix

# Rename to raw_counts for merging
raw_counts <- counts

# merge with IDC column
cleaned_counts <- cbind.data.frame(colnames(raw_counts), t(raw_counts))
colnames(cleaned_counts)[1] <- "Sample.Name"
cleaned_counts <- left_join(cleaned_counts, original_data[, c(1:2)], by = "Sample.Name")[, c(2085, 1:2084)]

# save cleaned raw counts for archive
write.xlsx(cleaned_counts, file.path(path, paste0("data_clean/archived/miRNA_GenR_cleaned_", format(Sys.Date(), "%d_%m_%Y"), ".xlsx")), 
           colNames = TRUE, rowNames = FALSE,
           overwrite = FALSE)

# remove surplus datasets, plates etc
rm(list = ls(pattern = "^Plate\\d+(_Concentration|_Haemolysis)?$"), 
   table, dataframes, concentration, concentration_file, 
   haemolysis, haemolysis_file, i, n_plates, name, 
   name_concentration, name_file, name_haemolysis, 
   df_QCres, counts, df, extracted_counts, 
   colnames_extracted_counts, covariate_cols, cols_tech_covs,
   Methyl450k3, Methyl450k3_Salas, MethylEPIC1, 
   GENR_EPICv1METH_birth_CellTypes_combined.data, 
   Methyl_450k_merged, Methyl_EPIC_merged, Methyl_merged,
   df_ethv3, df_gen)
# ==================================================================================================
# TMM-normalization
# ==================================================================================================
### Creating a DGE-list object & examining miRNA expression level filters
dge <- DGEList(counts = raw_counts)

# Filter out lowly expressed miRNAs using edgeR's `filterByExpr` function
# base settings of function are: min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7
keep <- filterByExpr(dge,
                     min.count = 10, 
                     min.total.count = 15, 
                     large.n = 10, 
                     min.prop = 0.7)
table(keep) # filters only 1 miRNA, miR.513b.5p 
# base settings of function are: min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7
keep2 <- filterByExpr(dge,
                      min.count = 15, 
                      min.total.count = 15, 
                      large.n = 10, 
                      min.prop = 0.7)
table(keep2) # 57 miRNAs filtered if min.count upped from 10 to 15 counts
dge_filtered <- dge[keep, keep.lib.sizes = FALSE] # Filter on expression level miRNAs
dge_filtered2 <- dge[keep2, keep.lib.sizes = FALSE] # Filter on expression level miRNAs

plot_grid(
  # Plot 1A: Histogram of total counts per sample before filtering to see initial expression distribution
  ggplot(data.frame(total_counts = colSums(raw_counts)), aes(x = total_counts)) +
    geom_histogram(bins = 150, color = "black", fill = "lightblue") +
    ggtitle("Distribution of Total Counts per Sample",
            subtitle = "Before Filtering") +
    xlab("Total Counts per Sample") +
    ylab("Frequency") +
    scale_x_continuous(limits = c(0, 1e7)) + 
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ),
  # Plot 1B: Histogram of total counts per microRNA before filtering to see initial expression distribution
  ggplot(data.frame(total_counts = rowSums(raw_counts)), aes(x = total_counts)) +
    geom_histogram(bins = 100, color = "black", fill = "lightblue") +
    ggtitle("Distribution of Total Counts per MicroRNA",
            subtitle = "Before Filtering") +
    xlab("Total Counts per MicroRNA") +
    ylab("Frequency") +
    scale_x_continuous(limits = c(0, 3e6)) + 
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
      #size = 22),
      #axis.title = element_text(size = 15),
      #legend.title = element_text(size = 15),
      #legend.text = element_text(size = 12)
    ),
  # Plot 2A: Histogram of total counts per sample after filtering to check effect of `filterByExpr`
  ggplot(data.frame(total_counts = colSums(dge_filtered$counts)), aes(x = total_counts)) +
    geom_histogram(bins = 150, color = "black", fill = "lightgreen") +
    ggtitle("Distribution of Total Counts per Sample",
            subtitle = "After Default Filtering") +
    xlab("Total Counts per Sample") +
    ylab("Frequency") +
    scale_x_continuous(limits = c(0, 1e7)) + 
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ),
  # Plot 2B: Histogram of total counts per microRNA after filtering to check effect of `filterByExpr`
  ggplot(data.frame(total_counts = rowSums(dge_filtered$counts)), aes(x = total_counts)) +
    geom_histogram(bins = 100, color = "black", fill = "lightgreen") +
    ggtitle("Distribution of Total Counts per MicroRNA",
            subtitle = "After Default Filtering") +
    xlab("Total Counts per MicroRNA") +
    ylab("Frequency") +
    scale_x_continuous(limits = c(0, 3e6)) + 
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ),
  # Plot 3A: Histogram of total counts per sample after filtering to check effect of `filterByExpr`
  ggplot(data.frame(total_counts = colSums(dge_filtered2$counts)), aes(x = total_counts)) +
    geom_histogram(bins = 150, color = "black", fill = "orange") +
    ggtitle("Distribution of Total Counts per Sample",
            subtitle = "After Stricter Filtering") +
    xlab("Total Counts per Sample") +
    ylab("Frequency") +
    scale_x_continuous(limits = c(0, 1e7)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ),
  # Plot 3B: Histogram of total counts per microRNA after filtering to check effect of `filterByExpr`
  ggplot(data.frame(total_counts = rowSums(dge_filtered2$counts)), aes(x = total_counts)) +
    geom_histogram(bins = 100, color = "black", fill = "orange") +
    ggtitle("Distribution of Total Counts per MicroRNA",
            subtitle = "After Stricter Filtering") +
    xlab("Total Counts per MicroRNA") +
    ylab("Frequency") +
    scale_x_continuous(limits = c(0, 3e6)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ),
  ncol = 2, axis = "b",
  byrow = TRUE, labels = c("1A", "1B", "2A", "2B", "3A", "3B"),
  rel_widths = c(3,2)
)

### Normalization using TMM & Log-transformation
dge <- calcNormFactors(dge, method = "TMM")

TMM_CPM_orig <- cpm(dge, log = FALSE) %>% t(.) %>% as.data.frame()
TMM_CPM_orig <- cbind.data.frame(rownames(TMM_CPM_orig), TMM_CPM_orig)
colnames(TMM_CPM_orig)[1] <- "Sample.Name"
TMM_CPM_orig <- left_join(TMM_CPM_orig, original_data[, c(1:2)], by = "Sample.Name")[, c(2085, 1:2084)]

logCPM_orig <- cpm(dge, log = TRUE) %>% t(.) %>% as.data.frame()
logCPM_orig <- cbind.data.frame(rownames(logCPM_orig), logCPM_orig)
colnames(logCPM_orig)[1] <- "Sample.Name"
logCPM_orig <- left_join(logCPM_orig, original_data[, c(1:2)], by = "Sample.Name")[, c(2085, 1:2084)]

# save data
write.xlsx(TMM_CPM_orig, file.path(path, paste0("data_clean/Midas/miRNA_GenR_counts_TMM_", format(Sys.Date(), "%d_%m_%Y"), ".xlsx")), 
           colNames = TRUE, rowNames = FALSE,
           overwrite = FALSE)

write.xlsx(logCPM_orig, file.path(path, paste0("data_clean/Midas/miRNA_GenR_counts_TMM_logCPM_", format(Sys.Date(), "%d_%m_%Y"), ".xlsx")), 
           colNames = TRUE, rowNames = FALSE,
           overwrite = FALSE)

rownames(TMM_CPM_orig) <- TMM_CPM_orig[, 2]
TMM_CPM_orig <- TMM_CPM_orig[, -c(1, 2)]

rownames(logCPM_orig) <- logCPM_orig[, 2]
logCPM_orig <- logCPM_orig[, -c(1, 2)]

### Winsorization to address outliers
#### Applying winsorization
# According to recommendations by https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03230-w

# Calculate TMM-normalized CPM that are not logtransformed for winsorization
CPM <- cpm(dge, log = FALSE)

# quantiles before TMM-normalization
quantile(dge$counts, probs = c(0.025, 0.975), na.rm = TRUE) # 7  - 9196

# 11 samples with extremely high concentrations (> 1000 uL) at plate 7, 11 / 1695 * 100 = ~0.65%
# However 2 are not included in the sample due to QC status, corresponds to ~0.53%
# 8 samples on plate 7 have also high but not as extreme as the others (100 - 700 uL), might be problematic but not sure
# If so, we'd have 17 samples on plate 7 that are problematic, which corresponds to ~1% of samples
CPM_winsorized <- t(apply(CPM, 1, winsorize, probs = c(0.025, 0.975), verbose = FALSE))

quantile(CPM_winsorized, probs = c(0.025, 0.975), na.rm = TRUE) 

# Examining first 6 miRNAs before and after winsorization
plot(log(CPM[1, ]))
plot(log(CPM_winsorized[1, ]))

plot((CPM[1, ]), (CPM_winsorized[1, ]))
plot((CPM[1, ]), df_tech_covs$Plaat)
plot((CPM_winsorized[1, ]), df_tech_covs$Plaat)
plot(log(CPM[1, ]), df_tech_covs$Plaat)
plot(log(CPM_winsorized[1, ]), df_tech_covs$Plaat)

plot((CPM[2, ]), (CPM_winsorized[2, ]))
plot((CPM[2, ]), df_tech_covs$Plaat)
plot((CPM_winsorized[2, ]), df_tech_covs$Plaat)
plot(log(CPM[2, ]), df_tech_covs$Plaat)
plot(log(CPM_winsorized[2, ]), df_tech_covs$Plaat)

plot((CPM[3, ]), (CPM_winsorized[3, ]))
plot((CPM[3, ]), df_tech_covs$Plaat)
plot((CPM_winsorized[3, ]), df_tech_covs$Plaat)
plot(log(CPM[3, ]), df_tech_covs$Plaat)
plot(log(CPM_winsorized[3, ]), df_tech_covs$Plaat)

plot((CPM[4, ]), (CPM_winsorized[4, ]))
plot((CPM[4, ]), df_tech_covs$Plaat)
plot((CPM_winsorized[4, ]), df_tech_covs$Plaat)
plot(log(CPM[4, ]), df_tech_covs$Plaat)
plot(log(CPM_winsorized[4, ]), df_tech_covs$Plaat)

plot((CPM[5, ]), (CPM_winsorized[5, ]))
plot((CPM[5, ]), df_tech_covs$Plaat)
plot((CPM_winsorized[5, ]), df_tech_covs$Plaat)
plot(log(CPM[5, ]), df_tech_covs$Plaat)
plot(log(CPM_winsorized[5, ]), df_tech_covs$Plaat)

plot((CPM[6, ]), (CPM_winsorized[6, ]))
plot((CPM[6, ]), df_tech_covs$Plaat)
plot((CPM_winsorized[6, ]), df_tech_covs$Plaat)
plot(log(CPM[6, ]), df_tech_covs$Plaat)
plot(log(CPM_winsorized[6, ]), df_tech_covs$Plaat)

#### Comparison unnormalized-winsorized dataset with original dataset and the renormalization & logtransformation of winsorized data
# Compute effective library sizes
effective_lib_sizes <- dge$samples$lib.size * dge$samples$norm.factors 

# Recover raw counts from CPM & rounding
raw_counts_winsor <- round(sweep(CPM_winsorized, 2, effective_lib_sizes, "*") / 1e6)

# Ensure no negative values (edge case rounding issue)
raw_counts_winsor[raw_counts_winsor < 0] <- 0  

# Extract original raw counts
original_counts <- dge$counts  

# Summary statistics
# Check correlation between original and reconstructed counts, 0.985031
cor_results <- cor(as.vector(original_counts), as.vector(raw_counts_winsor), method = "pearson")

as.matrix(table(original_counts[c(1:2083), c(1:1695)] == raw_counts_winsor[c(1:2083), c(1:1695)]))

# Descriptives original counts first 10 microRNAs
desc_orig_counts <- cbind.data.frame(
  as.data.frame(colnames(TMM_CPM_orig)[1:10]),
  psych::describe(TMM_CPM_orig)[c(1:10), c(3, 4, 13, 5, 8:9, 11:12)]  %>%
    round(., 2),
  sapply(TMM_CPM_orig[, c(1:10)], quantile, probs = c(.25, .75)) %>% 
    as.data.frame(.) %>% 
    round(., 2) %>% 
    t(.)
) %>% select(., c(1:5, 10, 11, 6:9))

colnames(desc_orig_counts)[1] <- "MicroRNA"

desc_orig_counts


# Log TMM-normaize and logtransform winsorized data again
dge_winsor <- DGEList(counts = raw_counts_winsor)
dge_winsor <- calcNormFactors(dge_winsor, method = "TMM")

TMM_CPM_winsor <- cpm(dge_winsor, log = FALSE) %>% t(.) %>% as.data.frame()
logCPM_winsor <- cpm(dge_winsor, log = TRUE) %>% t(.) %>% as.data.frame()

# Descriptives winsorized counts first 10 microRNAs
desc_winsor_counts <- cbind.data.frame(
  as.data.frame(colnames(TMM_CPM_winsor)[1:10]),
  psych::describe(TMM_CPM_winsor)[c(1:10), c(3, 4, 13, 5, 8:9, 11:12)]  %>%
    round(., 2),
  sapply(TMM_CPM_winsor[, c(1:10)], quantile, probs = c(.25, .75)) %>% 
    as.data.frame(.) %>% 
    round(., 2) %>% 
    t(.)
) %>% select(., c(1:5, 10, 11, 6:9))

colnames(desc_winsor_counts)[1] <- "MicroRNA"

desc_winsor_counts 

# Check correlation between original and winsorized logtransformed counts
cor_results_log <- cor(unlist(as.vector(logCPM_orig)), unlist(as.vector(logCPM_winsor)), method = "pearson")
cor_results_log

# scatterplots to compare original and winsorized logCPM values across samples, including visualization against plate variables
plot(logCPM_orig[, 1], logCPM_winsor[, 1])
plot(logCPM_orig[, 1], df_tech_covs$Plaat)
plot(logCPM_winsor[, 1], df_tech_covs$Plaat)

plot(logCPM_orig[, 2], logCPM_winsor[, 2])
plot(logCPM_orig[, 2], df_tech_covs$Plaat)
plot(logCPM_winsor[, 2], df_tech_covs$Plaat)

plot(logCPM_orig[, 3], logCPM_winsor[, 3])
plot(logCPM_orig[, 3], df_tech_covs$Plaat)
plot(logCPM_winsor[, 3], df_tech_covs$Plaat)

plot(logCPM_orig[, 4], logCPM_winsor[, 4])
plot(logCPM_orig[, 4], df_tech_covs$Plaat)
plot(logCPM_winsor[, 4], df_tech_covs$Plaat)

plot(logCPM_orig[, 5], logCPM_winsor[, 5])
plot(logCPM_orig[, 5], df_tech_covs$Plaat)
plot(logCPM_winsor[, 5], df_tech_covs$Plaat)

plot(logCPM_orig[, 6], logCPM_winsor[, 6])
plot(logCPM_orig[, 6], df_tech_covs$Plaat)
plot(logCPM_winsor[, 6], df_tech_covs$Plaat)

### Saving winsorized datasets
TMM_CPM_winsor_save <- TMM_CPM_winsor
TMM_CPM_winsor_save <- cbind.data.frame(rownames(TMM_CPM_winsor_save), TMM_CPM_winsor_save)
colnames(TMM_CPM_winsor_save)[1] <- "Sample.Name"
TMM_CPM_winsor_save <- left_join(TMM_CPM_winsor_save, original_data[, c(1:2)], by = "Sample.Name")[, c(2085, 1:2084)]

logCPM_winsor_save <- logCPM_winsor
logCPM_winsor_save <- cbind.data.frame(rownames(logCPM_winsor_save), logCPM_winsor_save)
colnames(logCPM_winsor_save)[1] <- "Sample.Name"
logCPM_winsor_save <- left_join(logCPM_winsor_save, original_data[, c(1:2)], by = "Sample.Name")[, c(2085, 1:2084)]

# saving datasets
write.xlsx(TMM_CPM_winsor_save, file.path(path, paste0("data_clean/Midas/miRNA_GenR_counts_TMM_WIN_", format(Sys.Date(), "%d_%m_%Y"), ".xlsx")), 
           colNames = TRUE, rowNames = FALSE,
           overwrite = FALSE)

write.xlsx(logCPM_winsor_save, file.path(path, paste0("data_clean/Midas/miRNA_GenR_counts_TMM_WIN_logCPM_", format(Sys.Date(), "%d_%m_%Y"), ".xlsx")), 
           colNames = TRUE, rowNames = FALSE,
           overwrite = FALSE)

rm(TMM_CPM_winsor_save, logCPM_winsor_save)

# ==================================================================================================
# Plate 7 issues addressed by winsorization
# ==================================================================================================
# MiRNA expression distribution before (A) and after winsorization (B) for let-7a-2-3p (1), let-7a-30 (2), and let-7a-5p
data_plotjes <- cbind.data.frame(logCPM_orig, df_tech_covs)
data_plotjes$concentration_high <- ifelse(data_plotjes$Concentration >= 10000, "high", "low")

data_plotjes_WIN <- cbind.data.frame(logCPM_winsor, df_tech_covs)
data_plotjes_WIN$concentration_high <- ifelse(data_plotjes_WIN$Concentration >= 10000, "high", "low")

# "let.7a.2.3p" "let.7a.3p"   "let.7a.5p" 
miRNA1_orig <- data_plotjes %>%
  ggplot(aes(x = let.7a.2.3p, y = Plaat, colour = concentration_high)) +
  geom_jitter(width = 0.1, alpha = 0.6, show.legend = FALSE) +  
  theme_minimal() +
  scale_color_manual(values = c("low" = "grey50", "high" = "red")) +  # Adjust colors as needed
  labs(title = "Expression let-7a-2-3p per plate",
       tag = "1A",
       x = "Expression let-7a-2-3p",
       y = "Plate",
       colour = "MiRNA concentration") +
  theme(plot.title = element_text(hjust = 0.5))

miRNA2_orig <- data_plotjes %>%
  ggplot(aes(x = let.7a.3p, y = Plaat, colour = concentration_high)) +
  geom_jitter(width = 0.1, alpha = 0.6, show.legend = FALSE) +  
  theme_minimal() +
  scale_color_manual(values = c("low" = "grey50", "high" = "red")) +  # Adjust colors as needed
  labs(title = "Expression let-7a-3p per plate",
       tag = "2A",
       x = "Expression let-7a-3p",
       y = "Plate",
       colour = "MiRNA concentration") +
  theme(plot.title = element_text(hjust = 0.5))

miRNA3_orig <- data_plotjes %>%
  ggplot(aes(x = let.7a.5p, y = Plaat, colour = concentration_high)) +
  geom_jitter(width = 0.1, alpha = 0.6, show.legend = FALSE) +  
  theme_minimal() +
  scale_color_manual(values = c("low" = "grey50", "high" = "red")) +  # Adjust colors as needed
  labs(title = "Expression let-7a-5p per plate",
       tag = "3A",
       x = "Expression let-7a-5p",
       y = "Plate",
       colour = "MiRNA concentration") +
  theme(plot.title = element_text(hjust = 0.5))

# data_plotjes_WIN "let.7a.2.3p" "let.7a.3p"   "let.7a.5p" 
miRNA1_WIN <- data_plotjes_WIN %>%
  ggplot(aes(x = let.7a.2.3p, y = Plaat, colour = concentration_high)) +
  geom_jitter(width = 0.1, alpha = 0.6, show.legend = FALSE) +  
  theme_minimal() +
  scale_color_manual(values = c("low" = "grey50", "high" = "red")) +  # Adjust colors as needed
  labs(title = "Winsorized expression let-7a-2-3p per plate",
       tag = "1B",
       x = "Expression let-7a-2-3p",
       y = "Plate",
       colour = "MiRNA concentration") +
  theme(plot.title = element_text(hjust = 0.5))

miRNA2_WIN <- data_plotjes_WIN %>%
  ggplot(aes(x = let.7a.3p, y = Plaat, colour = concentration_high)) +
  geom_jitter(width = 0.1, alpha = 0.6, show.legend = FALSE) +  
  theme_minimal() +
  scale_color_manual(values = c("low" = "grey50", "high" = "red")) +  # Adjust colors as needed
  labs(title = "Winsorized expression let-7a-3p per plate",
       tag = "2B",
       x = "Expression let-7a-3p",
       y = "Plate",
       colour = "MiRNA concentration") +
  theme(plot.title = element_text(hjust = 0.5))

miRNA3_WIN <- data_plotjes_WIN %>%
  ggplot(aes(x = let.7a.5p, y = Plaat, colour = concentration_high)) +
  geom_jitter(width = 0.1, alpha = 0.6, show.legend = FALSE) +  
  theme_minimal() +
  scale_color_manual(values = c("low" = "grey50", "high" = "red")) +  # Adjust colors as needed
  labs(title = "Winsorized expression let-7a-5p per plate",
       tag = "3B",
       x = "Expression let-7a-5p",
       y = "Plate",
       colour = "MiRNA concentration") +
  theme(plot.title = element_text(hjust = 0.5))

plot(miRNA1_orig)
plot(miRNA1_WIN)

plot(miRNA2_orig)
plot(miRNA2_WIN)

plot(miRNA3_orig)
plot(miRNA3_WIN)


# ==================================================================================================
# Sanity checks on the data
# ==================================================================================================
## Descriptive summary statistics on miRNAs
# Descriptives original (non-winsorized) TMM-normalized counts (per million) for the first 10 miRNAs in the dataset
desc_orig_counts
# Descriptives winsorized & TMM-normalized counts (per million) for the first 10 miRNAs in the dataset
desc_winsor_counts

## Constructing summary statistics}
# Descriptives original counts all microRNAs
desc_orig_counts_full <- cbind.data.frame(
  as.data.frame(colnames(TMM_CPM_orig)),
  psych::describe(TMM_CPM_orig)[, c(3, 4, 13, 5, 8:12)]  %>%
    round(., 2),
  sapply(TMM_CPM_orig, quantile, probs = c(.25, .75)) %>% 
    as.data.frame(.) %>% 
    round(., 2) %>% 
    t(.)
) %>% select(., c(1:5, 11, 12, 6:10))

colnames(desc_orig_counts_full)[1] <- "MicroRNA"
desc_orig_counts_full[, 1] <- gsub(".", "-", desc_orig_counts_full[, 1], fixed = TRUE)

# Descriptives winsorized counts all microRNAs
desc_winsor_counts_full <- cbind.data.frame(
  as.data.frame(colnames(TMM_CPM_winsor)),
  psych::describe(TMM_CPM_winsor)[, c(3, 4, 13, 5, 8:12)]  %>%
    round(., 2),
  sapply(TMM_CPM_winsor, quantile, probs = c(.25, .75)) %>% 
    as.data.frame(.) %>% 
    round(., 2) %>% 
    t(.)
) %>% select(., c(1:5, 11, 12, 6:10))
colnames(desc_winsor_counts_full)[1] <- "MicroRNA"
desc_winsor_counts_full[, 1] <- gsub(".", "-", desc_winsor_counts_full[, 1], fixed = TRUE)


miRNA_descriptives <- list(desc_orig_counts_full, 
                           desc_winsor_counts_full)

# save data
write.xlsx(miRNA_descriptives, file.path(path, paste0("data_clean/DataWiki/miRNA_GenR_descriptives_", format(Sys.Date(), "%d_%m_%Y"), ".xlsx")), 
           sheetName = c("counts_TMM", 
                         "counts_TMM_WIN"),
           colNames = TRUE, rowNames = FALSE,
           overwrite = TRUE)

## Density plots for normalized data

plotDensities(logCPM_orig[, c(1:20)], 
              main = "Density Plot of Log-CPM After Normalization",
              legend = "topright")

plotDensities(logCPM_winsor[, c(1:20)], 
              main = "Density Plot of Log-CPM After Normalization and Winsorization",
              legend = "topright")

## Principal component analysis
# Perform PCA for all datasets
# original and winsorized dataset
pca_orig <- prcomp(logCPM_orig, center = TRUE, scale. = TRUE)
PCs_orig <- pca_orig$x[, c(1:10)] %>% as.data.frame()
PCs_orig$Sample.Name <- rownames(PCs_orig)

pca_winsor <- prcomp(logCPM_winsor, center = TRUE, scale. = TRUE)
PCs_winsor <- pca_winsor$x[, c(1:10)] %>% as.data.frame()
PCs_winsor$Sample.Name <- rownames(PCs_winsor)



pca_orig_expl_var <- pca_orig$sdev^2  # Eigenvalues (variance for each PC)
pca_orig_prop_variance <- pca_orig_expl_var / sum(pca_orig_expl_var)  # Normalize

pca_winsor_expl_var <- pca_winsor$sdev^2  # Eigenvalues (variance for each PC)
pca_winsor_prop_variance <- pca_winsor_expl_var / sum(pca_winsor_expl_var)  # Normalize

PCs <- c("PC 1", "PC 2", "PC 3", "PC 4", "PC 5",
         "PC 6", "PC 7", "PC 8", "PC 9", "PC 10")
PCs_orig_propvar <- pca_orig_prop_variance[1:10] %>% round(., 2)
PCs_winsor_propvar <- pca_winsor_prop_variance[1:10] %>% round(., 2)

table_PCs <- cbind.data.frame(PCs, PCs_orig_propvar, PCs_winsor_propvar)
colnames(table_PCs) <- c("PCs", "PCs original data", "PCs winsorized data")

table_PCs 

### Multidimensional Scaling (MDS) plots for sample similarity
# MDS Plot for visualization of sample similarity, original dataset
average_expression <- rowMeans(logCPM_orig)
pca_orig_df <- data.frame(PC1 = pca_orig$x[, 1], 
                          PC2 = pca_orig$x[, 2], 
                          AvgExpr = average_expression, 
                          Sample = rownames(logCPM_orig))

explained_variance <- round(100 * (pca_orig$sdev^2 / sum(pca_orig$sdev^2)), 2)

x_label <- paste0("Principal Component 1 (", explained_variance[1], "% variance explained)")
y_label <- paste0("Principal Component 2 (", explained_variance[2], "% variance explained)")

fig_PC_orig <- ggplot(pca_orig_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = AvgExpr), size = 3) +
  scale_color_gradient(low = "blue", high = "red") +
  ggtitle("PCA Plot of Samples with Average Expression as Color, original logCPM") +
  xlab(x_label) +
  ylab(y_label) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color = "Avg Expression")

# MDS Plot for visualization of sample similarity, winsorized dataset
average_expression <- rowMeans(logCPM_winsor)
pca_winsor_df <- data.frame(PC1 = pca_winsor$x[, 1], 
                            PC2 = pca_winsor$x[, 2], 
                            AvgExpr = average_expression, 
                            Sample = rownames(logCPM_winsor))

explained_variance <- round(100 * (pca_winsor$sdev^2 / sum(pca_winsor$sdev^2)), 2)

x_label <- paste0("Principal Component 1 (", explained_variance[1], "% variance explained)")
y_label <- paste0("Principal Component 2 (", explained_variance[2], "% variance explained)")

fig_PC_winsor <- ggplot(pca_winsor_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = AvgExpr), size = 3) +
  scale_color_gradient(low = "blue", high = "red") +
  ggtitle("PCA Plot of Samples with Average Expression as Color, winsorized logCPM") +
  xlab(x_label) +
  ylab(y_label) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color = "Avg Expression")

plot_grid(fig_PC_orig,
          fig_PC_winsor,
          ncol = 1)#, #axis = "b",
# byrow = TRUE, labels = c("1A", "1B", "2A", "2B", "3A", "3B")
#)
# ==================================================================================================
# Covariates
# ==================================================================================================
### Descriptives
# Loop through plate numbers and create dichotomous variables
for (i in 1:18) {
  df_tech_covs[[paste0("Plate", i)]] <- ifelse(df_tech_covs$Plaat == i, 1, 0)
}

df_tech_covs$Preterm_birth <- ifelse(df_tech_covs$GESTBIR <= 37, 1, 0)
df_tech_covs$Postterm_birth <- ifelse(df_tech_covs$GESTBIR >= 42, 1, 0)

colnames(df_tech_covs)[c(23:25)] <- c("Gestational Age", "Gender", "Birth Weight")

# descriptives table
desc_covs <- cbind.data.frame(
  as.data.frame(colnames(df_tech_covs)[c(3, 4, 13:19, 21:22, 23, 25)]),
  psych::describe(df_tech_covs)[c(3, 4, 13:19, 21:22, 23, 25),  c(2, 3, 4, 13, 5, 8:12)]  %>%
    round(., 2) )


colnames(desc_covs)[1] <- "Variable"
desc_covs[2, 2:11] <- unlist(desc_covs[2, 2:11]) %>% round(., 0) %>% as.factor(.)


# number of males
n_males <- df_tech_covs %>%
  filter(., Gender == 1) %>%
  nrow(.)
# number of preterms
n_preterms <- df_tech_covs %>%
  filter(., Preterm_birth == 1) %>%
  nrow(.)
# number of postterms
n_postterms <- df_tech_covs %>%
  filter(., Postterm_birth == 1) %>%
  nrow(.)

desc_covs

### Correlation between covariates and miRNA principal components
# Correlations with PCs
make_cor_matrix(df_tech_covs[, c(3, 4, 13:19, 21:22, 26:27, 23:25, 34, 45, 46, 47)], PCs_orig[, c(1:10)],
                "Correlation matrix Covariates & PCs Original data")

make_cor_matrix(df_tech_covs[, c(3, 4, 13:19, 21:22, 26:27, 23:25, 34, 45, 46, 47)], PCs_winsor[, c(1:10)], 
                "Correlation matrix Covariates & PCs Winsorized data")

# ==================================================================================================
# Haemolysis check
# ==================================================================================================
haemo_miRNAs <- c("miR.451a", "miR.324.5p", "miR.425.5p", "miR.191.5p", "miR.23a.3p", "miR.150.5p", "miR.29a.3p") 
#haemo_miRNAs %in% colnames(logCPM_orig)

haemo_miRNAs_subset <- logCPM_orig[, c(haemo_miRNAs)] 
haemo_miRNAs_subset$Sample.Name <- rownames(haemo_miRNAs_subset)
haemo_miRNAs_subset <- haemo_miRNAs_subset %>%
  left_join(., df_tech_covs, by = "Sample.Name")

miRNA_haemo <- haemo_miRNAs_subset$miR.451a
miRNA <- haemo_miRNAs_subset$miR.150.5p # miR.451a


haemo_miRNAs_subset$miR23a.miR451.ratio <- haemo_miRNAs_subset$miR.23a.3p / miRNA_haemo
haemo_miRNAs_subset$miR23a.miR451.delta <- haemo_miRNAs_subset$miR.23a.3p - miRNA_haemo
haemo_miRNAs_subset$miR23a.miR451.delta.seven <- ifelse(abs(haemo_miRNAs_subset$miR23a.miR451.delta) >= 7, 1, 0)
haemo_miRNAs_subset$miR23a.miR451.delta.five <- ifelse(abs(haemo_miRNAs_subset$miR23a.miR451.delta) < 5, 1, 0)

# Haemolysis plots
haemo1 <- haemo_miRNAs_subset %>%
  ggplot(., mapping = aes(x = miRNA_haemo, y = miR.23a.3p, colour = as.factor(Haemolysis))) +
  geom_point(aes(alpha = as.factor(Haemolysis))) +  
  theme_minimal() +
  scale_color_manual(values = c("0" = "gray", "1" = "yellow", "2" = "green", 
                                "3" = "blue", "9" = "red")) +
  scale_alpha_manual(values = c("0" = 0.2, "1" = 0.8, "2" = 1, 
                                "3" = 1, "9" = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Haemolysis in miR-451a vs. miR-23a-3p",
       x = "MiR-451a",
       y = "MiR-23a-3p",
       colour = "Haemolysis",
       alpha = "Haemolysis")

haemo2 <- haemo_miRNAs_subset %>%
  ggplot(., mapping = aes(x = miRNA_haemo, y = miR.23a.3p, colour = as.factor(miR23a.miR451.delta.seven))) +
  geom_point(alpha = 0.6) +  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Absolute miR-Delta in miR-451a vs. miR-23a-3p, divided by |delta| cut-off >= 7",
       x = "MiR-451a",
       y = "MiR-23a-3p",
       colour = "|Delta| >= 7")

haemo3 <- haemo_miRNAs_subset %>%
  ggplot(aes(x = miRNA_haemo, y = miR.23a.3p, colour = miR23a.miR451.ratio)) +
  geom_point(alpha = 0.6) +  
  theme_minimal() +
  scale_color_gradient(low = "blue", high = "red") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "MiR-Ratio in miR-451a vs. miR-23a-3p",
       x = "MiR-451a",
       y = "MiR-23a-3p",
       colour = "MiR-Ratio")

# Weird plate shit
haemo4 <- haemo_miRNAs_subset %>%
  ggplot(aes(x = miRNA_haemo, y = Plaat, colour = as.factor(Haemolysis))) +
  geom_jitter(width = 0.1, aes(alpha = as.factor(Haemolysis))) +
  theme_minimal() +
  scale_color_manual(values = c("0" = "gray", "1" = "yellow", "2" = "green", 
                                "3" = "blue", "9" = "red")) +
  scale_alpha_manual(values = c("0" = 0.3, "1" = 0.8, "2" = 1, 
                                "3" = 1, "9" = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Haemolysis vs. Plate",
       x = "Expression levels of miR-451a",
       y = "Plate",
       colour = "Haemolysis",
       alpha = "Haemolysis")
haemo4

plot_grid(haemo1, haemo2, haemo3,
          ncol = 1, axis = "b",
          byrow = TRUE, labels = c("1", "2", "3")
          )

# table ratiodelta

a <- aggregate(miR23a.miR451.ratio ~ Haemolysis, data = haemo_miRNAs_subset, median) # ratio is lower for samples with any haemolysis, not linear
b <- aggregate(miR23a.miR451.delta ~ Haemolysis, data = haemo_miRNAs_subset, median) # delta is higher for samples with any haemolysis, not linear

c <- left_join(a, b, by = "Haemolysis")
colnames(c)[c(2, 3)] <- c("miR-ratio", "miR-delta")

c
