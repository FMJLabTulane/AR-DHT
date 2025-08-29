combined <- qread(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined_with_MQI_all_beta.qs)")
# Example: AR gene set (you can expand/replace with your curated list)
# Example: alternate AR gene set (stress, metabolism, mitochondria relevant)
ar_genes <- c("FKBP5","NDRG1","HSP90AA1","UBE2C",
              "CCND1","MYC","IGF1R","IRS2","PRKAA2")

# Make sure gene symbols are present
ar_genes <- intersect(ar_genes, rownames(combined))

# Add AR activation score
combined <- AddModuleScore(
  combined,
  features = list(ar_genes),
  name = "ARscore",
  nbin = 24,
  ctrl = 100,
  seed = 42
)

# The score will be stored in obj@meta.data$ARscore1
head(combined@meta.data$ARscore1)


# --- AR activation score: DHT[10nM] vs EtOH ---

library(Seurat)
library(dplyr)
library(ggplot2)

# 1) Pull per-cell data
df <- FetchData(
  combined,
  vars = c("ARscore1", "treatment"),
  slot = "data"
) %>%
  tibble::as_tibble(rownames = "cell") %>%
  dplyr::filter(!is.na(ARscore1),
                !is.na(treatment),
                treatment %in% c("DHT[10nM]", "EtOH")) %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("EtOH", "DHT[10nM]")))

# 2) Basic stats (Wilcoxon, group sizes, FC of means)
n_etoh <- sum(df$treatment == "EtOH")
n_dht  <- sum(df$treatment == "DHT[10nM]")

p_wilcox <- wilcox.test(ARscore1 ~ treatment, data = df)$p.value

mean_etoh <- mean(df$ARscore1[df$treatment == "EtOH"], na.rm = TRUE)
mean_dht  <- mean(df$ARscore1[df$treatment == "DHT[10nM]"], na.rm = TRUE)
fc_mean   <- if (is.finite(mean_etoh) && mean_etoh != 0) mean_dht / mean_etoh else NA_real_

label_txt <- sprintf(
  "Wilcoxon p = %.3g\nMean FC (DHT/EtOH) = %.3f\nn: EtOH=%d, DHT=%d",
  p_wilcox, fc_mean, n_etoh, n_dht
)

# 3) Plot
p <- ggplot(df, aes(x = treatment, y = ARscore1)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.25, size = 0.6) +
  labs(x = NULL, y = "AR activation score (module: ARscore1)",
       title = "AR activation score: DHT[10nM] vs EtOH") +
  annotate("text", x = 1.5, y = max(df$ARscore1, na.rm = TRUE), 
           label = label_txt, hjust = 0.5, vjust = -0.4, size = 3.4) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p)
label_txt

library(Seurat)
library(dplyr)
library(ggplot2)

# --- 0) pick the AR symbol that exists in your object (AR or NR3C4) ---
ar_symbol <- intersect(c("AR","NR3C4"), rownames(combined))
if (length(ar_symbol) == 0) stop("AR/NR3C4 not found in rownames(combined).")
ar_symbol <- ar_symbol[1]

# --- 1) flag AR+ cells using the normalized 'data' slot of the default assay ---
assay_use <- DefaultAssay(combined)
ar_expr <- GetAssayData(combined, assay = assay_use, slot = "data")[ar_symbol, ]
ARpos <- ar_expr > 0

message(sprintf("AR+ cells: %d / %d (%.1f%%)",
                sum(ARpos), length(ARpos), 100*mean(ARpos)))

# subset to AR+ only
combined_ARpos <- subset(combined, cells = colnames(combined)[ARpos])

# --- 2) pull ARscore1 and treatment just for AR+ cells ---
df <- FetchData(
  combined_ARpos,
  vars = c("ARscore1", "treatment")
) %>%
  tibble::as_tibble(rownames = "cell") %>%
  dplyr::filter(!is.na(ARscore1),
                !is.na(treatment),
                treatment %in% c("EtOH", "DHT[10nM]")) %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("EtOH", "DHT[10nM]")))

# basic stats
n_etoh <- sum(df$treatment == "EtOH")
n_dht  <- sum(df$treatment == "DHT[10nM]")
p_wilcox <- wilcox.test(ARscore1 ~ treatment, data = df)$p.value
mean_etoh <- mean(df$ARscore1[df$treatment == "EtOH"])
mean_dht  <- mean(df$ARscore1[df$treatment == "DHT[10nM]"])
fc_mean   <- if (mean_etoh != 0) mean_dht / mean_etoh else NA_real_

label_txt <- sprintf(
  "AR+ only (%s in %s)\nWilcoxon p = %.3g\nMean FC (DHT/EtOH) = %.3f\nn: EtOH=%d, DHT=%d",
  ar_symbol, assay_use, p_wilcox, fc_mean, n_etoh, n_dht
)

# --- 3) plot for AR+ cells ---
p <- ggplot(df, aes(x = treatment, y = ARscore1)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.25, size = 0.6) +
  labs(x = NULL, y = "AR activation score (ARscore1)",
       title = "AR activation score: DHT[10nM] vs EtOH (AR+ cells only)") +
  annotate("text", x = 1.5, y = max(df$ARscore1, na.rm = TRUE),
           label = label_txt, hjust = 0.5, vjust = -0.4, size = 3.4) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
print(p)
label_txt
