# Inflammatory proteins

library(tidyverse)
library(limma)

# Basic models
npx_df <- olink_dep %>%
  select(id, il8_24:csf1_24)

pheno_df <- olink_dep %>%
  select(id, ageyr_24, c_ageyr_24, sex, mat_ed, mat_soc, bmi_10, z_bmi_10,
         bmi_24, z_bmi_24, smoking_24, auditc_24, z_auditc_24, class,
         cisr_dep_24)

rownames(npx_df) <- npx_df$id
npx_df <- npx_df %>% select(-id)

design <- model.matrix(~0 + class + sex + mat_ed + mat_soc, pheno_df)
fit <- lmFit(t(npx_df), design)
contr <- makeContrasts(class2 - class1, class3 - class1,
                       class4 - class1, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2, proportion = 0.1)

top.table21 <- topTable(fit2, coef = 1, sort.by = "P", n = Inf,
                        adjust = "BH", confint = TRUE) %>%
  rownames_to_column(var = "protein") %>%
  mutate(FC = 2^(abs(logFC)))

top.table31 <- topTable(fit2, coef = 2, sort.by = "P", n = Inf,
                        adjust = "BH", confint = TRUE) %>%
  rownames_to_column(var = "protein") %>%
  mutate(FC = 2^(abs(logFC)))

top.table41 <- topTable(fit2, coef = 3, sort.by = "P", n = Inf,
                        adjust = "BH", confint = TRUE) %>%
  rownames_to_column(var = "protein") %>%
  mutate(FC = 2^(abs(logFC)))

res_comb <- rbind(top.table21, top.table31, top.table41)

# Adjusted models
design_adj <- model.matrix(~0 + class + sex + mat_ed + mat_soc + z_bmi_10, pheno_df)
fit_adj <- lmFit(t(npx_df), design_adj)
contr_adj <- makeContrasts(class2 - class1, class3 - class1,
                           class4 - class1, levels = design_adj)
fit2_adj <- contrasts.fit(fit_adj, contr_adj)
fit2_adj <- eBayes(fit2_adj, proportion = 0.1)

top.table21_adj <- topTable(fit2_adj, coef = 1, sort.by = "P", n = Inf,
                            adjust = "BH", confint = TRUE) %>%
  rownames_to_column(var = "protein") %>%
  mutate(FC = 2^(abs(logFC)))

top.table31_adj <- topTable(fit2_adj, coef = 2, sort.by = "P", n = Inf,
                            adjust = "BH", confint = TRUE) %>%
  rownames_to_column(var = "protein") %>%
  mutate(FC = 2^(abs(logFC)))

top.table41_adj <- topTable(fit2_adj, coef = 3, sort.by = "P", n = Inf,
                            adjust = "BH", confint = TRUE) %>%
  rownames_to_column(var = "protein") %>%
  mutate(FC = 2^(abs(logFC)))

res_comb_adj <- rbind(top.table21_adj, top.table31_adj, top.table41_adj)