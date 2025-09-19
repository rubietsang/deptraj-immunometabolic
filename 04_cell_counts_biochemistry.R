# Blood cell counts and clinical biochemistry

library(tidyverse)
library(broom)

# Basic models
m1 <- blood_dep %>%
  gather(marker, value, c(ln_wbc_24:ln_ast_alt_24)) %>%
  group_by(marker) %>%
  do(tidy(lm(value ~ class + sex + mat_ed + mat_soc, data = .),
          conf.int = TRUE, conf.level = 0.95)) %>%
  filter(term == "class2" | term == "class3" | term == "class4") %>%
  cbind(., adj.p = p.adjust(.$p.value, "BH")) %>%
  mutate(marker = sub("_24", "", marker),
         OR = exp(estimate) - 0.01,
         OR.lci = exp(conf.low) - 0.01,
         OR.uci = exp(conf.high) - 0.01,
         model = "Basic",
         label = case_when(adj.p < 0.1 ~ "*"))

# Adjusted models
m2 <- blood_dep %>%
  gather(marker, value, c(ln_wbc_24:ln_ast_alt_24)) %>%
  group_by(marker) %>%
  do(tidy(lm(value ~ class + sex + mat_ed + mat_soc + z_bmi_10, data = .),
          conf.int = TRUE, conf.level = 0.95)) %>%
  filter(term == "class2" | term == "class3" | term == "class4") %>%
  cbind(., adj.p = p.adjust(.$p.value, "BH")) %>%
  mutate(marker = sub("_24", "", marker),
         OR = exp(estimate) - 0.01,
         OR.lci = exp(conf.low) - 0.01,
         OR.uci = exp(conf.high) - 0.01,
         model = "Adjusted",
         label = case_when(adj.p < 0.1 ~ "*"))
