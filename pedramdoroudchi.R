# SAS Task
# by Pedram Doroudchi


library(tidyverse)
library(flextable)
library(officer)

setwd("/Users/pedramdoroudchi/Desktop")


# create data
FCD <- c(
  rep('Positive', 145), rep('Negative', 3), rep('Positive', 1), rep('Negative', 1), 
  rep('Positive', 5), rep('Negative', 1), rep('Positive', 2), rep('Negative', 142)
)

CCD1 <- c(rep('Positive', 150), rep('Negative', 150))

CCD2 <- c(rep('Positive', 148), rep('Negative', 2), rep('Positive', 6), rep('Negative', 144)) 

combined_dat <- data.frame(FCD = FCD, CCD1 = CCD1, CCD2 = CCD2)


# create confusion matrix
confusion_matrix_ccd_pos <- combined_dat %>% 
  filter(CCD1 == 'Positive') %>% 
  group_by(FCD) %>% 
  reframe(
    `CCD1+ & CCD2+` = sum(CCD2 == 'Positive'),
    `CCD1+ & CCD2-` = sum(CCD2 == 'Negative'),
    `Total CCD1+` = n()
  ) %>% 
  arrange(desc(FCD))

confusion_matrix_ccd_pos_tot <- combined_dat %>% 
  filter(CCD1 == 'Positive') %>% 
  reframe(
    FCD = 'Total',
    `CCD1+ & CCD2+` = sum(CCD2 == 'Positive'),
    `CCD1+ & CCD2-` = sum(CCD2 == 'Negative'),
    `Total CCD1+` = n()
  ) %>% 
  arrange(desc(FCD))

confusion_matrix_ccd_pos_w_tot <- rbind(
  confusion_matrix_ccd_pos,
  confusion_matrix_ccd_pos_tot
)

confusion_matrix_ccd_neg <- combined_dat %>% 
  filter(CCD1 == 'Negative') %>% 
  group_by(FCD) %>% 
  reframe(
    `CCD1- & CCD2+` = sum(CCD2 == 'Positive'),
    `CCD1- & CCD2-` = sum(CCD2 == 'Negative'),
    `Total CCD1-` = n()
  ) %>% 
  arrange(desc(FCD))

confusion_matrix_ccd_neg_tot <- combined_dat %>% 
  filter(CCD1 == 'Negative') %>% 
  reframe(
    FCD = 'Total',
    `CCD1- & CCD2+` = sum(CCD2 == 'Positive'),
    `CCD1- & CCD2-` = sum(CCD2 == 'Negative'),
    `Total CCD1-` = n()
  ) %>% 
  arrange(desc(FCD))

confusion_matrix_ccd_neg_w_tot <- rbind(
  confusion_matrix_ccd_neg,
  confusion_matrix_ccd_neg_tot
)

confusion_matrix <- merge(
  confusion_matrix_ccd_pos_w_tot,
  confusion_matrix_ccd_neg_w_tot
) %>% 
  mutate(
    FCD = factor(FCD, levels = c('Positive', 'Negative', 'Total'), labels = c('FCD+', 'FCD-', 'Total'))
  ) %>% 
  arrange(FCD) %>% 
  mutate(
    FCD = as.character(FCD)
  )


# zeta estimates and bootstrap confidence intervals
# is PPA between FCD and CCD non-inferior to that of two CCD replicates by a margin of .05?

# ZPPA1 (PPA_C1C2 - PPA_C1F)
zPPA1_hat <- (confusion_matrix[confusion_matrix$FCD == 'Total', "CCD1+ & CCD2+"] - confusion_matrix[confusion_matrix$FCD == 'FCD+', "Total CCD1+"]) / 
  confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1+"]

# 95% bootstrap confidence interval
set.seed(1)
zPPA1_hats <- c()
for (i in 1:1000) {
  resampled_dat <- combined_dat %>% 
    filter(CCD1 == 'Positive') %>% 
    slice_sample(n = 150, replace = T)
  zPPA1_hats[i] <- (sum(resampled_dat$CCD2 == 'Positive') - sum(resampled_dat$FCD == 'Positive')) / nrow(resampled_dat)
}
zPPA1_lower_CI = quantile(zPPA1_hats, .025)
zPPA1_upper_CI = quantile(zPPA1_hats, .975)


# ZPPA2 (PPA_C2C1 - PPA_C2F)
prev <- 0.3
zPPA2_hat <- (
  (prev*confusion_matrix[confusion_matrix$FCD == 'Total', "CCD1+ & CCD2+"]*confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1-"]) - 
    (prev*confusion_matrix[confusion_matrix$FCD == 'FCD+', "CCD1+ & CCD2+"]*confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1-"]) - 
    ((1-prev)*confusion_matrix[confusion_matrix$FCD == 'FCD+', "CCD1- & CCD2+"]*confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1+"])
) / 
  (
    (prev*confusion_matrix[confusion_matrix$FCD == 'Total', "CCD1+ & CCD2+"]*confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1-"]) + 
      ((1-prev)*confusion_matrix[confusion_matrix$FCD == 'Total', "CCD1- & CCD2+"]*confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1+"])
  )

# 95% bootstrap confidence interval
zPPA2_hats <- c()
for (i in 1:1000) {
  resampled_dat <- combined_dat %>% 
    slice_sample(n = 300, replace = T)
  zPPA2_hats[i] <- (
    (prev*sum(resampled_dat$CCD1 == 'Positive' & resampled_dat$CCD2 == 'Positive')*sum(resampled_dat$CCD1 == 'Negative')) - 
      (prev*sum(resampled_dat$CCD2 == 'Positive' & resampled_dat$CCD1 == 'Positive' & resampled_dat$FCD == 'Positive')*sum(resampled_dat$CCD1 == 'Negative')) - 
      ((1-prev)*sum(resampled_dat$CCD2 == 'Positive' & resampled_dat$CCD1 == 'Negative' & resampled_dat$FCD == 'Positive')*sum(resampled_dat$CCD1 == 'Positive'))
    ) / 
    (
      (prev*sum(resampled_dat$CCD1 == 'Positive' & resampled_dat$CCD2 == 'Positive')*sum(resampled_dat$CCD1 == 'Negative')) + 
        ((1-prev)*sum(resampled_dat$CCD1 == 'Negative' & resampled_dat$CCD2 == 'Positive')*sum(resampled_dat$CCD1 == 'Positive'))
  )
}
zPPA2_lower_CI = quantile(zPPA2_hats, .025)
zPPA2_upper_CI = quantile(zPPA2_hats, .975)


# ZNPA1 (NPA_C1C2 - NPA_C1F)
zNPA1_hat <- (confusion_matrix[confusion_matrix$FCD == 'Total', "CCD1- & CCD2-"] - confusion_matrix[confusion_matrix$FCD == 'FCD-', "Total CCD1-"]) / 
  confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1-"]

# 95% bootstrap confidence interval
zNPA1_hats <- c()
for (i in 1:1000) {
  resampled_dat <- combined_dat %>% 
    filter(CCD1 == 'Negative') %>% 
    slice_sample(n = 150, replace = T)
  zNPA1_hats[i] <- (sum(resampled_dat$CCD2 == 'Negative') - sum(resampled_dat$FCD == 'Negative')) / nrow(resampled_dat)
}
zNPA1_lower_CI = quantile(zNPA1_hats, .025)
zNPA1_upper_CI = quantile(zNPA1_hats, .975)


# ZNPA2 (NPA_C2C1 - NPA_C2F)
zNPA2_hat <- (
  ((1-prev)*confusion_matrix[confusion_matrix$FCD == 'Total', "CCD1- & CCD2-"]*confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1+"]) - 
    ((1-prev)*confusion_matrix[confusion_matrix$FCD == 'FCD-', "CCD1- & CCD2-"]*confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1+"]) - 
    (prev*confusion_matrix[confusion_matrix$FCD == 'FCD-', "CCD1+ & CCD2-"]*confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1-"])
) / 
  (
    ((1-prev)*confusion_matrix[confusion_matrix$FCD == 'Total', "CCD1- & CCD2-"]*confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1+"]) + 
      (prev*confusion_matrix[confusion_matrix$FCD == 'Total', "CCD1+ & CCD2-"]*confusion_matrix[confusion_matrix$FCD == 'Total', "Total CCD1-"])
  )

# 95% bootstrap confidence interval
zNPA2_hats <- c()
for (i in 1:1000) {
  resampled_dat <- combined_dat %>% 
    slice_sample(n = 300, replace = T)
  zNPA2_hats[i] <- (
    ((1-prev)*sum(resampled_dat$CCD1 == 'Negative' & resampled_dat$CCD2 == 'Negative')*sum(resampled_dat$CCD1 == 'Positive')) - 
      ((1-prev)*sum(resampled_dat$CCD2 == 'Negative' & resampled_dat$CCD1 == 'Negative' & resampled_dat$FCD == 'Negative')*sum(resampled_dat$CCD1 == 'Positive')) - 
      (prev*sum(resampled_dat$CCD2 == 'Negative' & resampled_dat$CCD1 == 'Positive' & resampled_dat$FCD == 'Negative')*sum(resampled_dat$CCD1 == 'Negative'))
  ) / 
    (
      ((1-prev)*sum(resampled_dat$CCD1 == 'Negative' & resampled_dat$CCD2 == 'Negative')*sum(resampled_dat$CCD1 == 'Positive')) + 
        (prev*sum(resampled_dat$CCD1 == 'Positive' & resampled_dat$CCD2 == 'Negative')*sum(resampled_dat$CCD1 == 'Negative'))
    )
}
zNPA2_lower_CI = quantile(zNPA2_hats, .025)
zNPA2_upper_CI = quantile(zNPA2_hats, .975)


# create analysis table
analysis_tbl <- data.frame(
  Prevalence = '30.0%',
  'ZPPA1 Estimate' = sprintf('%.3f', zPPA1_hat),
  'ZPPA1 95% CI' = paste0('( ', sprintf('%.3f', zPPA1_lower_CI), ', ', sprintf('%.3f', zPPA1_upper_CI), ')'),
  'ZPPA2 Estimate' = sprintf('%.3f', zPPA2_hat),
  'ZPPA2 95% CI' = paste0('( ', sprintf('%.3f', zPPA2_lower_CI), ', ', sprintf('%.3f', zPPA2_upper_CI), ')'),
  'ZNPA1 Estimate' = sprintf('%.3f', zNPA1_hat),
  'ZNPA1 95% CI' = paste0('( ', sprintf('%.3f', zNPA1_lower_CI), ', ', sprintf('%.3f', zNPA1_upper_CI), ')'),
  'ZNPA2 Estimate' = sprintf('%.3f', zNPA2_hat),
  'ZNPA2 95% CI' = paste0('( ', sprintf('%.3f', zNPA2_lower_CI), ', ', sprintf('%.3f', zNPA2_upper_CI), ')')
)


# flextable

# confusion matrix
confusion_matrix_ft <- flextable(confusion_matrix)

# add header
confusion_matrix_ft <- add_header_row(confusion_matrix_ft, values = c("", "CCD1+", "CCD1-"), colwidths = c(1, 3, 3))

# set alignment for headers, footer and body
confusion_matrix_ft <- align(confusion_matrix_ft, part = "header", align = "center")
confusion_matrix_ft <- align(confusion_matrix_ft, part = "body", align = "left")
confusion_matrix_ft <- valign(confusion_matrix_ft, part = "body", valign = "top")

# configure table fonts and styles
confusion_matrix_ft <- font(confusion_matrix_ft, fontname = "Times New Roman", part = "all")
confusion_matrix_ft <- fontsize(confusion_matrix_ft, size = 10, part = "all")

# set header background color
#confusion_matrix_ft <- bg(confusion_matrix_ft, i = ~ `%Agr` != "100.0", bg = "#FDE9D9", part = "body")

# merge cells with same values in certain columns
#confusion_matrix_ft <- merge_v(confusion_matrix_ft, j = 1:3)
confusion_matrix_ft <- merge_at(confusion_matrix_ft, part = "header", i=1:2, j=1)
confusion_matrix_ft <- compose(confusion_matrix_ft, i = 1:2, part = "header", j = 1, as_paragraph(as_chunk(c(''))))
confusion_matrix_ft <- compose(confusion_matrix_ft, i = 2, part = "header", j = 2, as_paragraph(as_chunk(c("CCD2+"))))
confusion_matrix_ft <- compose(confusion_matrix_ft, i = 2, part = "header", j = 3, as_paragraph(as_chunk(c('CCD2-'))))
confusion_matrix_ft <- compose(confusion_matrix_ft, i = 2, part = "header", j = 4, as_paragraph(as_chunk(c("Total"))))
confusion_matrix_ft <- compose(confusion_matrix_ft, i = 2, part = "header", j = 5, as_paragraph(as_chunk(c('CCD2+'))))
confusion_matrix_ft <- compose(confusion_matrix_ft, i = 2, part = "header", j = 6, as_paragraph(as_chunk(c("CCD2-"))))
confusion_matrix_ft <- compose(confusion_matrix_ft, i = 2, part = "header", j = 7, as_paragraph(as_chunk(c("Total"))))

# add cell and table borders
confusion_matrix_ft <- border_inner(confusion_matrix_ft)
confusion_matrix_ft <- border_outer(confusion_matrix_ft)

# set padding, autofit and print
confusion_matrix_ft <- padding(confusion_matrix_ft, padding.top = 72*.02, padding.bottom = 72*.02, padding.left = 72*.04, padding.right = 72*.04, part = 'all')
confusion_matrix_ft <- autofit(confusion_matrix_ft)
confusion_matrix_ft <- set_table_properties(confusion_matrix_ft, layout = "autofit")
print(confusion_matrix_ft)


# analysis table
analysis_ft <- flextable(analysis_tbl)

# add header
analysis_ft <- add_header_row(analysis_ft, values = c("", "ZPPA1", "ZPPA2", "ZNPA1", "ZNPA2"), colwidths = c(1, 2, 2, 2, 2))

# set alignment for headers, footer and body
analysis_ft <- align(analysis_ft, part = "all", align = "center")
analysis_ft <- valign(analysis_ft, part = "all", valign = "center")

# configure table fonts and styles
analysis_ft <- font(analysis_ft, fontname = "Times New Roman", part = "all")
analysis_ft <- fontsize(analysis_ft, size = 10, part = "all")

# set header background color
analysis_ft <- bg(analysis_ft, bg = "#D9D9D9", part = "header")

# merge cells with same values in certain columns
#analysis_ft <- merge_v(analysis_ft, j = 1:3)
analysis_ft <- compose(analysis_ft, i = 2, part = "header", j = 1, as_paragraph(as_chunk(c("Prevalence"))))
analysis_ft <- compose(analysis_ft, i = 2, part = "header", j = 2, as_paragraph(as_chunk(c('Estimate'))))
analysis_ft <- compose(analysis_ft, i = 2, part = "header", j = 3, as_paragraph(as_chunk(c("95% CI"))))
analysis_ft <- compose(analysis_ft, i = 2, part = "header", j = 4, as_paragraph(as_chunk(c('Estimate'))))
analysis_ft <- compose(analysis_ft, i = 2, part = "header", j = 5, as_paragraph(as_chunk(c("95% CI"))))
analysis_ft <- compose(analysis_ft, i = 2, part = "header", j = 6, as_paragraph(as_chunk(c('Estimate'))))
analysis_ft <- compose(analysis_ft, i = 2, part = "header", j = 7, as_paragraph(as_chunk(c("95% CI"))))
analysis_ft <- compose(analysis_ft, i = 2, part = "header", j = 8, as_paragraph(as_chunk(c('Estimate'))))
analysis_ft <- compose(analysis_ft, i = 2, part = "header", j = 9, as_paragraph(as_chunk(c("95% CI"))))

# add cell and table borders
analysis_ft <- border_inner(analysis_ft)
analysis_ft <- border_outer(analysis_ft)

# set padding, autofit and print
analysis_ft <- padding(analysis_ft, padding.top = 72*.02, padding.bottom = 72*.02, padding.left = 72*.04, padding.right = 72*.04, part = 'all')
analysis_ft <- autofit(analysis_ft)
analysis_ft <- set_table_properties(analysis_ft, layout = "autofit")
print(analysis_ft)


# export
table_doc <- read_docx() %>% 
  body_add_flextable(confusion_matrix_ft) %>% 
  body_add_break() %>% 
  body_add_flextable(analysis_ft)

print(table_doc, target = paste0("analysis_tables_R_PD", Sys.Date(), ".docx"))

