Genotype-phenotype distributions
================
Jaeda Patton
2023-05-10

## Loading data and functions

``` r
# load general functions
source(file.path("..", "scripts", "general_functions.R"))

# reading in cleaned data from binned sort experiment

# reading in complete fluorescence data
AncSR1_meanF_data <- read.csv(file.path("..", "results", 
                                        "mutation_effects_model",
                                        "AncSR1_complete_data.csv.gz"), 
                              stringsAsFactors = TRUE)
AncSR2_meanF_data <- read.csv(file.path("..", "results", 
                                        "mutation_effects_model",
                                        "AncSR2_complete_data.csv.gz"), 
                              stringsAsFactors = TRUE)
meanF_data <- bind_rows(AncSR1 = AncSR1_meanF_data, AncSR2 = AncSR2_meanF_data, 
                        .id = "bg")
rm(AncSR1_meanF_data, AncSR2_meanF_data)
```

## Analyses

### How many functional genotypes are on each background?

We define “functional” variants as those with fluorescence at least as
high as that of AncSR2-WT:SRE1. This is our reference variant.

``` r
# defining minimum fluorescence cutoff for active variants
AncSR2WT_SRE1_data <- meanF_data %>% 
  filter(AA_var == "GSKV", bg == "AncSR2", RE == "SRE1 (AA)")
```

To classify variants as functional vs. nonfunctional, we will test
whether their fluorescence is significantly lower than that of the
reference variant. If it is, then variants are classified as
nonfunctional, otherwise they are considered functional.

Since our fluorescence estimates come from different sources (observed
vs. inferred), we will use different tests to classify the variants
whose fluorescence is estimated from the binned sort data, the debulk
sort data, and the mutation effects models:

1.  Variants that are classified as having null fluorescence from the
    debulk sort data (see `data_cleaning.md`) are considered
    nonfunctional.
2.  Variants whose fluorescence is estimated from the binned sort data
    have been observed in at least two replicates (most in three
    replicates). We will use a t-test to test the null hypothesis that
    their fluorescence is greater than or equal to that of the reference
    variant.
3.  Variants whose fluorescence is estimated from the mutation effects
    model have only a single predicted fluorescence value. To generate a
    p-value for these variants, we will use a nonparametric bootstrap
    test to test the same null hypothesis as in the t-test. This
    approach takes advantage of the fact that we can use the prediction
    error observed during model cross-validation as an estimate of the
    error associated with each predicted variant’s fluorescence. For
    each variant with predicted fluorescence $F_{\text{test}}$, we will
    take the distribution of cross-validation residuals in the interval
    $F_{\text{test}} \pm 0.1$ as an estimate of the error distribution
    for that variant. We will then create bootstrap fluorescence samples
    from this error distribution and compute the p-value as the fraction
    of the bootstrap sample with fluorescence greater than or equal to
    the reference variant. This will be done separately for AncSR1
    vs. AncSR2 variants, since the models were fit separately for each
    protein background and thus have different error distributions.

A false-discovery rate threshold of 0.1 will be used to classify
variants as functional vs. nonfunctional.

``` r
# Compute p-value for calling variants functional that were observed in the 
# binned sort experiment.
#
# H0: Variant is at least as fluorescent as reference variant.
# HA: Variant is less fluorescent than reference variant.
# Use a t-test and FDR-adjusted p-value. Variants with padj < 0.1 are called nonfunctional.

# get active (not null) variants
# meanF_data_active <- meanF_data %>% 
#   filter(active) %>%
#   mutate(RE = factor(RE, levels = levels(REs[[1]])))

# test for normality of measurement error for observed active variants
# m.error.active <- meanF_data_active %>%
#   filter(type == "binned") %>%
#   select(avg_meanF, meanF_REP1:meanF_REP4) %>%
#   pivot_longer(meanF_REP1:meanF_REP4, names_to = "rep", values_to = "obs") %>%
#   mutate(e = obs - avg_meanF) %>%
#   pull(e)
# m.error.active <- m.error.active[!is.na(m.error.active)]
# hist(m.error.active, freq = FALSE)
# abline(v = 0)
# error.normal.fit <- fitdistr(m.error.active, "normal")
# lines(seq(-.2, .2, 0.01), dnorm(seq(-.2, .2, 0.01), error.normal.fit$estimate[1], error.normal.fit$estimate[2]))

# compute p-values with Mann-Whitney U test
# pfunctional <- apply(meanF_data_active %>% 
#                        filter(type == "binned") %>%
#                        select(meanF_REP1:meanF_REP4), 
#                      1, function(x) 
#                        wilcox.test(AncSR2WT_SRE1_data %>% 
#                                      select(meanF_REP1:meanF_REP3) %>%
#                                      as.numeric(), x,
#                                    alternative = "greater")$p.value)
# # FDR correction
# padjfunctional <- p.adjust(pfunctional, "fdr")

# pfunctionalt <- meanF_data_active %>%
#   filter(type == "binned") %>%
#   select(meanF_REP1:meanF_REP4) %>%
#   apply(1, function(x) t.test(as.numeric(x),
#                               AncSR2WT_SRE1_data %>%
#                                 select(meanF_REP1:meanF_REP3) %>%
#                                 as.numeric(),
#                               "less", na.action = "na.omit")$p.value)
# padjfunctionalt <- p.adjust(pfunctionalt, "fdr")
# data.frame(f = meanF_data_active %>% filter(type == "binned") %>% pull(avg_meanF),
#            sig = padjfunctionalt <= 0.1) %>%
#   ggplot(aes(x = f, fill = sig)) +
#   geom_histogram() +
#   geom_vline(xintercept = AncSR2WT_SRE1_data$avg_meanF)

if(!file.exists(file.path(results_dir, "pbinned.rda"))) {
  pbinned <- apply(select(data, meanF_REP1:meanF_REP4), 1, function(x)
      t.test(as.numeric(x), 
             as.numeric(select(AncSR2WT_SRE1_data, meanF_REP1:meanF_REP4)),
             "less", na.action = "na.omit")$p.value)
  save(pbinned, file = file.path(results_dir, "pbinned.rda"))
} else load(file.path(results_dir, "pbinned.rda"))

padjbinned <- p.adjust(pbinned, "fdr")
```

``` r
# For variants with predicted fluorescence from reference-free model, test for
# fluorescence less than the reference with a nonparametric bootstrap test,
# where bootstrapped fluorescence samples are taken from the cross-validation
# error distribution for variants with fluorescence similar to that of the test
# variant. This method helps to account for the non-normal error distributions
# observed for high predicted fluorescence values.
#
# H0: Variant is at least as fluorescent as reference variant.
# HA: Variant is less fluorescent than reference variant.
# Adjust p-value using FDR correction. Variants with padj < 0.1 are called 
# nonfunctional.

# first load the cross-validation fits from the model fitting
load(file.path("..", "results", "mutation_effects_model", 
               "AncSR1.cv.pred.fine.rda"))
load(file.path("..", "results", "mutation_effects_model", "AncSR1_foldid.rda"))
load(file.path("..", "results", "mutation_effects_model", 
               "AncSR1_model_data.rda"))
AncSR1.cv.pred <- lapply(AncSR1.cv.pred.fine, function(x) as.numeric(x[,10]))

AncSR1.cv.obs <- lapply(1:10, function(x) 
  AncSR1_model_data$avg_meanF[AncSR1_foldid == x])

rm(AncSR1.cv.pred.fine, AncSR1_foldid, AncSR1_model_data)

# concatenate fits across all 10 cross-validation folds
AncSR1.cv <- data.frame(pred = unlist(AncSR1.cv.pred, use.names = FALSE),
                        obs = unlist(AncSR1.cv.obs, use.names = FALSE))
# ggplot(AncSR1.cv, aes(x = pred, y = obs)) +
#   geom_bin_2d() +
#   scale_fill_continuous(trans = log10plus1) +
#   geom_hline(yintercept = min(meanF_data %>% filter(type == "predicted", bg == "AncSR1", active) %>% pull(avg_meanF)), linetype = 2) +
#   geom_vline(xintercept = min(meanF_data %>% filter(type == "predicted", bg == "AncSR1", active) %>% pull(avg_meanF)), linetype = 2) +
#   geom_hline(yintercept = AncSR2WT_SRE1_data$avg_meanF, linetype = 2, color = "gray") +
#   geom_vline(xintercept = AncSR2WT_SRE1_data$avg_meanF, linetype = 2, color = "gray")

# compute residuals
AncSR1.cv$res <- AncSR1.cv$pred - AncSR1.cv$obs

# Bootstrap predicted fluorescence based on CV error distributions. For each
# predicted variant, sample from the distribution of CV residuals concatenated
# across all 10 CV folds, within a range of +/- 0.1 of the predicted
# fluorescence. Sample 250 bootstrap samples per replicate.

if(!file.exists(file.path(results_dir, "bspred.AncSR1.rda"))) {
  # parallel processing
  cores <- 24
  cl <- parallel::makeCluster(cores, "FORK", outfile = "")
  registerDoParallel(cl)
  bspred.AncSR1 <- foreach(i = 1:cores, .combine = 'cbind') %dopar% {
    data <- meanF_data %>% filter(type == "predicted", bg == "AncSR1")
    
    # split data into chunks
    size <- nrow(data)
    chunksize <- ceiling(size / cores)
    chunk <- data[((i - 1) * chunksize + 1):min(i * chunksize, size),]
    
    # create bootstrap reps
    sapply(chunk$avg_meanF, function(x) {
      # residual distribution centered around meanF of test variant
      res <- AncSR1.cv %>% filter(pred > x - 0.1 & pred < x + 0.1) %>% pull(res)
      # bootstrap
      sample(x + res, 250, replace = TRUE)
    })
  }
  stopCluster(cl)
  save(bspred.AncSR1, file = file.path(results_dir, "bspred.AncSR1.rda"))
} else load(file.path(results_dir, "bspred.AncSR1.rda"))

# AncSR1.pred.bs <- sapply(meanF_data %>% 
#                            filter(type == "predicted", bg == "AncSR1") %>%
#                            pull(avg_meanF), function(x) {
#   res <- AncSR1.cv %>% filter(pred > x - 0.1 & pred < x + 0.1) %>% pull(res)
#   sample(x + res, 250, replace = TRUE)
# })

# Compute p-value for functional variants as fraction of bootstrap replicates
# that are greater than or equal to that of the mean AncSR2:SRE1 WT variant.
if(!file.exists(file.path(results_dir, "ppredicted.AncSR1.rda"))) {
  ppredicted.AncSR1 <- apply(bspred.AncSR1, 2, function(x) 
    sum(x >= AncSR2WT_SRE1_data$avg_meanF) / 250)
  save(ppredicted.AncSR1, file = file.path(results_dir, "ppredicted.AncSR1.rda"))
} else load(file.path(results_dir, "ppredicted.AncSR1.rda"))
# FDR correction
padjpredicted.AncSR1 <- p.adjust(ppredicted.AncSR1, "fdr")
```

``` r
# For variants with predicted fluorescence from reference-free model, test for
# fluorescence less than the reference with a nonparametric bootstrap test,
# where bootstrapped fluorescence samples are taken from the cross-validation
# error distribution for variants with fluorescence similar to that of the test
# variant. This method helps to account for the non-normal error distributions
# observed for high predicted fluorescence values.
#
# H0: Variant is at least as fluorescent as reference variant.
# HA: Variant is less fluorescent than reference variant.
# Adjust p-value using FDR correction. Variants with padj < 0.1 are called 
# nonfunctional.

# first load the cross-validation fits from the model fitting
load(file.path("..", "results", "mutation_effects_model", 
               "AncSR2.cv.pred.fine.rda"))
load(file.path("..", "results", "mutation_effects_model", "AncSR2_foldid.rda"))
load(file.path("..", "results", "mutation_effects_model", 
               "AncSR2_model_data.rda"))
AncSR2.cv.pred <- lapply(AncSR2.cv.pred.fine, function(x) as.numeric(x[,13]))

AncSR2.cv.obs <- lapply(1:10, function(x)
  AncSR2_model_data$avg_meanF[AncSR2_foldid == x])

rm(AncSR2.cv.pred.fine, AncSR2_foldid, AncSR2_model_data)

# concatenate fits across all 10 cross-validation folds
AncSR2.cv <- data.frame(pred = unlist(AncSR2.cv.pred, use.names = FALSE),
                        obs = unlist(AncSR2.cv.obs, use.names = FALSE))

# compute residuals
AncSR2.cv$res <- AncSR2.cv$pred - AncSR2.cv$obs
# Bootstrap predicted fluorescence based on CV error distributions. For each
# predicted variant, sample from the distribution of CV residuals concatenated
# across all 10 CV folds, within a range of +/- 0.1 of the predicted
# fluorescence. Sample 250 bootstrap samples per replicate.

if(!file.exists(file.path(results_dir, "bspred.AncSR2.rda"))) {
  # parallel processing
  cores <- 24
  cl <- parallel::makeCluster(cores, "FORK", outfile = "")
  registerDoParallel(cl)
  bspred.AncSR2 <- foreach(i = 1:cores, .combine = 'cbind') %dopar% {
    data <- meanF_data %>% filter(type == "predicted", bg == "AncSR2")
    
    # split data into chunks
    size <- nrow(data)
    chunksize <- ceiling(size / cores)
    chunk <- data[((i - 1) * chunksize + 1):min(i * chunksize, size),]
    
    # create bootstrap reps
    sapply(chunk$avg_meanF, function(x) {
      # residual distribution centered around meanF of test variant
      res <- AncSR2.cv %>% filter(pred > x - 0.1 & pred < x + 0.1) %>% pull(res)
      # bootstrap
      sample(x + res, 250, replace = TRUE)
    })
  }
  stopCluster(cl)
  save(bspred.AncSR2, file = file.path(results_dir, "bspred.AncSR2.rda"))
} else load(file.path(results_dir, "bspred.AncSR2.rda"))

# AncSR2.pred.bs <- sapply(meanF_data %>% 
#                            filter(type == "predicted", bg == "AncSR2") %>%
#                            pull(avg_meanF), function(x) {
#   res <- AncSR2.cv %>% filter(pred > x - 0.1 & pred < x + 0.1) %>% pull(res)
#   sample(x + res, 250, replace = TRUE)
# })

# Compute p-value for functional variants as fraction of bootstrap replicates
# that are greater than or equal to that of the mean AncSR2:SRE1 WT variant.
if(!file.exists(file.path(results_dir, "ppredicted.AncSR2.rda"))) {
  ppredicted.AncSR2 <- apply(bspred.AncSR2, 2, function(x) 
    sum(x >= AncSR2WT_SRE1_data$avg_meanF) / 250)
  save(ppredicted.AncSR2, file = file.path(results_dir, "ppredicted.AncSR2.rda"))
} else load(file.path(results_dir, "ppredicted.AncSR2.rda"))
# FDR correction
padjpredicted.AncSR2 <- p.adjust(ppredicted.AncSR2, "fdr")
```

``` r
# classify variants as functional if padj >= 0.1 (not significantly less 
# fluorescent than AncSR2:SRE1 WT)
meanF_data <- rbind(meanF_data %>% 
                      filter(type == "binned") %>%
                      mutate(functional = padjbinned >= 0.1),
                    meanF_data %>%
                      filter(type == "predicted", bg == "AncSR1") %>%
                      mutate(functional = padjpredicted.AncSR1 >= 0.1),
                    meanF_data %>%
                      filter(type == "predicted", bg == "AncSR2") %>%
                      mutate(functional = padjpredicted.AncSR2 >= 0.1),
                    meanF_data %>%
                      filter(type == "debulk") %>%
                      mutate(functional = FALSE)) %>%
  arrange(bg, AA_var, RE)

# export data
write.csv(meanF_data, file = gzfile(file.path(results_dir, "meanF_data_fxnal.csv.gz")),
          row.names = FALSE)

# plot histogram of fluorescence colored by functional vs. not functional
ggplot(meanF_data, aes(x = avg_meanF, fill = type, alpha = functional)) +
  geom_histogram(position = "stack") +
  scale_y_continuous(trans = log10plus1, name = "Count + 1") +
  scale_alpha_manual(values = c(0.6, 1)) +
  geom_vline(xintercept = AncSR2WT_SRE1_data$avg_meanF, 
             color = "gray30", linetype = 2) +
  facet_grid(cols = vars(bg)) +
  theme_classic() +
  labs(x = "Fluorescence")
```

![](genotype_phenotype_distributions_files/figure-gfm/classifyfunctional-1.png)<!-- -->

The histograms above show the distribution of fluorescence for variants
classified as functional vs. nonfunctional. The vertical dashed line
shows the fluorescence of the reference variant.

``` r
# print number of variants (protein:RE) classified as functional on each 
# background
meanF_data %>%
  group_by(bg, type) %>%
  filter(type != "debulk") %>%
  summarize(count = sum(functional)) %>%
  pivot_wider(names_from = type, values_from = count) %>%
  mutate(all = binned + predicted, fraction.total = all / 2560000) %>%
  knitr::kable(caption = "Number of functional protein:RE variants")
```

| bg     | binned | predicted |  all | fraction.total |
|:-------|-------:|----------:|-----:|---------------:|
| AncSR1 |    192 |       114 |  306 |      0.0001195 |
| AncSR2 |   3763 |      1270 | 5033 |      0.0019660 |

Number of functional protein:RE variants

``` r
# print number of protein variants classified as functional on each background
meanF_data %>%
  group_by(bg, AA_var) %>%
  summarize(functional = sum(functional) > 0) %>%
  group_by(bg) %>%
  summarize(count = sum(functional)) %>%
  mutate(fraction.total = count / 160000) %>%
  knitr::kable(caption = "Number of functional protein variants")
```

| bg     | count | fraction.total |
|:-------|------:|---------------:|
| AncSR1 |   259 |      0.0016187 |
| AncSR2 |  2391 |      0.0149438 |

Number of functional protein variants

``` r
# histogram of fluorescence colored by FDR <= 0.15
meanF_data_active %>%
  ggplot(aes(x = avg_meanF, fill = padjfunctional <= 0.15)) +
  geom_histogram() +
  geom_vline(xintercept = AncSR2WT_SRE1_data$avg_meanF) +  # fluorescence of reference
  theme_classic() +
  labs(title = "active variants", x = "fluorescence", fill = "significant")

meanF_data_functional <- meanF_data_active %>% filter(padjfunctional > 0.15)
```

How many functional genotypes are there in each ancestral background?
Which REs do they bind to?

``` r
fontsize <- 16

# number of functional protein:RE variants on each background
a <- meanF_data_functional %>%
  group_by(bg) %>%
  count() %>%
  ggplot(aes(x = bg, y = n, fill = bg)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = bg_color(), drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Ancestral\nbackground", 
       y = "number of functional\nprotein:response element variants") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        legend.position = "none") +
  guides(x = guide_axis(angle = 45))

# number of functional protein variants on each background
a2 <- meanF_data_functional %>%
  group_by(bg) %>%
  distinct(AA_var, .keep_all = TRUE) %>%
  count() %>%
  ggplot(aes(x = bg, y = n, fill = bg)) +
  geom_col(width = 0.9) +
  # geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = bg_color(), drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Ancestral\nbackground", 
       y = "number of bound protein variants") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        legend.position = "none") +
  guides(x = guide_axis(angle = 45))

# number of functional variants on each background that bind to each RE
b <- meanF_data_functional %>%
  group_by(bg, RE) %>%
  count() %>%
  mutate(RE = factor(RE, levels = levels(REs[[1]]))) %>%
  ggplot(aes(x = RE, y = n, fill = bg)) +
  geom_col(width = 0.9) +
  facet_grid(rows = vars(bg), scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = bg_color(), drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  # geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) +
  labs(x = "Response element", y = "Number of bound protein variants", 
       fill = "Ancestral\nbackground") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = c(0.9, 0.9)) +
  guides(x = guide_axis(angle = 45))

# same as b but plotting the backgrounds next to each other
b2 <- meanF_data_functional %>%
  group_by(bg, RE) %>%
  count() %>%
  ungroup() %>%
  complete(bg, RE, fill = list(n = 0)) %>%
  mutate(RE = factor(RE, levels = levels(REs[[1]]))) %>%
  ggplot(aes(x = RE, y = n, fill = bg)) +
  geom_col(width = 0.9, position = "dodge") +
  # geom_text(aes(label = n), vjust = -0.5, color = "black",
  #           size = 3.5, position = position_dodge(width = 0.9)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = bg_color(), drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Response element", y = "Number of bound protein variants", 
       fill = "Ancestral\nbackground") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = c(0.9, 0.9)) +
  guides(x = guide_axis(angle = 45))

a + plot_spacer() + b + plot_layout(widths = c(1, 0.5, 5))
a2 + plot_spacer() + b + plot_layout(widths = c(1, 0.5, 10))
a2 + plot_spacer() + b2 + plot_layout(widths = c(1, 0.5, 10))
```

## How much overlap is there between the variants that are functional on each background?

``` r
# get number of protein variants that bind each RE in only AncSR1 background,
# only AncSR2 background, or both backgrounds
variants_bound_by_bg <- meanF_data_functional %>%
  filter(type == "exp") %>%
  select(AA_var, RE, bg) %>%
  group_by(RE) %>%
  mutate(n = 1) %>%
  pivot_wider(names_from = bg, values_from = n) %>%
  mutate(both = AncSR1 + AncSR2 == 2) %>%
  summarize(AncSR1 = sum(AncSR1, na.rm = T),
            AncSR2 = sum(AncSR2, na.rm = T), 
            both = sum(both, na.rm = T)) %>%
  mutate(AncSR1 = AncSR1 - both,
         AncSR2 = AncSR2 - both,
         RE = factor(RE, levels = c("all", levels(REs[[1]])))) %>%
  rbind(c("all", colSums(select(., 2:4)))) %>%
  mutate(AncSR1 = as.numeric(AncSR1), AncSR2 = as.numeric(AncSR2),
         both = as.numeric(both))

# compute p-value for overlap by bootstrap sampling variants in each background 
# and computing amount of random overlap
variants_bound_by_bg$p <- apply(select(variants_bound_by_bg, 2:4), 1,
                                function(x) {
                                  overlaps <- numeric(length = 1000)
                                  for(i in 1:1000) {
                                    AncSR1_sample <- sample(160000, sum(x[c(1,3)]))
                                    AncSR2_sample <- sample(160000, sum(x[2:3]))
                                    overlaps[i] <- length(intersect(AncSR1_sample, AncSR2_sample))
                                  }
                                  p <- sum(overlaps >= x[3]) / 1000
                                  p
                                })
variants_bound_by_bg$padj <- p.adjust(variants_bound_by_bg$p, 
                                      method = "bonferroni")

# plot fraction of protein variants per RE that are functional in the AncSR1,
# AncSR2, or both backgrounds
a <- variants_bound_by_bg %>%
  mutate(AncSR1_total = AncSR1 + both) %>%
  pivot_longer(2:4, names_to = "background") %>%
  filter(background != "AncSR2") %>%
  ggplot(aes(x = RE, y = value, 
             fill = factor(background, levels = c("both", "AncSR1")))) +
  geom_col(position = "stack", width = 0.9) +
  scale_fill_manual(values = c(both = "deepskyblue", bg_color("AncSR1")),
                    labels = c("AncSR1 and AncSR2", "AncSR1 only"),
                    drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  geom_text(aes(y = AncSR1_total, label = ifelse(padj <= .05, "*", "")), 
            vjust = 0, color = "black", size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Variants functional on AncSR1", x = "Response element", 
       y = "Number of protein variants", fill = "Ancestral background") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size =12),
        legend.title = element_text(size = 14),
        legend.position = c(0.8, 0.9),
        plot.caption = element_text(hjust = 0)) +
  guides(x = guide_axis(angle = 45))

# plot the distribution of fluorescence for protein variants that bind in 
# both backgrounds
fluorescence_both_backgrounds <- meanF_data_functional %>%
  select(AA_var, RE, bg, avg_meanF) %>%
  group_by(RE) %>%
  pivot_wider(names_from = bg, values_from = avg_meanF) %>%
  drop_na()
fluorescence_bg_p <- fluorescence_both_backgrounds %>%
  ungroup() %>%
  summarize(p = wilcox.test(AncSR1, AncSR2, alternative = "less")$p.value)

b <- fluorescence_both_backgrounds %>%
  pivot_longer(3:4, names_to = "background") %>%
  ggplot(aes(x = factor(RE, levels = levels(REs[[1]])), y = value, 
             fill = background)) +
  geom_boxplot() +
  scale_fill_manual(values = c(bg_color()), drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Response element", y = "Fluorescence", 
       fill = "Ancestral\nbackground") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        legend.text = element_text(size = 14),
        axis.line.x = element_line(color = "black",),
        legend.position = "right",
        plot.caption = element_text(hjust = 0)) +
  guides(x = guide_axis(angle = 45))

c <- fluorescence_both_backgrounds %>%
  pivot_longer(3:4, names_to = "background") %>%
  ggplot(aes(x = background, y = value, fill = background)) +
  geom_boxplot() +
  scale_fill_manual(values = c(bg_color()), drop = FALSE) +
  labs(title = "Variants functional on\nboth backgrounds", 
       x = "Ancestral\nbackground", y = "Fluorescence") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line.x = element_line(color = "black",),
        legend.position = "none",
        plot.caption = element_text(hjust = 0)) +
  guides(x = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(c("AncSR1", "AncSR2")), 
              test = "wilcox.test",
              test.args = list(alternative = "less"),
              map_signif_level = TRUE)

a + plot_spacer() + c + plot_layout(widths = c(5, 0.5, 2))

# plot fluorescence in each background for variants that bind in both backgrounds
d <- fluorescence_both_backgrounds %>%
  ggplot(aes(x = AncSR1, y = AncSR2)) +
  geom_point()
d
```

## Do protein variants bind specifically or promiscuously to different RE variants?

``` r
# Calculate p-values for protein variants being worse than AncSR2-WT:SRE1 on 
# multiple RE variants. The probability that it is not significantly worse
# than the AncSR2-WT:SRE1 genotype on multiple REs is simply the product of
# those probabilities across all REs that it is measured on. We calculate this 
# for all "active" variants (those who are significantly better than null 
# variants) for all combinations of RE variants that it is measured on, and use
# a Benjamini-Hochberg FDR p-value correction for each [1,16] number of RE
# variants used in the p-value computation. We take variants with padj > 0.15
# as those who are not significantly worse than AncSR2-WT:SRE1 on multiple REs.

if(!file.exists(file.path(results_dir, "multiple_REs_bound_data.rda"))) {
  multiple_REs_bound_data <- list(AncSR1 = list(), AncSR2 = list())
  multiple_REs_bound_data[["AncSR1"]][[1]] <- meanF_data_functional %>% 
    filter(bg == "AncSR1") %>% 
    mutate(p = pfunctional, padj = padjfunctional, RE1 = RE) %>%
    select(AA_var, RE1, p, padj)
  multiple_REs_bound_data[["AncSR2"]][[1]] <- meanF_data_functional %>% 
    filter(bg == "AncSR2") %>% 
    mutate(p = pfunctional, padj = padjfunctional, RE1 = RE) %>%
    select(AA_var, RE1, p, padj)
  
  # initialize data frames
  for(background in c("AncSR1", "AncSR2")) {
    for(i in 2:16) {
      df <- data.frame(AA_var = NA)
      df[sapply(1:i, function(x) paste0("RE", x))] <- NA
      df$p <- NA
      multiple_REs_bound_data[[background]][[i]] <- df
    }
  }
  
  for(background in c("AncSR1", "AncSR2")) {
    AA_vars <- meanF_data_active %>%
        filter(bg == background) %>% pull("AA_var") %>% 
        unique() %>% as.character() %>% sort()
    
    for(aa in AA_vars) {
      print(aa)
      aa_df <- meanF_data_active %>% 
          filter(AA_var == aa, bg == background) %>%
          select(RE, meanF_REP1:meanF_REP4, pfunctional)
      i <- 2
      while(i <= nrow(aa_df)) {
        # compute all possible sets of i REs for variant aa
        RE_combs <- combn(as.character(aa_df$RE), i)
        
        # compute p-value of aa being worse than AncSR2-WT:SRE1 on any of the
        # set of REs
        p <- apply(RE_combs, 2, function(x) {
          df2 <- aa_df %>% filter(RE %in% x)
          prod(df2$pfunctional)
        })
        
        new_df <- data.frame(AA_var = rep(aa, ncol(RE_combs)))
        new_df[sapply(1:i, function(x) paste0("RE", x))] <- t(RE_combs)
        new_df$p <- p
        multiple_REs_bound_data[[background]][[i]] <- 
          multiple_REs_bound_data[[background]][[i]] %>%
          rbind(new_df)
        
        i <- i + 1
      }
    }
  }
  
  # merge data for AncSR1 & AncSR2
  multiple_REs_bound_data <- lapply(1:16, function(x) 
    rbind(multiple_REs_bound_data[[1]][[x]] %>% 
            slice(-1) %>% mutate(bg = "AncSR1"), 
          multiple_REs_bound_data[[2]][[x]] %>% 
            slice(-1) %>% mutate(bg = "AncSR2")))
  
  # compute adjusted p-values
  multiple_REs_bound_data <- lapply(multiple_REs_bound_data, function(x)
    mutate(x, padj = p.adjust(p, "fdr")))
  
  save(multiple_REs_bound_data, 
       file = file.path(results_dir, "multiple_REs_bound_data.rda"))
} else load(file = file.path(results_dir, "multiple_REs_bound_data.rda"))

# TODO: check what's up with AncSR1 AAKM


# identify protein variants that pass the significance threshold for binding to 
# k or more RE variants
multiple_REs_variants <- lapply(multiple_REs_bound_data, function(x)
  x %>% filter(padj > 0.15) %>% distinct(AA_var, bg))
multiple_REs_variants <- 
  multiple_REs_variants[sapply(multiple_REs_variants, nrow) > 0]
for(i in 1:(length(multiple_REs_variants)-1)) {
  multiple_REs_variants[[i]] <- anti_join(multiple_REs_variants[[i]], multiple_REs_variants[[i+1]])
}
# convert to data frame
multiple_REs_variants <- lapply(1:length(multiple_REs_variants), function(x)
  mutate(multiple_REs_variants[[x]], nREs = x))
multiple_REs_variants <- do.call(rbind, multiple_REs_variants)


# plot the number of REs each RH variant is active on
a1 <- multiple_REs_variants %>%
  mutate(specific = factor(ifelse(nREs == 1, "specific", "promiscuous"),
                           levels = c("specific", "promiscuous"))) %>%
  group_by(bg, specific) %>%
  count(name = "nvars") %>%
  ggplot(aes(x = as.factor(specific), y = nvars, fill = bg)) +
  geom_col(width=0.9) +
  geom_text(aes(label=nvars), vjust=-0.5, color="black", size=4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_grid(rows = vars(bg), scales = "free") +
  scale_fill_manual(values=bg_color()) +
  labs(title = "", x = "", 
       y = "Number of protein variants",
       fill = "Ancestral\nbackground") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_blank()) +
  guides(x = guide_axis(angle = 45))

a2 <- multiple_REs_variants %>%
  group_by(bg, nREs) %>% 
  count(name = "nvars") %>%
  ggplot(aes(x = as.factor(nREs), y = nvars, fill = bg)) +
  geom_col(width=0.9) +
  geom_text(aes(label=nvars), vjust=-0.5, color="black", size=4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_grid(rows = vars(bg), scales = "free") +
  scale_fill_manual(values=bg_color()) +
  labs(title = "",
       x = "Number of response elements bound", 
       y = "",
       fill = "Ancestral\nbackground") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = c(0.9, 0.9),
        strip.background = element_blank(),
        strip.text.y = element_blank())

a1 + a2 + plot_layout(widths = c(1, 5))

# plot as proportion promiscuous vs. specific
totals <- multiple_REs_variants %>%
  mutate(specific = factor(ifelse(nREs == 1, "specific", "promiscuous"),
                           levels = c("specific", "promiscuous"))) %>%
  group_by(bg, specific) %>%
  count(name = "nvars") %>%
  group_by(bg) %>%
  summarize(total = sum(nvars))
a3 <- multiple_REs_variants %>%
  mutate(specific = factor(ifelse(nREs == 1, "specific", "promiscuous"),
                           levels = c("specific", "promiscuous"))) %>%
  group_by(bg, specific) %>%
  count(name = "nvars") %>%
  ggplot(aes(x = bg, y = nvars, fill = bg, alpha = fct_rev(specific))) +
  geom_bar(stat = "identity", width=0.9, position = "fill") +
  geom_text(aes(label=nvars), position = position_fill(vjust = 0.5),
            color="black", size=4, alpha = 1) +
  # geom_text(data = totals, 
  #           aes(x = bg, label = total, y = 1, fill = NULL, alpha = NULL), 
  #           vjust = -0.5, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values=bg_color()) +
  scale_alpha_manual(values = c(0.5, 1)) +
  labs(title = "", x = "Ancestral\nbackground", 
       y = "Proportion of protein variants",
       fill = "Ancestral\nbackground",
       alpha = "Specificity") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_blank()) +
  guides(x = guide_axis(angle = 45))

# plot as proportion promiscuous
a4 <- multiple_REs_variants %>%
  mutate(specific = factor(ifelse(nREs == 1, "specific", "promiscuous"),
                           levels = c("specific", "promiscuous"))) %>%
  group_by(bg) %>%
  summarize(prop = sum(specific == "promiscuous") / n()) %>%
  ggplot(aes(x = bg, y = prop, fill = bg)) +
  geom_col(width=0.9) +
  scale_fill_manual(values=bg_color()) +
  labs(title = "", x = "Ancestral\nbackground", 
       y = "Fraction of promiscuous protein variants",
       fill = "Ancestral\nbackground",
       alpha = "Specificity") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_blank()) +
  guides(x = guide_axis(angle = 45))


# plot which REs each RH variant is active on by number of genotypes that it 
# binds to

# Find the RE variant(s) that is (are) bound by each protein variant. If a
# a protein variant binds to multiple sets of REs significantly but the union
# of those REs are not significantly bound, take the set with the highest
# p-value
multiple_REs_bound_data <- lapply(multiple_REs_bound_data, function(x)
  x %>% group_by(bg, AA_var) %>% arrange(desc(padj), .by_group = TRUE))
multiple_REs_variants$RE <- NA
for(i in 1:nrow(multiple_REs_variants)) {
  row <- multiple_REs_variants[i,]
  revars <- multiple_REs_bound_data[[as.numeric(row[3])]] %>% 
    ungroup() %>%
    semi_join(row, by = c("AA_var", "bg")) %>%
    slice(1) %>% 
    select(2:(as.numeric(row[3])+1)) %>% 
    t() %>% 
    str_c(collapse = ", ")
  multiple_REs_variants[i, "RE"] <- revars
}

# plot
b <- multiple_REs_variants %>%
  group_by(bg, nREs) %>%
  count(RE) %>%
  arrange(n) %>%
  ggplot(aes(x = factor(nREs), y = n, color = RE, fill = bg, label = RE)) +
  geom_col(width = 0.9, color = "black", position = "fill") +
  facet_grid(rows = vars(bg)) +
  geom_text(aes(color = bg), 
            position = position_fill(vjust = 0.5), size = 2.5) +
  geom_text(data = multiple_REs_variants %>% 
              group_by(bg, nREs) %>% count(name = "nvars"),
            aes(x = factor(nREs), label=nvars), y = 1, vjust=-0.5, color="black", size=4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values=bg_color()) +
  scale_color_manual(values = c("white", "black"), guide = "none") +
  labs(title = "",
       x = "Number of REs bound", 
       y = "Fraction of protein variants binding\neach RE combination",
       fill = "Ancestral\nbackground") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = c(0.9, 0.9),
        strip.background = element_blank(),
        strip.text.y = element_blank())
b


# for each RE, how many protein variants bind it specifically vs. promiscuously?
RE_frac_specific <- data.frame(bg = rep(NA, 32), RE = NA, 
                               specific = NA, promiscuous = NA)
i <- 1
for(revar in REs[[1]]) {
  for(background in c("AncSR1", "AncSR2")) {
    specific <- multiple_REs_variants %>%
      filter(bg == background, nREs == 1, RE == revar) %>%
      nrow()
    promiscuous <- multiple_REs_variants %>%
      filter(bg == background, nREs > 1, 
             grepl(sub("\\(..\\)", "", revar), RE)) %>%
      nrow()
    RE_frac_specific[i,] <- c(background, revar, specific, promiscuous)
    i <- i + 1
  }
}

c <- RE_frac_specific %>%
  pivot_longer(3:4, names_to = "specific", values_to = "n") %>%
  ggplot(aes(x = factor(RE, levels = levels(REs[[1]])), 
             y = as.numeric(n), 
             fill = bg, 
             alpha = factor(specific, levels = c("specific", "promiscuous")))) +
  geom_col(width = 0.9, position = "dodge") +
  facet_grid(rows = vars(bg), scales = "free") +
  scale_alpha_manual(values = c(1, 0.5)) +
  labs(x = "Response element", y = "Number of protein variants",
       fill = "Ancestral\nbackground", alpha = "Specificity") +
  scale_fill_manual(values=bg_color()) +
  theme_classic() +
  guides(x = guide_axis(angle = 45)) +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.y = element_blank())
c

# plot as proportion promiscuous vs. specific
totals <- RE_frac_specific %>%
  pivot_longer(3:4, names_to = "specific", values_to = "n") %>%
  mutate(n = as.numeric(n)) %>%
  group_by(bg, RE) %>%
  summarize(total = sum(n))
c2 <- RE_frac_specific %>%
  pivot_longer(3:4, names_to = "specific", values_to = "n") %>%
  mutate(n = as.numeric(n)) %>%
  ggplot(aes(x = factor(RE, levels = levels(REs[[1]])), 
             y = n, 
             fill = bg, 
             alpha = factor(specific, levels = c("promiscuous", "specific")))) +
  geom_bar(stat = "identity", width = 0.9, position = "fill") +
  # scale_y_continuous(expand = c(0, 0.1)) +
  facet_grid(rows = vars(bg), scales = "free") +
  scale_alpha_manual(values = c(0.5, 1)) +
  # geom_text(data = totals, 
  #           aes(x = RE, label = total, y = 1, fill = NULL, alpha = NULL), 
  #           vjust = -0.5) +
  labs(x = "Response element", y = "",
       fill = "Ancestral\nbackground", alpha = "Specificity") +
  scale_fill_manual(values=bg_color()) +
  theme_classic() +
  guides(x = guide_axis(angle = 45)) +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.y = element_blank())

# plot as proportion promiscuous
c3 <- RE_frac_specific %>%
  mutate(specific = as.numeric(specific),
         promiscuous = as.numeric(promiscuous),
         prop = promiscuous / (specific + promiscuous)) %>%
  ggplot(aes(x = factor(RE, levels = levels(REs[[1]])), y = prop, fill = bg)) +
  geom_col(width=0.9, position = "dodge") +
  scale_fill_manual(values=bg_color()) +
  labs(title = "", x = "Ancestral\nbackground", 
       y = "Fraction of promiscuous protein variants",
       fill = "Ancestral\nbackground",
       alpha = "Specificity") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_blank()) +
  guides(x = guide_axis(angle = 45))

a3 + plot_spacer() + c2 + plot_layout(widths = c(1.5, 0.5, 6))

# plot as number promiscuous vs. specific
c4 <- RE_frac_specific %>%
  pivot_longer(3:4, names_to = "specific", values_to = "n") %>%
  mutate(n = as.numeric(n)) %>%
  ggplot(aes(x = factor(RE, levels = levels(REs[[1]])), 
             y = n, 
             fill = bg, 
             alpha = factor(specific, levels = c("promiscuous", "specific")))) +
  geom_bar(stat = "identity", width = 0.9, position = "stack") +
  # scale_y_continuous(expand = c(0, 0.1)) +
  facet_grid(rows = vars(bg), scales = "free") +
  scale_alpha_manual(values = c(0.5, 1)) +
  # geom_text(data = totals, 
  #           aes(x = RE, label = total, y = 1, fill = NULL, alpha = NULL), 
  #           vjust = -0.5) +
  labs(x = "Response element", y = "Number of functional protein variants",
       fill = "Ancestral\nbackground", alpha = "Specificity") +
  scale_fill_manual(values=bg_color()) +
  theme_classic() +
  guides(x = guide_axis(angle = 45)) +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.y = element_blank())
c4

# plot as number promiscuous vs. specific side by side
c5 <- RE_frac_specific %>%
  pivot_longer(3:4, names_to = "specific", values_to = "n") %>%
  mutate(n = as.numeric(n)) %>%
  ggplot(aes(x = bg, 
             y = n, 
             fill = bg, 
             alpha = factor(specific, levels = c("promiscuous", "specific")))) +
  geom_bar(stat = "identity", width = 0.9, position = "stack") +
  facet_grid(cols = vars(factor(RE, levels = levels(REs[[1]])))) +
  # scale_y_continuous(expand = c(0, 0.1)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  # geom_text(data = totals, 
  #           aes(x = RE, label = total, y = 1, fill = NULL, alpha = NULL), 
  #           vjust = -0.5) +
  labs(x = "Response element", y = "Number of functional protein variants",
       fill = "Ancestral\nbackground", alpha = "Specificity") +
  scale_fill_manual(values=bg_color()) +
  theme_classic() +
  guides(x = guide_axis(angle = 45)) +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.x = element_blank())
c5
```

## What is the distribution of fluorescence for protein variants that bind specifically vs. promiscously?

``` r
meanF_data_functional <- meanF_data_functional %>%
  left_join(multiple_REs_variants %>% select(-RE)) %>%
  mutate(nREs = ifelse(nREs == 1, "specific", "promiscuous")) %>% 
  rename(specific = nREs)
a <- meanF_data_functional %>%
  drop_na(specific) %>%
  mutate(specific = factor(specific, levels = c("specific", "promiscuous"))) %>%
  ggplot(aes(x = specific, 
             y = avg_meanF, fill = bg, alpha = specific)) +
  geom_boxplot() +
  scale_alpha_manual(values = c(1, 0.5), guide = "none") +
  labs(y = "Fluorescence", fill = "Ancestral\nbackground", x = "") +
  scale_fill_manual(values=bg_color()) +
  geom_signif(aes(group = specific),
              comparisons = list(c("specific", "promiscuous")), 
              test = "wilcox.test",
              map_signif_level = TRUE) +
  facet_grid(cols = vars(bg)) +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.x = element_blank())
a
```

``` r
a2 <- multiple_REs_variants %>%
  group_by(bg, nREs) %>% 
  count(name = "nvars") %>%
  ggplot(aes(x = as.factor(nREs), y = nvars, fill = bg)) +
  geom_col(width=0.9) +
  geom_text(aes(label=nvars), vjust=-0.5, color="black", size=4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_grid(rows = vars(bg), scales = "free") +
  scale_fill_manual(values=bg_color()) +
  labs(title = "",
       x = "Number of response elements bound", 
       y = "Number of protein variants",
       fill = "Ancestral\nbackground") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = c(0.9, 0.9),
        strip.background = element_blank(),
        strip.text.y = element_blank())
a2

c <- fluorescence_both_backgrounds %>%
  pivot_longer(3:4, names_to = "background") %>%
  ggplot(aes(x = background, y = value, fill = background)) +
  geom_boxplot() +
  scale_fill_manual(values = c(bg_color()), drop = FALSE) +
  labs(title = "Protein:RE variant combos\nfunctional in both backgrounds", 
       x = "Ancestral background", y = "Fluorescence", 
       fill = "Ancestral\nbackground") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line.x = element_line(color = "black",),
        legend.position = "none",
        plot.caption = element_text(hjust = 0)) +
  guides(x = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(c("AncSR1", "AncSR2")), 
              test = "wilcox.test",
              test.args = list(alternative = "less"),
              map_signif_level = function(p) sprintf("p = %.2g", p))
c

c + a2 + plot_layout(widths = c(1.5, 3))
```

## What is the genetic basis of specificity vs. promiscuity?

``` r
# logo plots for all functional variants in the AncSR2 background
a <- ggseqlogo(meanF_data_functional %>% 
                 filter(bg == "AncSR2") %>%
                 pull(AA_var) %>% unique() %>% as.character(),
                 method = "bits") +
  theme(text=element_text(size=fontsize),
        plot.title = element_text(size = fontsize),
        legend.position = "none") +
  labs(x="AA position", title="All AncSR2 protein variants")
a

# logo plots for variants that bind specifically to one RE in the AncSR2 background
AncSR2_specific_vars <- lapply(REs[[1]], function(x) multiple_REs_variants %>%
                                 filter(bg == "AncSR2", nREs == 1, RE == x) %>% 
                                 pull(AA_var) %>%
                                 as.character())
names(AncSR2_specific_vars) <- REs[[1]]
AncSR2_specific_vars <- AncSR2_specific_vars[sapply(AncSR2_specific_vars, length) >= 50]
b <- ggseqlogo(AncSR2_specific_vars, 
          method="bits", ncol = length(meanF_data_functional)) +
  theme(text=element_text(size=fontsize),
        strip.text.x = element_text(size = fontsize),
        legend.position = "none") +
  labs(x="AA position")
b


# logo plots for variants that bind promiscuously to more than one RE in the AncSR2 background
c <- ggseqlogo(meanF_data_functional %>% 
                 filter(specific == "promiscuous", bg == "AncSR2") %>%
                 pull(AA_var) %>% unique() %>% as.character(), 
          method="bits") +
  theme(text=element_text(size=fontsize),
        legend.position = "none") +
  labs(x="AA position", title="Promiscuous AncSR2 protein variants")
c
```

## do variants switch between functional/nonfunctional/promiscuous/specific between AncSR1 and AnSR2?

``` r
alluvial_data <- multiple_REs_variants %>%
  mutate(specificity = factor(ifelse(nREs > 1, "promiscuous", RE), 
                           levels = c(levels(REs[[1]]), 
                                      "same specificity", "switched specificity", 
                                      "promiscuous", "nonfunctional"))) %>%
  select(AA_var, bg, specificity) %>%
  pivot_wider(names_from = bg, values_from = specificity) %>%
  replace_na(list(AncSR1 = "nonfunctional", AncSR2 = "nonfunctional")) %>%
  filter(AncSR1 != "nonfunctional")
for(i in 1:nrow(alluvial_data)) {
  spec1 <- as.character(alluvial_data[i,2] %>% pull(1))
  spec2 <- as.character(alluvial_data[i,3] %>% pull(1))
  if(spec2 %in% as.character(REs[[1]])) {
    if(spec1 %in% as.character(REs[[1]])) {
      if(spec1 == spec2) alluvial_data[i,3] <- "same specificity"
      else alluvial_data[i,3] <- "switched specificity"
    } else {
      alluvial_data[i,3] <- "switched specificity"
    }
  }
}

alluvial_data <- alluvial_data %>%
  mutate(AncSR1 = factor(ifelse(AncSR1 %in% REs[[1]], "specific", "promiscuous"),
                         levels = c("promiscuous", "specific")),
         AncSR2 = factor(AncSR2, levels = c("promiscuous", "same specificity", 
                                            "switched specificity", "nonfunctional")))


a <- alluvial_data %>%
  filter(AncSR2 != "nonfunctional") %>%
  count(AncSR1, AncSR2) %>%
  ggplot(aes(axis1 = AncSR1, 
             axis2 = AncSR2, 
             y = n)) +
  geom_alluvium(aes(fill = AncSR2), width = 3/5) +
  geom_stratum(width = 3/5) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("AncSR1\nbackground", "AncSR2\nbackground"),
                   expand = c(0.15, 0.05),
                   position = "top") +
  # scale_fill_manual(values = viridis(4)) +
  labs(title = "Protein variants functional\nin both backgrounds", 
       x = "", y = "Number of protein variants",
       fill = "AncSR2 phenotype") +
  theme_classic() +
  guides(fill = "none") +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = fontsize))
a


# mean fluorescence for things that became promiscuous, 
# retained specificity, or retained promiscuity
b <- alluvial_data %>%
  filter(AncSR2 %in% c("promiscuous", "same specificity")) %>%
  unite(col = "phechange", AncSR1, AncSR2, sep = " to ") %>%
  mutate(phechange = factor(ifelse(phechange == "specific to same specificity",
                            "retained specificity", 
                            ifelse(phechange == "promiscuous to promiscuous",
                                   "retained promiscuity",
                                   phechange)),
                            levels = c("retained specificity",
                                       "specific to promiscuous",
                                   "retained promiscuity" ))) %>%
  inner_join(meanF_data_functional) %>%
  filter(bg == "AncSR1") %>%
  ggplot(aes(x = phechange, y = avg_meanF, fill = bg)) +
  geom_boxplot() +
  scale_fill_manual(values = bg_color("AncSR1")) +
  labs(x = "AncSR1 to AncSR2\nphenotype change", y = "Fluorescence", 
       fill = "", title = "Protein:RE variant combos\nfunctional in AncSR1 background") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line.x = element_line(color = "black",),
        legend.position = "none",
        plot.caption = element_text(hjust = 0)) +
  guides(x = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(c("retained specificity", "specific to promiscuous"),
                                 c("specific to promiscuous", "retained promiscuity")), 
              test = "wilcox.test",
              test.args = list(alternative = "less"),
              map_signif_level = function(p) sprintf("p = %.2g", p))

c <- fluorescence_both_backgrounds %>%
  pivot_longer(3:4, names_to = "background") %>%
  ggplot(aes(x = background, y = value, fill = background)) +
  geom_boxplot() +
  scale_fill_manual(values = c(bg_color()), drop = FALSE) +
  labs(title = "Protein:RE variant combos\nfunctional in both backgrounds", 
       x = "Ancestral background", y = "", fill = "Ancestral\nbackground") +
  theme_classic() +
  theme(text = element_text(size = fontsize),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line.x = element_line(color = "black",),
        legend.position = "right",
        plot.caption = element_text(hjust = 0)) +
  guides(x = guide_axis(angle = 45)) +
  geom_signif(comparisons = list(c("AncSR1", "AncSR2")), 
              test = "wilcox.test",
              test.args = list(alternative = "less"),
              map_signif_level = function(p) sprintf("p = %.2g", p))

b + c + plot_layout(widths = c(3, 3))
```

# Network analysis

Let’s create the genotype networks, where each node is a protein variant
and their associated traits are the REs that they bind to.

``` r
# create data frame of which protein variants bind to which REs
REs_bound <- multiple_REs_variants %>%
  separate_longer_delim(cols = RE, delim = ", ") %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = RE, values_from = present)
  

## AncSR1
AncSR1_comb <- combn(REs_bound %>% 
                       filter(bg == "AncSR1") %>% 
                       pull(AA_var) %>% 
                       as.character(), 2)
AncSR1_edges <- t(AncSR1_comb)
AncSR1_connected <- apply(AncSR1_edges, 1, connected)
AncSR1_connected <- AncSR1_edges[AncSR1_connected,]
  

# create undirected graph where edges represent single step mutations between
# RH genotypes permissible by the genetic code; node attributes are boolean
# variables indicating binding to each RE, and one integer variable for number
# of REs bound
AncSR1_graph <- graph_from_data_frame(AncSR1_connected, directed=FALSE, 
                                      vertices = REs_bound %>% 
                                        filter(bg == "AncSR1") %>% select(-bg) %>%
                                        mutate(AA_var = as.character(AA_var)))
write_graph(AncSR1_graph, 
            file.path(results_dir, "AncSR1_RH_network.graphml"),
            format="graphml")

## AncSR2
AncSR2_comb <- combn(REs_bound %>% 
                       filter(bg == "AncSR2") %>% 
                       pull(AA_var) %>% 
                       as.character(), 2)
AncSR2_edges <- t(AncSR2_comb)
AncSR2_connected <- apply(AncSR2_edges, 1, connected)
AncSR2_connected <- AncSR2_edges[AncSR2_connected,]
  

# create undirected graph where edges represent single step mutations between
# RH genotypes permissible by the genetic code; node attributes are boolean
# variables indicating binding to each RE, and one integer variable for number
# of REs bound
AncSR2_graph <- graph_from_data_frame(AncSR2_connected, directed=FALSE, 
                                      vertices = REs_bound %>% 
                                        filter(bg == "AncSR2") %>% select(-bg) %>%
                                        mutate(AA_var = as.character(AA_var)))
write_graph(AncSR2_graph, 
            file.path(results_dir, "AncSR2_RH_network.graphml"),
            format="graphml")
```