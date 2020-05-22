library(tidyverse)
library(ggtext)


# Read data
source("tidy_data.R")


# Calculate Mahalanobis distances and determine rank by Normalisation,
# Concentration and Quadrant
df_res <- df %>%
    pivot_wider(names_from = Variable, values_from = Value) %>%
    group_by(Normalisation, Concentration) %>%
    nest() %>%
    mutate(data = map(data, function(tbl) {
        mat <- tbl %>% select(starts_with("average")) %>% as.matrix()
        tbl %>% mutate(
            Mahalanobis_dist = mahalanobis(mat, colMeans(mat), cov(mat)),
            Quadrant = case_when(
                average.p53.intensity > 0 & average.cell.count > 0 ~ "p53_up.cellcount_up",
                average.p53.intensity < 0 & average.cell.count > 0 ~ "p53_down.cellcount_up",
                average.p53.intensity < 0 & average.cell.count < 0 ~ "p53_down.cellcount_down",
                TRUE ~ "p53_up.cellcount_down"))
        })) %>%
    unnest(data) %>%
    group_by(Normalisation, Concentration, Quadrant) %>%
    mutate(Rank = rank(-Mahalanobis_dist, ties.method = "first")) %>%
    ungroup()


# Plot Mahalanobis distances
dir.create(file.path(".", "01_plots"), showWarnings = FALSE)
quad_name = c(
    "p53_up.cellcount_up" = "Increased p53 and\nincreased cell count",
    "p53_down.cellcount_up" = "Decreased p53 and\nincreased cell count",
    "p53_down.cellcount_down" = "Decreased p53 and\ndecreased cell count",
    "p53_up.cellcount_down" = "Increased p53 and\ndecreased cell count")
df_res %>%
    mutate(Quadrant = factor(Quadrant, levels = names(quad_name))) %>%
    filter(Normalisation == "RPS19") %>%
    ggplot(aes(Mahalanobis_dist)) +
    geom_histogram(bins = 50, position = "nudge") +
    facet_grid(
        Quadrant ~ Concentration,
        labeller = labeller(Quadrant = as_labeller(quad_name))) +
    scale_x_log10() +
    theme_minimal() +
    labs(
        title = strwrap(
            "Distribution of Mahalanobis distances at different drug doses
            relative to DMSO-treated RPS19-KD",
            width = 10000,
            simplify = TRUE),
        x = "Mahalanobis distance", y = "Count")
ggsave("01_plots/distr_Mahalanobis_dist.pdf", height = 6, width = 10)
ggsave("01_plots/distr_Mahalanobis_dist.png", height = 6, width = 10)


# Store as Excel sheets
dir.create(file.path(".", "02_tables"), showWarnings = FALSE)
walk(names(norm_name), ~dir.create(
    file.path(".", "02_tables", .x), showWarnings = FALSE))
wd <- getwd()
df_res %>%
    group_by(Normalisation, Concentration, Quadrant) %>%
    nest() %>%
    filter(str_detect(Quadrant, "p53_up")) %>%
    mutate(data = map(
        data, ~.x %>% arrange(Rank) %>% select(Rank, everything()))) %>%
    mutate(tmp = pmap(
        list(Normalisation, Concentration, Quadrant, data),
        function(norm, conc, quad, df) {
            setwd(sprintf("./02_tables/%s", norm))
            fn <- sprintf("%s.csv", paste(quad, conc, sep = "_"))
            write_csv(df, fn)
            setwd(wd)
        }
    ))


df_stan <- df %>% filter(Concentration == "1.04um")
stan_model <- stan_model("model.stan")
fit <- sampling(
	stan_model,
	data = list(
		N = nrow(df_stan),
		y = df_stan %>% select(contains("RPS19"))))


df_ppc <- fit %>%
	summary("y_rep") %>%
	pluck("summary") %>%
	as_tibble() %>%
	select(mean) %>%
	mutate(
		idx = rep(1:nrow(df_stan), each = 2),
		col = rep(
			paste(
				"PPC",
				df_stan %>% select(contains("RPS19")) %>% names(),
				sep = "_"),
			nrow(df_stan))) %>%
	pivot_wider(names_from = "col", values_from = "mean")

df_stan %>%
	rename_at(vars(contains("RPS19")), ~paste("Data", .x, sep = "_")) %>%
	bind_cols(df_ppc) %>%
	select(contains("RPS19")) %>%
	pivot_longer(
		everything(),
		names_to = c(".value", "Variable"),
		names_pattern = "(PPC|Data)_*(average.*$)") %>%
	ggplot(aes(Data, PPC)) +
	geom_point() +
	facet_wrap(~ Variable) +
	theme_minimal()
