library(tidyverse)
library(ggtext)
library(ggrepel)


# Read data
source("tidy_data.R")


# Calculate Mahalanobis distances by Compound, Concentration and Quadrant;
# only consider data normalised to DMSO-treated RPS19-KD

df_res <- df %>%
    pivot_wider(names_from = Variable, values_from = Value) %>%
	filter(Normalisation == "RPS19") %>%
	select(-Normalisation) %>%
    group_by(Concentration) %>%
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
    group_by(Concentration, Quadrant) %>%
    mutate(Rank = rank(-Mahalanobis_dist, ties.method = "min")) %>%
	ungroup()


# Combine Mahalanobis distances for every compound by calculating the 2-norm,
# rank compounds by decreasing 2-norm and save as CSV table
df_final <- df_res %>%
	select(name, Quadrant, Concentration, Rank) %>%
	complete(name, Quadrant, Concentration) %>%
	group_by(Quadrant) %>%
	mutate(Rank = replace_na(Rank, max(Rank, na.rm = TRUE))) %>%
	group_by(name, Quadrant) %>%
	summarise(
		n_doses = length(Rank),
		Rank_at_dose = toString(paste(Rank, Concentration, sep = "@")),
		Rank_product_metric = prod(Rank)^(1/n_doses)) %>%
	group_by(Quadrant) %>%
	mutate(Combined_rank = rank(Rank_product_metric, ties.method = "min")) %>%
	ungroup()


# Store as Excel sheet
dir.create(file.path(".", "02_tables"), showWarnings = FALSE)
dir.create(file.path(".", "02_tables", "final_ranking"), showWarnings = FALSE)
wd <- getwd()
df_final %>%
	left_join(
		df_res %>% distinct(name, Drug.Name, Main.Target),
		by = "name") %>%
	group_by(Quadrant) %>%
	select(-n_doses) %>%
	nest() %>%
	mutate(data = map(
		data,
		~ .x %>%
			select(Combined_rank, everything()) %>%
			arrange(Combined_rank))) %>%
	mutate(tmp = pmap(
		list(Quadrant, data),
		function(quad, df) {
			setwd("./02_tables/final_ranking")
			fn <- sprintf("final_%s.csv", quad)
			write_csv(df, fn)
            setwd(wd)
		}))


class_name <- c(
	"Protein", "Transcription - Pol I and II",
	"Transcription - Pol II", "Other")
data <- df_res %>%
	left_join(df_final, by = c("name", "Quadrant")) %>%
	mutate(Class = case_when(
		str_detect(Main.Target, "HSP90") ~ "Protein",
		str_detect(Main.Target, "opoisomerase") ~ "Transcription - Pol I and II",
		str_detect(Drug.Name, "(Flavopiridol|Dinaciclib)") ~ "Transcription - Pol II",
		TRUE ~ "Other")) %>%
	mutate(Class = factor(Class, levels = class_name))
y_limits <- c(1, NA)
ggplot(data, aes(average.p53.intensity, average.cell.count)) +
	geom_point() +
	facet_wrap(~ Concentration) +
	theme_minimal() +
	geom_point(
		data = data %>%
			filter(Combined_rank <= 30 & Quadrant == "p53_up.cellcount_down"),
		aes(average.p53.intensity, average.cell.count, colour = Class),
		size = 3, alpha = 0.5) +
	geom_text_repel(
		data = data %>%
			filter(Combined_rank <= 30 & Quadrant == "p53_up.cellcount_down") %>%
			mutate(label = if_else(
				str_detect(Drug.Name, "(Flavopiridol|Dinaciclib)"),
				Drug.Name, "")),
		aes(average.p53.intensity, average.cell.count, label = label),
		size = 2, colour = "red", force = 10, ylim = y_limits) +
	scale_colour_manual(values = setNames(
		c("blue", "orange", "red", "grey"),
		class_name))
ggsave(
	"01_plots/p53_intensity_vs_cell_count_final_ranking_highlight.png",
	height = 10,
	width = 10)
ggsave(
	"01_plots/p53_intensity_vs_cell_count_final_ranking_highlight.pdf",
	height = 10,
	width = 10)
