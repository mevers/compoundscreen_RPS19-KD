library(tidyverse)
library(ggtext)


# Read data
source("tidy_data.R")


# Plot
dir.create(file.path(".", "01_plots"), showWarnings = FALSE)
data <- df %>%
    pivot_wider(names_from = Variable, values_from = Value)
ggplot(data, aes(average.p53.intensity, average.cell.count)) +
    geom_point(alpha = 0.2) +
    geom_point(
        data = data %>% filter(str_detect(Drug.Name, "Flavopir")),
        aes(average.p53.intensity, average.cell.count),
        colour = "#498eaf", alpha = 0.6, size = 3) +
    geom_point(
        data = data %>% filter(str_detect(Drug.Name, "Dinaci")),
        aes(average.p53.intensity, average.cell.count),
        colour = "#d44e28", alpha = 0.6, size = 3) +
    facet_grid(
        Normalisation ~ Concentration,
        labeller = labeller(Normalisation = as_labeller(norm_name))) +
    labs(
        x = "Log2 p53 intensity",
        y = "Log2 cell count",
        title = strwrap(
            "p53 intensity vs. cell count at different drug doses and using
            different normalisation strategies",
            width = 10000, simplify = TRUE),
        subtitle =
            toString(c(
                "Highlighted: <span style='color:#498eaf'>Flavopiridol</span>",
                "<span style='color:#d44e28'>Dinaciclib</span>"))) +
    theme_minimal() +
    theme(plot.subtitle = element_markdown())
ggsave("01_plots/p53_intensity_vs_cell_count.pdf", height = 6, width = 10)
ggsave("01_plots/p53_intensity_vs_cell_count.png", height = 6, width = 10)
