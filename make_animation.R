library(tidyverse)
library(ggtext)
library(animation)
library(tweenr)


# Read data
source("tidy_data.R")


# Prepare data for tween'ed animation
data <- df %>%
    pivot_wider(names_from = Variable, values_from = Value) %>%
    arrange(name, Concentration) %>%
    select(starts_with("average"), everything()) %>%
    as.data.frame()
data_tween <- filter(data, Concentration == "0.04um") %>%
    keep_state(20) %>%
    tween_state(
        filter(data, Concentration == "0.21um"),
        ease = "cubic-in-out", nframes = 100) %>%
    keep_state(20) %>%
    tween_state(
        filter(data, Concentration == "1.04um"),
        ease = "cubic-in-out", nframes = 100) %>%
    keep_state(20) %>%
    tween_state(
        filter(data, Concentration == "5um"),
        ease = "cubic-in-out", nframes = 100) %>%
    keep_state(20) %>%
    group_split(.frame)


# Make tween'ed animation
oopt <- ani.options(interval = 1 / 20)
i <- 1
setwd("01_plots")
saveGIF({
    for (df in data_tween) {
        cat(sprintf("Processing %i/%i\n", i, length(data_tween)))
        gg <- ggplot(df, aes(average.p53.intensity, average.cell.count)) +
            geom_point(size = 3, alpha = 0.2) +
            geom_point(
                data = df %>% filter(str_detect(Drug.Name, "Flavopir")),
                aes(average.p53.intensity, average.cell.count),
                colour = "#498eaf", alpha = 0.8, size = 4) +
            geom_point(
                data = df %>% filter(str_detect(Drug.Name, "Dinaci")),
                aes(average.p53.intensity, average.cell.count),
                colour = "#d44e28", alpha = 0.8, size = 4) +
            geom_point(
                data = df %>% filter(str_detect(Main.Target, "CDK")),
                aes(average.p53.intensity, average.cell.count),
                colour = "#e5bb4b", size = 6, shape = 1) +
            theme_minimal() +
            facet_wrap(
                ~ Normalisation,
                ncol = 1,
                labeller = as_labeller(norm_name)) +
            xlim(-4, 5) + ylim(-6, 3) +
            labs(
                x = "Log2 p53 intensity", y = "Log2 cell count",
                title = "p53 intensity vs. cell count at different drug doses",
                subtitle = sprintf(
                    toString(c(
                        "Dose = %s",
                        "<span style='color:#498eaf'>Flavopiridol</span>",
                        "<span style='color:#d44e28'>Dinaciclib</span>",
                        "<span style='color:#e5bb4b'>CDK inhibitor</span>")),
                    unique(df$Concentration))) +
            theme(
                text = element_text(size = 28),
                plot.subtitle = element_markdown())
        plot(gg)
        i <- i + 1
    }},
    movie.name = "animation.gif",
    ani.width = 1024, ani.height = 1024)
