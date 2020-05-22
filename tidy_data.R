library(tidyverse)

# Screening data
df_data <- read_csv(
    "00_original_data/data_merge_RPS19_compoundscreen.csv")


# If there are multiple entries by unique `Compound.Alias` (see e.g.
# "SN00750429") we keep only the first row.
df_compound_details <- read_csv(
    "00_original_data/WEHI FDA LIbrary Compound Details_WS2.csv") %>%
    group_by(Compound.Alias) %>%
    filter(row_number() == 1)


# Left-join of data with compound details based on `Compound.Alias` and
# reshape from wide to long
df <- df_data %>%
    left_join(df_compound_details, by = c("name" = "Compound.Alias")) %>%
    select(name, starts_with("average"), Drug.Name, Main.Target) %>%
    pivot_longer(
        starts_with("average"),
        names_to = c(".value", "Concentration"),
        names_pattern = "(average.+[9T])\\.(.+um)") %>%
    mutate_if(is.character, iconv, from = "latin1", to = "UTF-8") %>%
    mutate_at(vars(starts_with("average")), log2) %>%
    pivot_longer(
        starts_with("average"),
        names_to = c(".value", "Normalisation"),
        names_pattern = "(.+).normalised.to.(.+)") %>%
	pivot_longer(
		starts_with("average"),
		names_to = "Variable",
		values_to = "Value")


# Set more meaningful names for the two types of normalisation
# These names will be used in plots
norm_name <- c(
    "OTP.NT" = "normalised to non-targeting control",
    "RPS19" = "normalised to RPS19-KD")
