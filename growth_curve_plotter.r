library(pacman)

p_load(data.table, tidyverse)

# Custom function to make names unique by appending numbers starting from 1
make.unique.reps <- function(names) {
  counts <- table(names) # Count occurrences of each name
  new_names <- names # Start with original names
  for (name in names(counts)) {
    if (counts[name] > 1) { # If name appears more than once
      indices <- which(names == name)
      for (i in seq_along(indices)) {
        new_names[indices[i]] <- paste0(name, "--", i) # Append suffix
      }
    }
  }
  new_names
}
data1 <- fread("growth_curves/Data 1.txt")

# Ensure column names are unique
colnames(data1) <- make.unique.reps(colnames(data1))

plate_reader_data <- fread("growth_curves/plate reader data.txt")

# Ensure column names are unique
colnames(plate_reader_data) <- make.unique.reps(colnames(plate_reader_data))

# Now rename V1 to time
data1 <- data1 %>% rename(time = V1)

plate_reader_data <- plate_reader_data %>% rename(time = V1)

# define color for each strain
# ispD = red,
# Non-targeting Control = blue,
# Wild Type = black

color_map <- c(
  "ispD" = "#e31a1c",
  "pgsA Hypomorph" = "#6a3d9a",
  "Non-targeting Control" = "#33a02c",
  "Non-targeting (mRFP) Control" = "#33a02c",
  "Wild Type" = "black"
)

strain_order <- c("Wild Type", "Non-targeting Control", "Non-targeting (mRFP) Control", "pgsA Hypomorph", "ispD")
# Make long

long_data1 <- data1 %>%
  melt(id.vars = "time", variable.name = "label", value.name = "OD600") %>%
  separate(label, into = c("strain", "replicate"), sep = "--") %>%
  mutate(minutes = round(time * 60, 0)) %>%
  mutate(hours = minutes / 60)

long_plate_reader_data <- plate_reader_data %>%
  melt(id.vars = "time", variable.name = "label", value.name = "OD600") %>%
  separate(label, into = c("strain", "replicate"), sep = "--") %>%
  mutate(minutes = round(time * 60, 0)) %>%
  mutate(hours = minutes / 60)

ggplot(
  long_data1 %>% filter(strain != "ispG") %>%
    mutate(strain = case_when(
      strain == "aRFP control" ~ "Non-targeting (mRFP) Control",
      strain == "WT" ~ "Wild Type",
      TRUE ~ strain
    )) %>%
    mutate(strain = factor(strain, levels = strain_order)),
  aes(x = hours, y = OD600, group = strain)
) +
  stat_summary(fun = mean, geom = "line", aes(color = factor(strain))) +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = factor(strain)), alpha = 0.3) +
  theme_minimal() +
  labs(x = "Time (hours)", y = expression("OD"[600]), color = "Strain", fill = "Strain") +
  theme(
    legend.title = element_text(face = "bold", size = 12, color = "black", angle = 0),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map)


ggplot(
  long_plate_reader_data %>% mutate(strain = case_when(
    strain == "5-P1" ~ "pgsA Hypomorph",
    strain == "Mci Control" ~ "Non-targeting Control",
    TRUE ~ NA_character_
  )) %>% filter(!is.na(strain)) %>% mutate(strain = factor(strain, levels = strain_order)),
  aes(x = hours, y = OD600, group = strain)
) +
  stat_summary(fun = mean, geom = "line", aes(color = factor(strain))) +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = factor(strain)), alpha = 0.4) +
  theme_minimal() +
  labs(x = "Time (hours)", y = expression("OD"[600]), color = "Strain", fill = "Strain") +
  theme(
    legend.title = element_text(face = "bold", size = 12, color = "black", angle = 0),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map)
