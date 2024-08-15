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


# Make long
long_data1 <- data1 %>%
  melt(id.vars = "time", variable.name = "label", value.name = "OD600") %>%
  separate(label, into = c("strain", "replicate"), sep = "--") %>%
  mutate(minutes = round(time * 60, 0))

long_plate_reader_data <- plate_reader_data %>%
  melt(id.vars = "time", variable.name = "label", value.name = "OD600") %>%
  separate(label, into = c("strain", "replicate"), sep = "--") %>%
  mutate(minutes = round(time * 60, 0))

ggplot(
  long_data1,
  aes(x = minutes, y = OD600, group = strain)
) +
  # scale_y_log10() +
  stat_summary(fun = mean, geom = "line", aes(color = factor(strain))) +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = factor(strain)), alpha = 0.4) +
  theme_minimal()

ggplot(
  long_plate_reader_data,
  aes(x = minutes, y = OD600, group = strain)
) +
  # scale_y_log10() +
  stat_summary(fun = mean, geom = "line", aes(color = factor(strain))) +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = factor(strain)), alpha = 0.4) +
  theme_minimal()
