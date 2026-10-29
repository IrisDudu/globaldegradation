#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript plot_theilsen_sdg1531_forecasts.R <input.csv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_columns <- c("SOC", "year", "degraded_area_km2", "improved_area_km2")
raw <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(required_columns %in% names(raw)))
stopifnot(identical(sort(unique(raw$year)), 2016:2024))
stopifnot(!anyNA(raw[, c("degraded_area_km2", "improved_area_km2")]))

region_key <- tibble::tribble(
  ~region,          ~soc,   ~file_stem,
  "China",         "CHN", "china",
  "Russia",        "RUS", "russia",
  "United States", "USA", "united_states",
  "Canada",        "CAN", "canada",
  "Brazil",        "BRA", "brazil",
  "Australia",     "AUS", "australia",
  "India",         "IND", "india"
)

sum_areas <- function(data) {
  data %>%
    group_by(year) %>%
    summarise(
      degradation = sum(degraded_area_km2, na.rm = TRUE) / 1e6,
      improvement = sum(improved_area_km2, na.rm = TRUE) / 1e6,
      .groups = "drop"
    )
}

observed <- bind_rows(
  sum_areas(raw) %>% mutate(region = "Global", file_stem = "global"),
  bind_rows(lapply(seq_len(nrow(region_key)), function(i) {
    key <- region_key[i, ]
    selected <- raw %>% filter(SOC == key$soc)
    if (nrow(selected) == 0L) stop("No records found for ", key$region)
    sum_areas(selected) %>% mutate(region = key$region, file_stem = key$file_stem)
  }))
) %>%
  select(region, file_stem, year, degradation, improvement)

expected_regions <- c("Global", region_key$region)
check <- observed %>% count(region)
stopifnot(setequal(check$region, expected_regions), all(check$n == 9L))

theil_sen <- function(year, area) {
  x <- year - 2020
  pair_slopes <- numeric(0)
  for (i in seq_len(length(x) - 1L)) {
    for (j in (i + 1L):length(x)) {
      pair_slopes <- c(pair_slopes, (area[j] - area[i]) / (x[j] - x[i]))
    }
  }
  slope <- median(pair_slopes, na.rm = TRUE)
  intercept_2020 <- median(area - slope * x, na.rm = TRUE)
  list(slope = slope, intercept_2020 = intercept_2020)
}

logistic_target_path <- function(years, base_value, target_value,
                                 base_year, target_year, steepness) {
  target_fraction <- target_value / base_value
  delta <- target_year - base_year
  q <- (1 - target_fraction) /
    (target_fraction * exp(steepness * delta) - 1)
  base_value * (1 + q) /
    (1 + q * exp(steepness * (years - base_year)))
}

model_rows <- list()
forecast_rows <- list()
target_rows <- list()
plot_rows <- list()

series_levels <- c(
  "Degradation obs.",
  "Degradation forecast",
  "Improvement obs.",
  "Improvement forecast",
  "KM-GBF target-consistent\nreference pathways",
  "G20 target-consistent\nreference pathways"
)

for (region_i in expected_regions) {
  z <- observed %>% filter(region == region_i) %>% arrange(year)
  stem_i <- unique(z$file_stem)

  deg_fit <- theil_sen(z$year, z$degradation)
  imp_fit <- theil_sen(z$year, z$improvement)

  model_rows[[length(model_rows) + 1L]] <- tibble(
    region = region_i,
    series = c("degradation", "improvement"),
    slope_million_km2_per_year = c(deg_fit$slope, imp_fit$slope),
    intercept_at_2020_million_km2 = c(deg_fit$intercept_2020, imp_fit$intercept_2020)
  )

  future_years <- 2025:2050
  fc <- tibble(
    region = region_i,
    year = future_years,
    degradation_forecast = pmax(0, deg_fit$intercept_2020 + deg_fit$slope * (future_years - 2020)),
    improvement_forecast = pmax(0, imp_fit$intercept_2020 + imp_fit$slope * (future_years - 2020))
  )
  forecast_rows[[length(forecast_rows) + 1L]] <- fc

  deg_2022 <- z$degradation[z$year == 2022]
  deg_2024 <- z$degradation[z$year == 2024]

  km_years <- 2022:2050
  g20_years <- 2024:2050
  km_target <- 0.70 * deg_2022
  g20_target <- 0.50 * deg_2024

  targets <- bind_rows(
    tibble(
      region = region_i,
      pathway = "KM-GBF target-consistent\nreference pathways",
      year = km_years,
      area = logistic_target_path(
        km_years, deg_2022, km_target,
        base_year = 2022, target_year = 2030, steepness = 0.18
      ),
      target_year = 2030,
      target_area = km_target
    ),
    tibble(
      region = region_i,
      pathway = "G20 target-consistent\nreference pathways",
      year = g20_years,
      area = logistic_target_path(
        g20_years, deg_2024, g20_target,
        base_year = 2024, target_year = 2040, steepness = 0.28
      ),
      target_year = 2040,
      target_area = g20_target
    )
  )
  target_rows[[length(target_rows) + 1L]] <- targets

  last_obs <- z %>% filter(year == 2024)
  plot_data <- bind_rows(
    z %>% transmute(year, area = degradation, series = "Degradation obs."),
    bind_rows(
      last_obs %>% transmute(year, area = degradation, series = "Degradation forecast"),
      fc %>% transmute(year, area = degradation_forecast, series = "Degradation forecast")
    ),
    z %>% transmute(year, area = improvement, series = "Improvement obs."),
    bind_rows(
      last_obs %>% transmute(year, area = improvement, series = "Improvement forecast"),
      fc %>% transmute(year, area = improvement_forecast, series = "Improvement forecast")
    ),
    targets %>% transmute(year, area, series = pathway)
  ) %>%
    mutate(series = factor(series, levels = series_levels))

  target_points <- targets %>%
    distinct(pathway, target_year, target_area) %>%
    transmute(
      year = target_year,
      area = target_area,
      series = factor(pathway, levels = series_levels)
    )

  plot_rows[[length(plot_rows) + 1L]] <- plot_data %>% mutate(region = region_i)

  colours <- c(
    "Degradation obs." = "#171717",
    "Degradation forecast" = "#D94F50",
    "Improvement obs." = "#009E5B",
    "Improvement forecast" = "#6BCB4A",
    "KM-GBF target-consistent\nreference pathways" = "#8A8F70",
    "G20 target-consistent\nreference pathways" = "#5F7593"
  )
  linetypes <- c(
    "Degradation obs." = "solid",
    "Degradation forecast" = "dotted",
    "Improvement obs." = "solid",
    "Improvement forecast" = "longdash",
    "KM-GBF target-consistent\nreference pathways" = "dotted",
    "G20 target-consistent\nreference pathways" = "dotted"
  )
  widths <- c(
    "Degradation obs." = 2.25,
    "Degradation forecast" = 2.05,
    "Improvement obs." = 2.25,
    "Improvement forecast" = 2.05,
    "KM-GBF target-consistent\nreference pathways" = 1.20,
    "G20 target-consistent\nreference pathways" = 1.20
  )

  ymax <- max(plot_data$area, na.rm = TRUE)
  y_upper <- max(pretty(c(0, ymax * 1.05), n = 5))
  y_breaks <- pretty(c(0, y_upper), n = 5)
  y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= y_upper]
  y_accuracy <- min(diff(y_breaks))

  p <- ggplot(plot_data, aes(x = year, y = area, colour = series,
                             linetype = series, linewidth = series, group = series)) +
    geom_line(lineend = "round", na.rm = TRUE) +
    geom_point(
      data = target_points,
      aes(x = year, y = area, colour = series),
      inherit.aes = FALSE,
      size = 4.0,
      shape = 16,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = colours, breaks = series_levels, drop = FALSE) +
    scale_linetype_manual(values = linetypes, breaks = series_levels, drop = FALSE) +
    scale_linewidth_manual(values = widths, breaks = series_levels, drop = FALSE) +
    scale_x_continuous(
      breaks = c(2016, 2024, 2030, 2040, 2050),
      limits = c(2016, 2050),
      expand = expansion(mult = c(0.015, 0.055))
    ) +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = y_breaks,
      labels = label_number(accuracy = y_accuracy, trim = TRUE),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL,
      colour = NULL,
      linetype = NULL,
      linewidth = NULL
    ) +
    theme_bw(base_size = 20, base_family = "sans") +
    theme(
      panel.border = element_rect(colour = "#808080", fill = NA, linewidth = 0.75),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "#DDDDDD", linewidth = 0.55),
      axis.ticks = element_line(colour = "#222222", linewidth = 0.75),
      axis.ticks.length = grid::unit(2.5, "mm"),
      axis.text = element_text(colour = "#222222", size = 19),
      axis.title = element_text(colour = "#111111", face = "bold", size = 21),
      plot.title = element_text(face = "bold", size = 22, hjust = 0),
      legend.position = "none",
      legend.text = element_text(size = 17, lineheight = 1.05),
      legend.key.width = grid::unit(15, "mm"),
      legend.key.height = grid::unit(7, "mm"),
      legend.spacing.y = grid::unit(0.5, "mm"),
      plot.margin = margin(7, 12, 7, 7)
    )

  ragg::agg_png(
    file.path(output_dir, paste0("sdg1531_theilsen_", stem_i, "_2016_2050.png")),
    width = 7.00,
    height = 4.67,
    units = "in",
    res = 300,
    background = "white"
  )
  print(p)
  dev.off()
}

model_summary <- bind_rows(model_rows)
forecast_data <- bind_rows(forecast_rows)
target_data <- bind_rows(target_rows)
all_plot_data <- bind_rows(plot_rows)

write.csv(model_summary, file.path(output_dir, "theil_sen_model_summary.csv"), row.names = FALSE)
write.csv(forecast_data, file.path(output_dir, "theil_sen_forecasts_2025_2050.csv"), row.names = FALSE)
write.csv(target_data, file.path(output_dir, "target_consistent_reference_pathways.csv"), row.names = FALSE)
write.csv(all_plot_data, file.path(output_dir, "figure_source_data.csv"), row.names = FALSE)

writeLines(
  c(
    "Methods and QA notes",
    "Observed areas were aggregated by year; global values sum all geographic records, and country values sum records sharing an ISO3 code.",
    "Theil-Sen fits use the median of all pairwise slopes for 2016-2024 and the median residual intercept.",
    "Forecasts cover 2025-2050 and are constrained to non-negative areas.",
    "KM-GBF pathway: smooth logistic decline from the mapped 2022 degradation area through 70% of that area in 2030.",
    "G20 pathway: smooth logistic decline from the mapped 2024 degradation area through 50% of that area in 2040.",
    "The target statement defines endpoints but not the annual trajectory; logistic paths provide smooth target-consistent reference curves.",
    "Units: million km2. Output PNG dimensions: 7.00 x 4.67 inches at 300 dpi (2100 x 1401 pixels)."
  ),
  file.path(output_dir, "README_methods_and_QA.txt")
)

cat("Created output in:", output_dir, "\n")
print(list.files(output_dir))
