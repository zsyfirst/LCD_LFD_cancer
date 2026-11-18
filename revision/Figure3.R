library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)

fig_dir <- "lowcarb_cancer/revision/figure"
res <- read_tsv("lowcarb_cancer/revision/tableS7.res.txt", show_col_types = FALSE)
casepy <- read_csv("lowcarb_cancer/revision/table/casePY_pool_can_obese.csv", show_col_types = FALSE)


# ============================================================
# Labels
# ============================================================

score_labs <- c( 
      #"lcssvq"  = "Overall LCD",
      "alcssvq" = "Animal LCD",
      "ulcssvq" = "Unhealthy LCD",
      "vlcssvq" = "Vegetable LCD",
      "hlcssvq" = "Healthy LCD",
      # "lfssvq"  = "Overall LFD",
      "alfssvq" = "Animal LFD",
      "ulfssvq" = "Unhealthy LFD",
      "vlfssvq" = "Vegetable LFD",
      "hlfssvq" = "Healthy LFD"
)

trend_labs <- c(
      #"lcssvtrend"  = "Overall LCD",
      "alcssvtrend" = "Animal LCD",
      "ulcssvtrend" = "Unhealthy LCD",
      "vlcssvtrend" = "Vegetable LCD",
      "hlcssvtrend" = "Healthy LCD",
      # "lfssvtrend"  = "Overall LFD",
      "alfssvtrend" = "Animal LFD",
      "ulfssvtrend" = "Unhealthy LFD",
      "vlfssvtrend" = "Vegetable LFD",
      "hlfssvtrend" = "Healthy LFD"
)

score_order <- c(
      # "Overall LCD", 
	"Animal LCD", "Unhealthy LCD", "Vegetable LCD", "Healthy LCD",
      # "Overall LFD", 
	"Animal LFD", "Unhealthy LFD", "Vegetable LFD", "Healthy LFD"
)

primary_prev <- c(
  "alcssvq", "ulcssvq", "vlcssvq", "hlcssvq",
  "alfssvq", "ulfssvq", "vlfssvq", "hlfssvq"
)

primary_trend <- c(
  "alcssvtrend", "ulcssvtrend", "vlcssvtrend", "hlcssvtrend",
  "alfssvtrend", "ulfssvtrend", "vlfssvtrend", "hlfssvtrend"
)

# ============================================================
# Categorical HR data: Q2-Q5 vs Q1
# ============================================================

hr_cat <- res %>%
      filter(analysis == "avestop.adj2") %>%
      mutate(
            postv = as.numeric(postv),
            diet_score = recode(prev, !!!score_labs),
            HR_CI = sprintf("%.2f (%.2f, %.2f)", hazardratio, lcl, ucl)
      ) %>%
      select(prev, diet_score, postv, hazardratio, lcl, ucl, HR_CI)

# Add Q1 reference rows
hr_ref <- hr_cat %>%
      distinct(prev, diet_score) %>%
      mutate(
            postv = 1,
            hazardratio = 1,
            lcl = 1,
            ucl = 1,
            HR_CI = "1.00"
      )

hr_dat <- bind_rows(hr_ref, hr_cat)

# ============================================================
# P for trend
# ============================================================

trend_dat <- res %>%
      filter(analysis == "trend.adj2") %>%
      mutate(
            diet_score = recode(variable, !!!trend_labs),
            p_trend_num = suppressWarnings(as.numeric(probchisq)),
            p_trend = case_when(
                  probchisq == "<.0001" ~ "<0.001",
                  p_trend_num < 0.001 ~ "<0.001",
                  TRUE ~ sprintf("%.3f", p_trend_num)
            )
      ) %>%
      select(diet_score, p_trend)

# ============================================================
# Cases/person-years
# Assumes casePY has index and columns _1 to _5
# ============================================================

casepy_long <- casepy %>%
      pivot_longer(
            cols = matches("^_?[1-5]$"),
            names_to = "postv",
            values_to = "cases_py"
      ) %>%
      mutate(
            postv = as.numeric(str_remove(postv, "^_")),
            prev = if_else(str_ends(index, "q"), index, paste0(index, "q"))
      ) %>%
      select(prev, postv, cases_py)

casepy_long <- casepy_long %>%
      mutate(
            cases = str_extract(as.character(cases_py), "^[0-9,]+")
      )

# ============================================================
# Combine
# ============================================================

plot_dat <- hr_dat %>%
      left_join(casepy_long, by = c("prev", "postv")) %>%
      left_join(trend_dat, by = "diet_score") %>%
      mutate(
            group = if_else(
                  str_detect(diet_score, "LCD"),
                  "Low-carbohydrate diet",
                  "Low-fat diet"
            ),
            category = paste0("Q", postv),
            p_trend_show = if_else(postv == 1, p_trend, ""),
            significant = lcl > 1 | ucl < 1,
            direction = case_when(
                  significant & hazardratio > 1 ~ "Higher risk",
                  significant & hazardratio < 1 ~ "Lower risk",
                  TRUE ~ "Not significant"
            ),
            color_value = abs(log(hazardratio))
      ) |>
  filter(
    prev %in% primary_prev
  )

# ============================================================
# Manual row positions
# ============================================================
score_base_y <- c(
  "Animal LCD"     = 21,
  "Unhealthy LCD"  = 15.5,
  "Vegetable LCD"  = 10,
  "Healthy LCD"    = 4.5,
  
  "Animal LFD"     = -4.5,
  "Unhealthy LFD"  = -10,
  "Vegetable LFD"  = -15.5,
  "Healthy LFD"    = -21
)

plot_dat <- plot_dat %>%
  mutate(
    base_y = score_base_y[diet_score],
    row_id = base_y - (postv - 1),
    score_label = if_else(postv == 1, diet_score, ""),
    category_label = paste0("  ", category)
  )

header_dat <- tibble(
  group = c("Low-carbohydrate diet", "Low-fat diet"),
  row_id = c(23, -2.5)
)

y_min <- -25
y_max <- 25.5
header_y <- 25

# ============================================================
# Left columns
# ============================================================

p_left <- ggplot() +
      geom_text(
            data = header_dat,
            aes(x = 0, y = row_id, label = group),
            hjust = 0, fontface = "bold", size = 4.2
      ) +
      geom_text(
            data = plot_dat %>% filter(postv == 1),
            aes(x = 0, y = row_id, label = diet_score),
            hjust = 0, size = 4
      ) +
      geom_text(
            data = plot_dat,
            aes(x = 0.7, y = row_id, label = category_label),
            hjust = 0, size = 4.2
      ) +
      annotate(
            "text",
            x = 0, y = header_y,
            label = "Diet score",
            hjust = 0, fontface = "bold", size = 4.4
      ) +
      annotate(
            "text",
            x = 0.7, y = header_y,
            label = "Category",
            hjust = 0, fontface = "bold", size = 4.3
      ) +
      scale_y_continuous(limits = c(y_min, y_max)) +
      coord_cartesian(xlim = c(0, 1.1), clip = "off") +
      theme_void() +
      theme(plot.margin = margin(5, 0, 5, 5))

# ============================================================
# Cases/person-years column
# ============================================================

p_cases <- ggplot() +
      geom_text(
            data = plot_dat,
            aes(x = 0, y = row_id, label = cases),
            hjust = 0, size = 4
      ) +
      annotate(
            "text",
            x = 0, y = header_y,
            # label = "Cases/person-years",
            label = "No. of cases",
            hjust = 0, fontface = "bold", size = 4.2
      ) +
      scale_y_continuous(limits = c(y_min, y_max)) +
      coord_cartesian(xlim = c(0, 1), clip = "off") +
      theme_void() +
      theme(plot.margin = margin(5, -3, 5, 0))

# ============================================================
# Middle forest plot
# ============================================================

p_mid <- ggplot(plot_dat, aes(y = row_id)) +
      geom_vline(
            xintercept = 1,
            linetype = "dashed",
            color = "grey55",
            linewidth = 0.45
      ) +
      geom_errorbarh(
            aes(xmin = lcl, xmax = ucl),
            height = 0,
            linewidth = 0.55,
            color = "grey35"
      ) +
      geom_point(
            aes(
                  x = hazardratio,
                  shape = direction,
                  fill = color_value,
                  color = direction
            ),
            size = 3.1,
            stroke = 0.85
      ) +
      scale_shape_manual(
            values = c(
                  "Higher risk" = 21,
                  "Lower risk" = 23,
                  "Not significant" = 22
            )
      ) +
      scale_color_manual(
            values = c(
                  "Higher risk" = "#08306B",
                  "Lower risk" = "#006D2C",
                  "Not significant" = "grey35"
            )
      ) +
      scale_fill_gradientn(
            colours = c("#deebf7", "#9ecae1", "#3182bd", "#08519c"),
            guide = "none"
      ) +
      scale_x_continuous(
            limits = c(0.76, 1.22),
            breaks = seq(0.8, 1.2, 0.1)
      ) +
      scale_y_continuous(limits = c(y_min, y_max)) +
      labs(
            x = "Hazard ratio (95% CI)",
            y = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
            legend.position = "none",
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 11),
            axis.title.x = element_text(size = 12),
            plot.margin = margin(5, -2, 5, -2)
      )

# ============================================================
# Right text columns
# ============================================================

p_right <- ggplot() +
      geom_text(
            data = plot_dat,
            aes(x = 0, y = row_id, label = HR_CI),
            hjust = 0, size = 4
      ) +
      geom_text(
            data = plot_dat,
            aes(x = 1.15, y = row_id, label = p_trend_show),
            hjust = 0.5, size = 4
      ) +
      annotate(
            "text",
            x = 0, y = header_y,
            label = "HR (95% CI)",
            hjust = 0, fontface = "bold", size = 4.2
      ) +
      annotate(
            "text",
            x = 1.15, y = header_y,
            label = "P for trend",
            hjust = 0.5, fontface = "bold", size = 4.2
      ) +
      scale_y_continuous(limits = c(y_min, y_max)) +
      coord_cartesian(xlim = c(0, 1.45), clip = "off") +
      theme_void() +
      theme(plot.margin = margin(5, 5, 5, 0))

# ============================================================
# Combine
# ============================================================

p_final <- p_left + p_cases + p_mid + p_right +
      plot_layout(widths = c(1.35, 0.6, 1.25, 1.15))

# p_final

# ============================================================
# Save
# ============================================================

ggsave(
      file.path(fig_dir, "Table23_Q1_to_Q5_forest.pdf"),
      p_final,
      width = 10,
      height = 12,
      device = cairo_pdf
)

ggsave(
      file.path(fig_dir, "Table23_Q1_to_Q5_forest.png"),
      p_final,
      width = 10,
      height = 12,
      dpi = 600
)

