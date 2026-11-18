library(data.table)
library(tidyverse)
library(ggcorrplot)
library(corrplot)
library(patchwork)

fig_dir <- "lowcarb_cancer/revision/food_items"
pool <- fread("lowcarb_cancer/data/pool_allcancer.csv")

diet_vars <- c(
  "ulcssvcon", "alcssvcon",  "vlcssvcon", "hlcssvcon","lcssvcon",
  "ulfssvcon",  "alfssvcon", "vlfssvcon", "hlfssvcon", "lfssvcon"
)

diet_labs <- c(
  "Unhealthy LCD", "Animal LCD",  "Vegetable LCD", "Healthy LCD", "Overall LCD", 
  "Unhealthy LFD", "Animal LFD",  "Vegetable LFD", "Healthy LFD",  "Overall LFD"
)

diet_order <- c(
      "Unhealthy LCD", "Animal LCD",  "Vegetable LCD", "Healthy LCD", "Overall LCD",
      "Unhealthy LFD", "Animal LFD",  "Vegetable LFD", "Healthy LFD", "Overall LFD"
)

################ Macronutrients correlation heatmap ################

nutrient_vars <- c("carb", "sumref_carb", "hcarb", "prot", "aprot", "vprot", "tfat", "afat", "vfat")
nutrient_labs <- c("Total Carb", "Unhealthy Carb", "Healthy Carb",
                   "Total Protein", "Animal Protein", "Plant Protein",
                   "Total Fat", "Animal Fat", "Plant Fat")

cordata_nutrient <- pool |>
      select(all_of(c(nutrient_vars, diet_vars))) |>
      drop_na()
names(cordata_nutrient) <- c(nutrient_labs, diet_labs)
cal_cor_nutrient <- round(cor(cordata_nutrient, method = "spearman"), 3)
cal_cor_long <- as.data.frame(cal_cor_nutrient[diet_labs, nutrient_labs]) |>
      rownames_to_column("Diet") |>
      pivot_longer(-Diet, names_to = "Nutrient", values_to = "Correlation") |>
      mutate(
            Diet = factor(Diet, levels = diet_order),
            Nutrient = factor(Nutrient, levels = rev(nutrient_labs))
      )

p_nutrient <- ggplot(cal_cor_long, aes(x = Diet, y = Nutrient, fill = Correlation)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f", Correlation)), size = 2.8) +
      scale_fill_gradient2(
            low = "#6194BF", mid = "white", high = "#E3732A",
            midpoint = 0, limits = c(-1, 1)
      ) +
      theme_minimal() +
      theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
            axis.text.y = element_text(size = 9),
            axis.title = element_blank(),
            panel.grid = element_blank()
      ) +
      labs(fill = "Spearman r")

ggsave(file.path(fig_dir, "cor_cum_nutrients.pdf"),p_nutrient, width = 8, height = 4)

################### Food item correlation heatmap #####################
food_vars <- c("whgrn_svcon", "rgrain_svcon", "pota_svcon", "veg_noleg_svcon",
               "frui_svcon", "leg_svcon", "fish_svcon", "rpmeats_svcon",
               "poul_svcon", "ssb_svcon", "dairy_svcon", "swtsnochoc_svcon")
food_labs <- c("Whole grain", "Refined grain", "Potato", "Vegetable",
               "Fruit", "Legume", "Seafood", "Red meat",
               "Poultry", "SSB", "Dairy", "Sweets")
cordata_food <- pool |>
      select(all_of(c(food_vars, diet_vars))) |>
      drop_na()
names(cordata_food) <- c(food_labs, diet_labs)
cal_cor_food <- round(cor(cordata_food, method = "spearman"), 3)
cal_cor_long_food <- as.data.frame(cal_cor_food[diet_labs, food_labs]) |>
      rownames_to_column("Diet") |>
      pivot_longer(-Diet, names_to = "Food Items", values_to = "Correlation") |>
      mutate(
            Diet = factor(Diet, levels = diet_order),
            `Food Items` = factor(`Food Items`, levels = rev(food_labs)) )

p_food <- ggplot(cal_cor_long_food, aes(x = Diet, y = `Food Items`, fill = Correlation)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f", Correlation)), size = 2.8) +
      scale_fill_gradient2(
            low = "#6194BF", mid = "white", high = "#E3732A",
            midpoint = 0, limits = c(-1, 1)
      ) +
      theme_minimal() +
      theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
            axis.text.y = element_text(size = 9),
            axis.title = element_blank(),
            panel.grid = element_blank()
      ) +
      labs(fill = "Spearman r")

ggsave(file.path(fig_dir, "cor_cum_food_item.pdf"),
      p_food, width = 8, height = 4)


#######################################################################
############# Combined Figure: LCD/LFD and nutritents & food items #############
p_nutrient_top <- p_nutrient +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(tag = "Macronutrients")
p_food_bottom <- p_food +labs(tag = "Food Items")

p_combined <- p_nutrient_top / p_food_bottom +
      plot_layout(guides = "collect", heights = c(9, 12)) &
      theme(legend.position = "right")

ggsave(file.path(fig_dir, "cor_cum_combined.pdf"),
      p_combined, width = 10, height = 10)


##########################################################################################
######################## LCD/LFD and AHEI correlation plot ###############################
ahei_vars <- c(diet_vars, "aheimsvcon")
ahei_labs <- c(diet_labs, "AHEI-2010")

cordata_AHEI <- pool |>
      select(all_of(ahei_vars)) |>
      drop_na()

names(cordata_AHEI) <- ahei_labs

cal_cor_AHEI <- round(cor(cordata_AHEI, method = "spearman"), 3)
cal_cor_p_AHEI <- round(cor_pmat(cordata_AHEI, method = "spearman"), 3)

pdf(file.path(fig_dir, "cor_cum_AHEI.pdf"), width = 9, height = 9)

corrplot(
      corr = cal_cor_AHEI,
      p.mat = cal_cor_p_AHEI,
      method = "color",
      type = "upper",
      tl.pos = "lt",
      tl.cex = 1.3,
      tl.col = "black",
      tl.srt = 45,
      tl.offset = 0.5,
      insig = "label_sig",
      sig.level = c(0.01, 0.05),
      pch.cex = 1.5,
      pch.col = "black"
)

darkblue80 <- rgb(0, 0, 139, maxColorValue = 255)

corrplot(
      corr = cal_cor_AHEI,
      type = "lower",
      add = TRUE,
      method = "number",
      tl.pos = "n",
      cl.pos = "n",
      diag = FALSE,
      number.digits = 2,
      number.cex = 1.0,
      col = colorRampPalette(c(darkblue80, darkblue80))(200)
)

dev.off()
