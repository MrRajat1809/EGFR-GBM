# Locate shared paths for Rscript or source().
.egfr_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.egfr_script)) {
  .egfr_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(.egfr_arg)) stop("Run with Rscript or source this file.")
  .egfr_script <- sub("^--file=", "", .egfr_arg[[1]])
}
.egfr_src <- dirname(normalizePath(.egfr_script, mustWork = TRUE))
while (!file.exists(file.path(.egfr_src, "common", "02_project_paths.R"))) {
  .egfr_parent <- dirname(.egfr_src)
  if (.egfr_parent == .egfr_src) stop("Cannot locate src/common/02_project_paths.R")
  .egfr_src <- .egfr_parent
}
source(file.path(.egfr_src, "common", "02_project_paths.R"))
use_working_directory("outputs/figures")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 1. Parse the raw SOPMA.txt file
lines <- readLines("SOPMA.txt")
data_lines <- lines[grepl("^[CHET]\\s+\\d+", lines)]

# Create empty dataframe
df <- data.frame(
  Position = integer(length(data_lines)),
  Consensus = character(length(data_lines)),
  Helix = numeric(length(data_lines)),
  Sheet = numeric(length(data_lines)),
  Turn = numeric(length(data_lines)),
  Coil = numeric(length(data_lines)),
  AA = character(length(data_lines)),
  stringsAsFactors = FALSE
)

# Loop and extract columns
for (i in seq_along(data_lines)) {
  parts <- strsplit(trimws(data_lines[i]), "\\s+")[[1]]
  df$Position[i] <- i
  df$Consensus[i] <- parts[1]
  df$Helix[i] <- as.numeric(parts[2])
  df$Sheet[i] <- as.numeric(parts[3])
  df$Turn[i] <- as.numeric(parts[4])
  df$Coil[i] <- as.numeric(parts[5])
  df$AA[i] <- parts[6]
}

# 2. Format Data
df <- df %>%
  mutate(
    Structure = case_when(
      Consensus == "H" ~ "Helix",
      Consensus == "E" ~ "Sheet",
      Consensus == "T" ~ "Turn",
      Consensus == "C" ~ "Coil"
    ),
    Structure = factor(Structure, levels = c("Helix", "Sheet", "Turn", "Coil")),
    ymin = ifelse(Structure == "Coil", 0.45, 0.1),
    ymax = ifelse(Structure == "Coil", 0.55, 0.9)
  )

plot_colors <- c(
  "Helix" = "#4DBBD5", 
  "Sheet" = "#9C27B0", 
  "Turn"  = "#00A087", 
  "Coil"  = "#F39C12"  
)

# 3. Create Panel 1: Top Structure Map
p_top <- ggplot(df) +
  geom_hline(yintercept = 0.5, color = plot_colors["Coil"], linewidth = 0.3) +
  geom_rect(aes(xmin = Position - 0.425, xmax = Position + 0.425, 
                ymin = ymin, ymax = ymax, fill = Structure)) +
  scale_fill_manual(values = plot_colors) +
  theme_bw() +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1.2, fill = NA), # ADDED: Bolder border
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.text.x = element_text(margin = margin(t = 8), size = 12, color = "black"),
    legend.position = "none" 
  )

# 4. Prepare Data for Panel 2
df_long <- df %>%
  select(Position, Helix, Sheet, Turn, Coil) %>%
  pivot_longer(cols = -Position, names_to = "Structure", values_to = "Score") %>%
  mutate(Structure = factor(Structure, levels = c("Helix", "Sheet", "Turn", "Coil")))

# 5. Create Panel 2: The Line Graph
p_bottom <- ggplot(df_long, aes(x = Position, y = Score, color = Structure)) +
  geom_line(linewidth = 0.5, alpha = 0.85) + 
  scale_color_manual(values = plot_colors) +
  theme_bw() +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.05, 0.05)) + 
  theme(
    panel.border = element_rect(color = "black", linewidth = 1.2, fill = NA), # ADDED: Bolder border
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.text.x = element_text(margin = margin(t = 8), size = 12, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black", margin = margin(r = 15, l = 5)),
    legend.key = element_rect(fill = NA, color = NA),
    legend.key.width = unit(1.5, "cm") # Slightly reduced so thicker lines don't look overly long
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) # ADDED: Thicker legend lines

# 6. Combine the Plots
combined_plot <- p_top / p_bottom + 
  plot_layout(heights = c(1, 2.5)) # Slightly adjusted ratio to match the sample image

# 7. Export as High-Resolution Image
# REDUCED: width changed from 12 to 9, height from 6 to 5 for a tighter aspect ratio
ggsave("SOPMA_Customized_Format.png", plot = combined_plot, width = 9, height = 5, dpi = 600)

print("Success! The visually customized version has been saved.")