# DE for CRC_GSE166555
library(ggplot2)
library(ggrepel)
library(patchwork)

data_1 <- read.csv("CRC-p.csv")
data_2 <- read.csv("CRC-c.csv")

# Extracting Data Rows
a5 <- data_1[1, 2:17]
b5 <- data_2[1, 2:13]
c5 <- data_1[2, 2:17]
d5 <- data_2[2, 2:13]
S7a <- data_1[3, 2:17]
S7b <- data_1[4, 2:17]
S7c <- data_1[5, 2:17]
S7d <- data_1[6, 2:17]

# Unified style parameters
base_size <- 10
axis_text_size <- 7
axis_title_size <- 7
legend_text_size <- 7
bar_text_size <- 2
axis_line_size <- 0.1
margin_top <- 3
margin_right <- 3
y_expand <- 0.3

# Color Scheme
group_colors <- c("IntOGen" = "#C6B1B2",
                  "DriverDBv4" = "#BEE5BF", 
                  "Our Method" = "#FFD1BA",
                  "GSE146771" = "#C6B1B2", 
                  "EMTAB8107" = "#BEE5BF", 
                  "GSE166555" = "#FFD1BA")

# Fix issue with period in column names
fix_column_names <- function(names) {
  names <- gsub("\\.", " ", names)
  names <- gsub("Our Method", "Our Method", names)
  return(names)
}

# Create theme settings for legend spacing adjustments
legend_theme <- theme(
  legend.margin = margin(t = -5, unit = "pt"),  # Reduce legend top spacing
  legend.box.margin = margin(b = -10)           # Reduce legend bottom spacing
)

# 1. Prepare and draw a5 graphics
plot_a5 <- data.frame(
  Method = fix_column_names(names(a5)),  # Fix column names
  Value = as.numeric(a5[1, ]),
  Group = factor(
    c(rep("IntOGen", 7), 
      rep("DriverDBv4", 8),
      "Our Method"),
    levels = c("IntOGen", "DriverDBv4", "Our Method")
  )
)
plot_a5$Method <- factor(plot_a5$Method, levels = fix_column_names(names(a5)))  # Repair factor level

plot_5a <- ggplot(plot_a5, aes(x = Method, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Value), vjust = -0.5, size = bar_text_size, color = "black") +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(x = "Methods", y = "Number of predicted genes") +
  theme_classic(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.title.x = element_text(margin = margin(t = margin_top), size = axis_title_size),
    axis.title.y = element_text(margin = margin(r = margin_right), size = axis_title_size),
    legend.position = "top",
    legend.text = element_text(size = legend_text_size),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_blank(),
    axis.line = element_line(size = axis_line_size, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid.major.y = element_line(size = 0.25, colour = "grey90"),
    text = element_text(family = "sans", face = "plain")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, y_expand))) +
  legend_theme  # 添加图例间距调整

# 2. b5
plot_b5 <- data.frame(
  Position = 1:ncol(b5),
  Label = fix_column_names(gsub("\\.\\d+$", "", names(b5))),  # Fix column names
  Value = as.numeric(b5[1, ]),
  Group = factor(
    c(rep("GSE146771", 3), rep("EMTAB8107", 3), rep("GSE166555", 6)),
    levels = c("GSE146771", "EMTAB8107", "GSE166555")
  )
)

plot_5b <- ggplot(plot_b5, aes(x = Position, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Value), vjust = -0.5, size = bar_text_size, color = "black") +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(x = "Subtypes", y = "Number of predicted genes") +
  theme_classic(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.title.x = element_text(margin = margin(t = margin_top), size = axis_title_size),
    axis.title.y = element_text(margin = margin(r = margin_right), size = axis_title_size),
    legend.position = "top",
    legend.text = element_text(size = legend_text_size),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_blank(),
    axis.line = element_line(size = axis_line_size, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid.major.y = element_line(size = 0.25, colour = "grey90"),
    text = element_text(family = "sans", face = "plain")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, y_expand))) +
  scale_x_continuous(breaks = plot_b5$Position, labels = plot_b5$Label) +
  legend_theme  # 添加图例间距调整

# 3. c5
plot_c5 <- data.frame(
  Method = fix_column_names(names(c5)),  # 
  Value = as.numeric(c5[1, ]),
  Group = factor(
    c(rep("IntOGen", 7), rep("DriverDBv4", 8), "Our Method"),
    levels = c("IntOGen", "DriverDBv4", "Our Method")
  )
)
plot_c5$Method <- factor(plot_c5$Method, levels = fix_column_names(names(c5)))  # 

plot_5c <- ggplot(plot_c5, aes(x = Method, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Value), vjust = -0.5, size = bar_text_size, color = "black") +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(x = "Methods", y = "Number of predicted genes") +
  theme_classic(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.title.x = element_text(margin = margin(t = margin_top), size = axis_title_size),
    axis.title.y = element_text(margin = margin(r = margin_right), size = axis_title_size),
    legend.position = "top",
    legend.text = element_text(size = legend_text_size),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_blank(),
    axis.line = element_line(size = axis_line_size, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid.major.y = element_line(size = 0.25, colour = "grey90"),
    text = element_text(family = "sans", face = "plain")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, y_expand))) +
  legend_theme  # Add legend spacing adjustment

# 4. d5
plot_d5 <- data.frame(
  Position = 1:ncol(d5),
  Label = fix_column_names(gsub("\\.\\d+$", "", names(d5))),  # Fix column names
  Value = as.numeric(d5[1, ]),
  Group = factor(
    c(rep("GSE146771",3), rep("EMTAB8107", 3), rep("GSE166555", 6)),
    levels = c("GSE146771", "EMTAB8107", "GSE166555")
  )
)

plot_5d <- ggplot(plot_d5, aes(x = Position, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Value), vjust = -0.5, size = bar_text_size, color = "black") +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(x = "Subtypes", y = "Number of predicted genes") +
  theme_classic(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.title.x = element_text(margin = margin(t = margin_top), size = axis_title_size),
    axis.title.y = element_text(margin = margin(r = margin_right), size = axis_title_size),
    legend.position = "top",
    legend.text = element_text(size = legend_text_size),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_blank(),
    axis.line = element_line(size = axis_line_size, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid.major.y = element_line(size = 0.25, colour = "grey90"),
    text = element_text(family = "sans", face = "plain")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, y_expand))) +
  scale_x_continuous(breaks = plot_d5$Position, labels = plot_d5$Label) +
  legend_theme  # 添加图例间距调整

# 5. S7a
plot_S7a <- data.frame(
  Method = fix_column_names(names(S7a)),  # 
  Value = as.numeric(S7a[1, ]),
  Group = factor(
    c(rep("IntOGen", 7), rep("DriverDBv4", 8), "Our Method"),
    levels = c("IntOGen", "DriverDBv4", "Our Method")
  )
)
plot_S7a$Method <- factor(plot_S7a$Method, levels = fix_column_names(names(S7a)))  # 

plot_S7a <- ggplot(plot_S7a, aes(x = Method, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Value), vjust = -0.5, size = bar_text_size, color = "black") +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(x = "Methods", y = "Number of predicted genes") +
  theme_classic(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.title.x = element_text(margin = margin(t = margin_top), size = axis_title_size),
    axis.title.y = element_text(margin = margin(r = margin_right), size = axis_title_size),
    legend.position = "top",
    legend.text = element_text(size = legend_text_size),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_blank(),
    axis.line = element_line(size = axis_line_size, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid.major.y = element_line(size = 0.25, colour = "grey90"),
    text = element_text(family = "sans", face = "plain")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, y_expand))) +
  legend_theme  # 

# 6. S7b图形
plot_S7b <- data.frame(
  Method = fix_column_names(names(S7b)),  # 
  Value = as.numeric(S7b[1, ]),
  Group = factor(
    c(rep("IntOGen", 7), rep("DriverDBv4", 8), "Our Method"),
    levels = c("IntOGen", "DriverDBv4", "Our Method")
  )
)
plot_S7b$Method <- factor(plot_S7b$Method, levels = fix_column_names(names(S7b)))  # 

plot_S7b <- ggplot(plot_S7b, aes(x = Method, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Value), vjust = -0.5, size = bar_text_size, color = "black") +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(x = "Methods", y = "Number of predicted genes") +
  theme_classic(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.title.x = element_text(margin = margin(t = margin_top), size = axis_title_size),
    axis.title.y = element_text(margin = margin(r = margin_right), size = axis_title_size),
    legend.position = "top",
    legend.text = element_text(size = legend_text_size),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_blank(),
    axis.line = element_line(size = axis_line_size, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid.major.y = element_line(size = 0.25, colour = "grey90"),
    text = element_text(family = "sans", face = "plain")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, y_expand))) +
  legend_theme  # 

# 7. S7c
plot_S7c <- data.frame(
  Method = fix_column_names(names(S7c)),  # 
  Value = as.numeric(S7c[1, ]),
  Group = factor(
    c(rep("IntOGen", 7), rep("DriverDBv4", 8), "Our Method"),
    levels = c("IntOGen", "DriverDBv4", "Our Method")
  )
)
plot_S7c$Method <- factor(plot_S7c$Method, levels = fix_column_names(names(S7c)))  # 

plot_S7c <- ggplot(plot_S7c, aes(x = Method, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Value), vjust = -0.5, size = bar_text_size, color = "black") +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(x = "Methods", y = "Number of predicted genes") +
  theme_classic(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.title.x = element_text(margin = margin(t = margin_top), size = axis_title_size),
    axis.title.y = element_text(margin = margin(r = margin_right), size = axis_title_size),
    legend.position = "top",
    legend.text = element_text(size = legend_text_size),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_blank(),
    axis.line = element_line(size = axis_line_size, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid.major.y = element_line(size = 0.25, colour = "grey90"),
    text = element_text(family = "sans", face = "plain")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, y_expand))) +
  legend_theme  # 添加图例间距调整

# 8. S7d
plot_S7d <- data.frame(
  Method = fix_column_names(names(S7d)),  # 
  Value = as.numeric(S7d[1, ]),
  Group = factor(
    c(rep("IntOGen", 7), rep("DriverDBv4", 8), "Our Method"),
    levels = c("IntOGen", "DriverDBv4", "Our Method")
  )
)
plot_S7d$Method <- factor(plot_S7d$Method, levels = fix_column_names(names(S7d)))  #

plot_S7d <- ggplot(plot_S7d, aes(x = Method, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Value), vjust = -0.5, size = bar_text_size, color = "black") +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(x = "Methods", y = "Number of predicted genes") +
  theme_classic(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.title.x = element_text(margin = margin(t = margin_top), size = axis_title_size),
    axis.title.y = element_text(margin = margin(r = margin_right), size = axis_title_size),
    legend.position = "top",
    legend.text = element_text(size = legend_text_size),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_blank(),
    axis.line = element_line(size = axis_line_size, colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid.major.y = element_line(size = 0.25, colour = "grey90"),
    text = element_text(family = "sans", face = "plain")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, y_expand))) +
  legend_theme  # 添加图例间距调整

# print(plot_5a)
# print(plot_5b)
# print(plot_5c)
# print(plot_5d)
# print(plot_S7a)
# print(plot_S7b)
# print(plot_S7c)
# print(plot_S7d)


# # Create a graphic combination of 4 rows and 2 columns, and move 5b and 5d to the last row
combined_4x2 <- (
  # First row: 5a and 5c (originally from the left side of the first row and the left side of the second row)
  (plot_5a + theme(plot.margin = margin(0, 20, 0, 20)) | 
     plot_5c + theme(plot.margin = margin(0, 0, 0, 20))) /
    
    # Second row: S7a and S7b (formerly third row)
    (plot_S7a + theme(plot.margin = margin(0, 20, 0, 20)) | 
       plot_S7b + theme(plot.margin = margin(0, 0, 0, 20))) /
    
    # Third row: S7c and S7d (formerly fourth row)
    (plot_S7c + theme(plot.margin = margin(0, 20, 0, 20)) | 
       plot_S7d + theme(plot.margin = margin(0, 20, 0, 20))) /
    
    # Fourth row (new): 5b and 5d (originally on the right side of the first row and the right side of the second row)
    (plot_5b + theme(plot.margin = margin(0, 20, 0, 20)) | 
       plot_5d + theme(plot.margin = margin(0, 20, 0, 20)))
) +
  patchwork::plot_annotation(
    tag_levels = 'a',
    theme = ggplot2::theme(
      plot.tag.position = "topleft",
      plot.tag = ggplot2::element_text(size = 12)
    )
  ) +
  patchwork::plot_layout(widths = c(1, 1))  # Make sure both columns are equal width

# 2. Force the theme to be applied to all sub-images
combined_4x2 <- combined_4x2 & 
  theme(plot.tag = element_text(size = 7.2))

# print
print(combined_4x2)

# save
ggsave(filename = "methods_comparison.pdf",  
       plot = combined_4x2,             
       width = 12,  #                
       height = 10, #                   
       dpi = 300              
)








