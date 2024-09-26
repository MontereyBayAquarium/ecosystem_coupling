#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

#required packages
librarian::shelf(tidyverse, here)

#set directories 
localdir <- here::here("output")
figdir <- here::here("figures")

# Get rocky intertidal position data
load(file.path(localdir, "processed/rocky_intertidal/position_data.rdata"))

# Calculate the range of points for each site across all years and identify points unique to 'after'
#We use !any(ssw_period == "before") to check if there are no rows with 'ssw_period'
#equal to "before" within each group of 'intertidal_sitename,' 'transect,' and 'location.' 
#If there are no such rows, it means the points are unique to the 'after' period, and we label them as "Only After."

heatmap_data <- mus_pos_build1 %>%
  group_by(intertidal_sitename, transect, location) %>%
  summarize(
    Period = factor(!any(ssw_period == "Before"), levels = c(TRUE, FALSE), labels = c("Post-SSW", "Pre-SSW"))
  ) %>%
  arrange(intertidal_sitename, transect) %>%
  group_by(intertidal_sitename) %>%
  mutate(transect_seq = dense_rank(transect),
         intertidal_sitename = factor(intertidal_sitename,levels = c("Hopkins","Point Pinos","Stillwater","Point Lobos")),
         Period = factor(Period, levels = c("Pre-SSW","Post-SSW")))


################################################################################
#Step 1 - convert frequency data to long format and calculate mean

#expand frequency
DatU <- vcdExtra::expand.dft(mus_size_build1, freq="total")

# Calculate the mean size by ssw_period
mean_size_by_period <- DatU %>%
  group_by(ssw_period) %>%
  summarize(mean_size = mean(as.numeric(size_bin)))

################################################################################
#Step 2 - test for significant change in mussel depth distribtion


depth_values <- heatmap_data %>%
  group_by(Period) %>%
  summarise(mean_depth = mean(location, na.rm = TRUE),
            depth_sd = sd(location, na.rm = TRUE))

print(depth_values)

# Filter data for Pre-SSW and Post-SSW periods
pre_ssw <- heatmap_data %>% filter(Period == "Pre-SSW") %>% pull(location)
post_ssw <- heatmap_data %>% filter(Period == "Post-SSW") %>% pull(location)

# t-test
t_test_result <- t.test(pre_ssw, post_ssw)

print(t_test_result)

################################################################################
#Step 3 - test for significant change in mussel size

size_values <- DatU %>%
  group_by(ssw_period) %>%
  summarise(mean_size = mean(size_bin, na.rm = TRUE),
            size_sd = sd(size_bin, na.rm = TRUE))

print(size_values)

# Filter data for Pre-SSW and Post-SSW periods
pre_ssw <- DatU %>% filter(ssw_period == "Before") %>% pull(size_bin)
post_ssw <- DatU %>% filter(ssw_period == "After") %>% pull(size_bin)

# t-test
t_test_result <- t.test(pre_ssw, post_ssw)

print(t_test_result)


################################################################################
#Step 4 - test for significant change in mussel cover


cov_values <- mus_cov_period %>%
  group_by(ssw_period) %>%
  summarise(mean_cov = mean(percent_cover, na.rm = TRUE),
            cov_sd = sd(percent_cover, na.rm = TRUE))

print(cov_values)

# Filter data for Pre-SSW and Post-SSW periods
pre_ssw <- mus_cov_period %>% filter(ssw_period == "Before") %>% pull(percent_cover)
post_ssw <- mus_cov_period %>% filter(ssw_period == "After") %>% pull(percent_cover)

# t-test
t_test_result <- t.test(pre_ssw, post_ssw)

print(t_test_result)



################################################################################
#plot heatmap by site


tile_theme <-  theme(axis.text=element_text(size=7, color = "black"),
                     axis.title=element_text(size=8,color = "black"),
                     plot.tag=element_text(size=8,color = "black"),
                     plot.title=element_text(size=8,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=5,color = "black"),
                     legend.title=element_text(size=6,color = "black"),
                     #legend.key.height = unit(0.1, "cm"),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=7, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())


# Create a heatmap with equal tile size and color points observed only in 'after' differently
p <- ggplot(heatmap_data, aes(x = transect_seq, y = location, fill = Period)) +
  geom_tile(width=1, height=1) +  # Set width and height to 1 for equal tile size
  scale_fill_manual(values = c("Pre-SSW" = "navyblue", "Post-SSW" = "indianred")) +
  labs(title = "",
       x = "Transect number",
       y = "Distance from high to \nlow intertidal (meters from baseline)",
       fill = NULL,
       tag = "A") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_reverse()+  # Flip the y-axis
  #ggforce::facet_col(vars(intertidal_sitename), scales = "free") +
  facet_grid(.~intertidal_sitename) +
  coord_fixed()+
  theme_bw()+tile_theme  + theme(legend.position = "top")+
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) #reduce right margin
p


# Create arrow scheme
arrow_data <- data.frame(x = 0.5, y = 0.85)

scheme1 <- ggplot(arrow_data, aes(x = x, y = y)) +
  geom_segment(aes(xend = x, y = 0.87, yend = y - 0.15),
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               size = 1) +
  labs(x = NULL, y = NULL) +  # Remove axis labels
  theme_void() +  # Remove axis lines and ticks
  annotate("text", x = 0.5, y = 0.92, label = "High \nintertidal", hjust = 0.5, vjust=5.5, size = 2) +
  annotate("text", x = 0.5, y = 0.68, label = "Low \nintertidal", hjust = 0.5, vjust=-0.2, size = 2)+
  theme(plot.margin = margin(0, 1, 0, 0, "cm")) #reduce left margin


# Merge
layout_matrix <- matrix(c(1,2), ncol=2)
g1_full <- gridExtra::grid.arrange(p, scheme1,
                                   layout_matrix=layout_matrix,
                                  widths=c(0.9, 0.1))
g1_full



#ggsave(g1_full, filename = file.path(figdir, "FigX_mussel_expansion.png"), 
 #      width = 7, height = 4, units = "in", dpi = 600)




# Theme
my_theme <-  theme(axis.text=element_text(size=8,color = "black"),
                   axis.title=element_text(size=9,color = "black"),
                   plot.tag=element_text(size=8,color = "black"),
                   plot.title=element_text(size=8,color = "black", face = "bold"),
                   # Gridlines
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key.size = unit(0.3, "cm"), 
                   #legend.key = element_rect(fill = "white"), # Set it to transparent
                   legend.spacing.y = unit(0.1, "cm"),  
                   legend.text=element_text(size=7,color = "black"),
                   legend.title=element_blank(),
                   #legend.key.height = unit(0.1, "cm"),
                   #legend.background = element_rect(fill=alpha('blue', 0)),
                   #facets
                   strip.text = element_text(size=7, face = "bold",color = "black", hjust=0),
                   strip.background = element_blank())

# Create a KDE plot for the mean size frequency distribution and overlay them
g2 <- ggplot(DatU, aes(x = as.numeric(size_bin), fill = ssw_period, color = ssw_period)) +
  geom_density(alpha = 0.8, adjust = 1.5) +
  labs(x = "Size \n(mm)", y = "Frequency", title = "Size frequency", tag = "C") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 120, by = 20)) +
  scale_fill_manual(values=c("indianred","navyblue"))+
  scale_color_manual(values=c("indianred","navyblue"))+
  #geom_vline(data = mean_size_by_period, aes(xintercept = mean_size, color = ssw_period), 
  #          linetype = "dotted", size = 1) +
  theme_bw() + my_theme + theme(axis.title.y = element_blank())
g2



# Create a KDE plot for the mean size frequency distribution and overlay them
g3 <- ggplot(mus_cov_period, aes(x = as.numeric(percent_cover), fill = ssw_period, color = ssw_period)) +
  geom_density(alpha = 0.8, adjust = 1.5) +
  labs(x = "Percent \ncover", y = "Frequency", title = "Percent cover", tag = "D") +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by = 10)) +
  scale_fill_manual(values=c("indianred","navyblue"))+
  scale_color_manual(values=c("indianred","navyblue"))+
  #geom_vline(data = mean_size_by_period, aes(xintercept = mean_size, color = ssw_period), 
  #          linetype = "dotted", size = 1) +
  theme_bw() + my_theme + theme(axis.title.y = element_blank())
g3



# Create a KDE plot for the mean size frequency distribution and overlay them
g4 <- ggplot(mus_pos_build1, aes(x = as.numeric(location), fill = ssw_period, color = ssw_period)) +
  geom_density(alpha = 0.8, adjust = 1.5) +
  labs(x = "Mussel depth \n(distance high to low intertidal)", y = "Frequency", title = "Depth distribution", tag = "B") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 10)) +
  scale_fill_manual(values=c("indianred","navyblue"))+
  scale_color_manual(values=c("indianred","navyblue"))+
  #geom_vline(data = mean_size_by_period, aes(xintercept = mean_size, color = ssw_period), 
  #          linetype = "dotted", size = 1) +
  theme_bw() + my_theme
g4

m <- ggpubr::ggarrange(g4, g2,g3, common.legend = TRUE, nrow=1, legend = "none")



n <- gridExtra::grid.arrange(g1_full, m, nrow=2, heights = c(0.65,0.25))



ggsave(n, filename = file.path(figdir, "Fig3_mussel_expansion.png"), 
      width = 7, height = 7.5, units = "in", dpi = 600) #last write 26 Sept 2024


################################################################################
#plot boxplot of mussel distribution

base_theme <-  theme(axis.text=element_text(size=11, color = "black"),
                     axis.title=element_text(size=12,color = "black"),
                     plot.tag=element_text(size=9,color = "black"),
                     plot.title=element_text(size=10,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.5, "cm"), 
                     #legend.key = element_rect(fill = "white"), # Set it to transparent
                     legend.spacing.y = unit(0.5, "cm"),  
                     legend.text=element_text(size=10,color = "black"),
                     legend.title=element_blank(),
                     #legend.key.height = unit(0.1, "cm"),
                     #legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=10, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())


# plot depth distribution
s1 <- ggplot(heatmap_data, aes(x = Period, y = location, fill = Period)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = c("Pre-SSW" = "navyblue", "Post-SSW" = "indianred")) +
  labs(title = "",
       x = "Period",
       y = "Distance (meters from high intertidal)",
       fill = "Period",
       tag = "A") +
  theme_minimal() +
  ggsignif::geom_signif(comparisons = list(c("Pre-SSW", "Post-SSW")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  theme_bw() + base_theme

# plot size
s2 <- ggplot(DatU %>%
               rename(Period = ssw_period) %>%
               mutate(Period = factor(Period, 
                                      levels = c("Before", "After"), 
                                      labels = c("Pre-SSW", "Post-SSW"))), 
             aes(x = Period, y = size_bin, fill = Period)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = c("Pre-SSW" = "navyblue", "Post-SSW" = "indianred")) +
  labs(title = "",
       x = "Period",
       y = "Mussel size (mm)",
       fill = "Period",
       tag = "B") +
  theme_minimal() +
  ggsignif::geom_signif(comparisons = list(c("Pre-SSW", "Post-SSW")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  theme_bw() + base_theme

s2

s <- ggpubr::ggarrange(s1, s2, nrow=1)

ggsave(s, filename = file.path(figdir, "FigS1_mussel_boxplots.png"), 
       width =7, height = 5, units = "in", dpi = 600, bg = "white") #last write 26 Sept 2024










