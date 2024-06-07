#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

#required packages
librarian::shelf(tidyverse, here)

#set directories 
localdir <- here::here("output")
figdir <- here::here("figures")

# Get rocky intertidal position data
load(file.path(localdir, "processed/rocky_intertidal/position_datav2.rdata"))


################################################################################
#Step 1 - convert frequency data to long format and calculate mean

#expand frequency
DatU <- vcdExtra::expand.dft(mus_size_build1, freq="total")

# Calculate the mean size by ssw_period
mean_size_by_period <- DatU %>%
  group_by(ssw_period) %>%
  summarize(mean_size = mean(as.numeric(size_bin)))


################################################################################
#plot heatmap by site


tile_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                     axis.title=element_text(size=7,color = "black"),
                     plot.tag=element_text(size=7,color = "black"),
                     plot.title=element_text(size=7,color = "black", face = "bold"),
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
                     strip.text = element_text(size=6, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())

# Calculate the range of points for each site across all years and identify points unique to 'after'
#We use !any(ssw_period == "before") to check if there are no rows with 'ssw_period'
#equal to "before" within each group of 'intertidal_sitename,' 'transect,' and 'location.' 
#If there are no such rows, it means the points are unique to the 'after' period, and we label them as "Only After."

heatmap_data <- mus_pos_build1 %>%
  group_by(intertidal_sitename, transect, location) %>%
  summarize(
    location_range = max(location) - min(location),
    Period = factor(!any(ssw_period == "Before"), levels = c(TRUE, FALSE), labels = c("Post-SSW", "Pre-SSW"))
  ) %>%
  arrange(intertidal_sitename, transect) %>%
  group_by(intertidal_sitename) %>%
  mutate(transect_seq = dense_rank(transect),
         intertidal_sitename = factor(intertidal_sitename,levels = c("Hopkins","Point Pinos","Stillwater","Point Lobos")),
         Period = factor(Period, levels = c("Pre-SSW","Post-SSW")))


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
my_theme <-  theme(axis.text=element_text(size=6,color = "black"),
                   axis.title=element_text(size=6,color = "black"),
                   plot.tag=element_text(size=7,color = "black"),
                   plot.title=element_text(size=7,color = "black", face = "bold"),
                   # Gridlines
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key.size = unit(0.3, "cm"), 
                   #legend.key = element_rect(fill = "white"), # Set it to transparent
                   legend.spacing.y = unit(0.1, "cm"),  
                   legend.text=element_text(size=6,color = "black"),
                   legend.title=element_blank(),
                   #legend.key.height = unit(0.1, "cm"),
                   #legend.background = element_rect(fill=alpha('blue', 0)),
                   #facets
                   strip.text = element_text(size=6, face = "bold",color = "black", hjust=0),
                   strip.background = element_blank())

# Create a KDE plot for the mean size frequency distribution and overlay them
g2 <- ggplot(DatU, aes(x = as.numeric(size_bin), fill = ssw_period, color = ssw_period)) +
  geom_density(alpha = 0.8, adjust = 1.5) +
  labs(x = "Size (mm)", y = "Frequency", title = "Size frequency", tag = "C") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 120, by = 10)) +
  scale_fill_manual(values=c("indianred","navyblue"))+
  scale_color_manual(values=c("indianred","navyblue"))+
  #geom_vline(data = mean_size_by_period, aes(xintercept = mean_size, color = ssw_period), 
  #          linetype = "dotted", size = 1) +
  theme_bw() + my_theme + theme(axis.title.y = element_blank())
g2



# Create a KDE plot for the mean size frequency distribution and overlay them
g3 <- ggplot(mus_cov_period, aes(x = as.numeric(percent_cover), fill = ssw_period, color = ssw_period)) +
  geom_density(alpha = 0.8, adjust = 1.5) +
  labs(x = "Percent cover", y = "Frequency", title = "Percent cover", tag = "D") +
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
  labs(x = "Mussel depth (distance high to low intertidal)", y = "Frequency", title = "Depth distribution", tag = "B") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 10)) +
  scale_fill_manual(values=c("indianred","navyblue"))+
  scale_color_manual(values=c("indianred","navyblue"))+
  #geom_vline(data = mean_size_by_period, aes(xintercept = mean_size, color = ssw_period), 
  #          linetype = "dotted", size = 1) +
  theme_bw() + my_theme
g4

m <- ggpubr::ggarrange(g4, g2,g3, common.legend = TRUE, nrow=1, legend = "none")



n <- gridExtra::grid.arrange(g1_full, m, nrow=2, heights = c(0.65,0.25))



ggsave(n, filename = file.path(figdir, "Fig2_mussel_expansion.png"), 
      width = 7, height = 7.5, units = "in", dpi = 600)



