#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, plotly)

#set directories 
localdir <- here::here("output")
figdir <- here::here("figures")

# Get rocky intertidal position data
load(file.path(localdir, "processed/rocky_intertidal/position_data.rdata"))


################################################################################
#Step 1 - filter to focal study area


DatU <- vcdExtra::expand.dft(mus_size_build1, freq="total")

# Calculate the mean size by ssw_period
mean_size_by_period <- DatU %>%
  group_by(ssw_period) %>%
  summarize(mean_size = mean(as.numeric(size_bin)))


###############################################################################
#3D plot


# Calculate the mean of mean_z_rock_height for each combination of intertidal_sitename, transect, and y_planar
heatmap_data <- mus_elev %>%
  group_by(intertidal_sitename, transect, y_planar, ssw_period
           ) %>%
  summarize(
    MeanZRockHeight = mean(mean_z_rock_height, na.rm=TRUE),
    MytilusPresence = as.numeric(any(final_classification == "mytilus californianus"))
  ) %>%
  arrange(intertidal_sitename, transect, y_planar) %>%
  mutate(ssw_period = factor(ssw_period, levels = c("After","Before")))


hopkins <- plot_ly(scene = 'scene1') %>%
  add_trace(data = as.data.frame(heatmap_data) %>% filter(intertidal_sitename == "Hopkins" & MytilusPresence > 0),
            x = ~transect,
            y = ~y_planar,
            z = ~MeanZRockHeight + 0.05,
            color = ~ssw_period,
            mode = "markers",
            colors = c("indianred","navyblue"),
            type = "scatter3d",
            showlegend = FALSE
  ) %>%
  add_trace(data = heatmap_data %>% filter(intertidal_sitename == "Hopkins",
                                           !(is.na(y_planar))),
            x = ~transect,
            y = ~y_planar, 
            z = ~MeanZRockHeight,
            colorscale = 'Viridis',
            type = "mesh3d",
            intensity = ~MeanZRockHeight) %>% hide_colorbar()%>%
  #test aesthetics
  layout(scene = list(
    yaxis = list(title = "High to <br> low intertidal (m)"),
    zaxis = list(title = "Elevation (m)"),
    xaxis = list(title = "Distance across <br> site (m)" #, autorange = "reversed"
    ),
    showlegend = TRUE,
    aspect.mode = "manual",
    aspectratio = list(x = 1.4, y = 1, z = 1) 
  ),
  title = "Hopkins",
  margin = list(t = 40) ,
  legend = list(title = list(
    text = "<br>Mussel \npresence"))
  ) %>%
  colorbar(title = "Elevation (m)")

hopkins


pinos <- plot_ly(scene = 'scene2') %>%
  add_trace(data = as.data.frame(heatmap_data) %>% filter(intertidal_sitename == "Point Pinos" & MytilusPresence > 0),
            x = ~transect,
            y = ~y_planar,
            z = ~MeanZRockHeight + 0.05,
            color = ~ssw_period,
            mode = "markers",
            colors = c("indianred","navyblue"),
            type = "scatter3d"
  ) %>%
  add_trace(data = heatmap_data %>% filter(intertidal_sitename == "Point Pinos",
                                           !(is.na(y_planar))),
            x = ~transect,
            y = ~y_planar, 
            z = ~MeanZRockHeight,
            colorscale = 'Viridis',
            type = "mesh3d",
            intensity = ~MeanZRockHeight) %>% hide_colorbar()%>%
  #test aesthetics
  layout(scene = list(
    yaxis = list(title = "High to <br> low intertidal (m)"),
    zaxis = list(title = "Elevation (m)"),
    xaxis = list(title = "Distance across <br> site (m)" #, autorange = "reversed"
    ),
    showlegend = TRUE,
    aspect.mode = "manual",
    aspectratio = list(x = 1.4, y = 1, z = 1) 
  ),
  title = "Point Pinos",
  margin = list(t = 40) ,
  legend = list(title = list(
    text = "<br>Mussel \npresence"))
  ) %>%
  colorbar(title = "Elevation (m)")

pinos


stillwater <- plot_ly(scene = 'scene3') %>%
  add_trace(data = as.data.frame(heatmap_data) %>% filter(intertidal_sitename == "Stillwater" & MytilusPresence > 0),
            x = ~transect,
            y = ~y_planar,
            z = ~MeanZRockHeight + 0.05,
            color = ~ssw_period,
            mode = "markers",
            colors = c("indianred","navyblue"),
            type = "scatter3d",
            showlegend=FALSE
  ) %>%
  add_trace(data = heatmap_data %>% filter(intertidal_sitename == "Stillwater"),
            x = ~transect,
            y = ~y_planar, 
            z = ~MeanZRockHeight,
            colorscale = 'Viridis',
            type = "mesh3d",
            intensity = ~MeanZRockHeight) %>% hide_colorbar()%>%
  #test aesthetics
  layout(scene3 = list(
    yaxis = list(title = "High to <br> low intertidal (m)"),
    zaxis = list(title = "Elevation (m)"),
    xaxis = list(title = "Distance across <br> site (m)" #, autorange = "reversed"
                 ),
    showlegend = TRUE,
    aspect.mode = "manual",
    aspectratio = list(x = 1.4, y = 1, z = 1) 
  ),
  title = "Stillwater",
  margin = list(t = 40) ,
  legend = list(title = list(
    text = "<br>Mussel \npresence"))
  ) %>%
  colorbar(title = "Elevation (m)")

stillwater


lobos <- plot_ly(scene='scene4') %>%
  add_trace(data = as.data.frame(heatmap_data) %>% filter(intertidal_sitename == "Point Lobos" & MytilusPresence > 0),
            x = ~transect,
            y = ~y_planar,
            z = ~MeanZRockHeight + 0.05,
            color = ~ssw_period,
            mode = "markers",
            colors = c("indianred","navyblue"),
            type = "scatter3d",
            showlegend=FALSE
  ) %>%
  add_trace(data = heatmap_data %>% filter(intertidal_sitename == "Point Lobos",
                                           !(is.na(y_planar))),
            x = ~transect,
            y = ~y_planar, 
            z = ~MeanZRockHeight,
            colorscale = 'Viridis',
            type = "mesh3d",
            intensity = ~MeanZRockHeight,
            showlegend=FALSE
            ) %>% hide_colorbar()%>%
  #test aesthetics
  layout(scene = list(
    yaxis = list(title = "High to <br> low intertidal (m)"),
    zaxis = list(title = "Elevation (m)"),
    xaxis = list(title = "Distance across <br> site (m)" #, autorange = "reversed"
    ),
    showlegend = TRUE,
    aspect.mode = "manual",
    aspectratio = list(x = 1.4, y = 1, z = 1) 
  ),
  title = "Point Lobos",
  margin = list(t = 40) ,
  legend = list(title = list(
    text = "<br>Mussel \npresence"))
  ) %>%
  colorbar(title = "Elevation (m)") 

lobos


#make scene


# Define custom grid style
axx <- list(
  gridcolor='rgb(255, 255, 255)',
  zerolinecolor='rgb(255, 255, 255)',
  showbackground=TRUE,
  backgroundcolor='rgb(230, 230,230)'
)

# Create subplot with separate scenes
fig <- subplot(hopkins, pinos, stillwater, lobos)


# Set up subplot layout with additional margin space
fig <- fig %>% layout(
  title = "",
  scene = list(
    domain = list(x = c(0, 0.5), y = c(0.5, 1)),  # Top-left panel
    xaxis = axx, yaxis = axx, zaxis = axx,
    aspectmode = 'cube',
    margin = list(l = 50, r = 50, t = 50, b = 100)  # Add margin space
  ),
  scene2 = list(
    domain = list(x = c(0.5, 1), y = c(0.5, 1)),  # Top-right panel
    xaxis = axx, yaxis = axx, zaxis = axx,
    aspectmode = 'cube',
    margin = list(l = 50, r = 50, t = 50, b = 100)  # Add margin space
  ),
  scene3 = list(
    domain = list(x = c(0, 0.5), y = c(0, 0.5)),  # Bottom-left panel
    xaxis = axx, yaxis = axx, zaxis = axx,
    aspectmode = 'cube',
    margin = list(l = 50, r = 50, t = 50, b = 100)  # Add margin space
  ),
  scene4 = list(
    domain = list(x = c(0.5, 1), y = c(0, 0.5)),  # Bottom-right panel
    xaxis = axx, yaxis = axx, zaxis = axx,
    aspectmode = 'cube',
    margin = list(l = 50, r = 50, t = 50, b = 100)  # Add margin space
  )
)

# Show the final subplot
fig













###############################################################################


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



ggsave(n, filename = file.path(figdir, "Fig2_mussel_expansionv2.png"), 
      width = 7, height = 7.5, units = "in", dpi = 600)



