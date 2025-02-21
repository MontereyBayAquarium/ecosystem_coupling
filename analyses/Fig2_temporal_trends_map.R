
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

#required packages
librarian::shelf(tidyverse, sf, zoo, minpack.lm, ggsignif)

#set directories 
localdir <- here::here("output")
remdir <- "/Volumes/seaotterdb$/kelp_recovery/data/" #note: this path is to a MBA server
figdir <- here::here("figures")

#read census data
census_orig <- read_csv(file.path(localdir, "processed/census/census_data_processed.csv")) 

# Get rocky intertidal data
load(file.path(localdir, "processed/rocky_intertidal/pisaster_mytilus_processed.rdata"))

#read foraging data
forage_orig <- read_csv(file.path(localdir,"processed/foraging_data/foraging_data_2016_2023.csv"))

#read bathy
#Note: this file is too large for GitHub and is currently stored on a MBA server. 
#authorized users can access the server, or the data can be obtained directly from
#https://geodata.mit.edu/catalog/stanford-pb497kd3664
#bathy_5m <- st_read(file.path(remdir, "gis_data/raw/bathymetry/contours_5m/contours_5m.shp")) %>% filter(CONTOUR == "-5")

#read state
ca_counties_orig <- st_read(file.path(localdir, "raw/gis_data/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################

################################################################################
#Step 1 - prep data summaries 

#calculate total number of sea otters for the study area by year
#the data extent is already filtered to the study area
census_join_summary <- census_orig %>%
  group_by(year) %>%
  summarize(total_n_indep = sum(n_indep)) %>% 
  pivot_longer(cols = c(total_n_indep), 
               names_to = "category", 
               values_to = "response")

#mussels
mus_build2 <- mus_build1 %>% dplyr::select(year, site = marine_site_name, response = percent_cover) %>%
  #calculate mean percent cover
  group_by(year) %>% summarize(response = mean(response)) %>% mutate(category = "mussels") 

#stars
pis_build2 <- pis_build1 %>% dplyr::select(year, site = marine_site_name, response = density_per_m2) %>%
  #calculate mean density
  group_by(year) %>% summarize(response = mean(response)) %>%   mutate(category = "P. ochraceus") 

# Join the data
combined_data <- bind_rows(census_join_summary, mus_build2, pis_build2) %>%
                  #assign identifier for the three variables
                  mutate(category = factor(category)) %>%
                  filter(year >=2000)%>%
                  filter(!(category %in% c("total_n_pup"))) %>%
                  mutate(category = case_when(
                    category == "total_n_indep" ~ "Total independent otters",
                    category == "mussels" ~ "Mussels",
                    TRUE ~ category
                  ),
                  #set order
                  category = factor(category, levels = c("Total independent otters","Mussels","P. ochraceus")),
                  period = ifelse(year <= 2013, "Before","After"),
                  period = factor(period, levels = c("Before","After"))
                  )
                  
################################################################################
#Step 2 - extract mussel dives from foraging data 

forage_build1 <- forage_orig %>% 
  #filter mussel dives
  filter(prey == "mus") %>%
  #select bout locations
  dplyr::select(year, month, day, bout, lat, long) %>% distinct() %>%
  filter(!is.na(lat))%>%
  st_as_sf(coords = c("long","lat"), crs = 4326)


################################################################################
#plot

# Theme
base_theme <-  theme(axis.text=element_text(size=10, color = "black"),
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


# Build inset
g1_inset <-  ggplotGrob(
  ggplot() +
    # Plot land
    geom_sf(data=foreign, fill="grey80", color="white", lwd=0.3) +
    geom_sf(data=usa, fill="grey80", color="white", lwd=0.3) +
    # Plot box
    annotate("rect", xmin=-122.6, xmax=-121, ymin=36.2, ymax=37.1, color="black", fill=NA, lwd=0.6) +
    # Label regions
    #geom_text(data=region_labels, mapping=aes(y=lat_dd, label=region), x= -124.4, hjust=0, size=2) +
    # Labels
    labs(x="", y="") +
    # Crop
    coord_sf(xlim = c(-124.5, -117), ylim = c(32.5, 42)) +
    # Theme
    theme_bw() + base_theme +
    theme( plot.margin = unit(rep(0, 4), "null"),
           panel.margin = unit(rep(0, 4), "null"),
           panel.background = element_rect(fill='transparent'), #transparent panel bg
           # plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
           axis.ticks = element_blank(),
           axis.ticks.length = unit(0, "null"),
           axis.ticks.margin = unit(0, "null"),
           axis.text = element_blank(),
           axis.title=element_blank(),
           axis.text.y = element_blank())
)

# Create the "Monterey" text label
monterey_label <- data.frame(
  x = c(-121.9, -121.97), # x-coordinate for the upper right corner
  y = c(36.64, 36.54),  # y-coordinate for the upper right corner
  label = c("Monterey \nBay", "Carmel \nBay")
)

#rocky intertidal sites
rocky_sites <- pis_build1 %>% dplyr::select(monitoring_site = marine_site_name, latitude, longitude)%>%distinct()


g <- ggplot() + 
  # Add landmarks
  geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
            size = 3, fontface = "bold") +
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, 
                    xmax = -121.96,
                    ymin = 36.625) +
  # Add foraging bouts for mussels
  geom_sf(
    data = forage_build1, #%>% filter(year == 2017),
    aes(fill = "Mussel forage \nbout"),
    size=0.1,
    shape = 19,
    show.legend = FALSE
  ) + 
  guides(
    fill = guide_legend(
      override.aes = list(size = 3),  # Adjust the legend point size as needed
      title = NULL  # Remove the legend title
    )
  ) + 
  # Add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") + 
  # Add rocky intertidal sites
  geom_point(
    data = rocky_sites,
    aes(x = longitude, y = ifelse(monitoring_site == "Point Lobos",latitude+.001,latitude)),
    shape = 24,  
    size = 3,
    fill = "orange",
    show.legend = FALSE
  ) + 
  # Add scale bar
  ggsn::scalebar(x.min = -121.99, x.max = -121.88, 
                 y.min = 36.519, y.max = 36.645,
                 location="bottomright",
                 dist = 2, dist_unit = "km",
                 transform=TRUE, 
                 model = "WGS84",
                 st.dist=0.02,
                 st.size=3,
                 border.size=.5,
                 height=.02
  ) + 
  # Add north arrow
  ggsn::north(x.min = -121.99, x.max = -121.88, 
              y.min = 36.519, y.max = 36.65,
              location = "topright", 
              scale = 0.05, 
              symbol = 10) + 
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.515, 36.645), crs = 4326) + 
  labs(title = "", tag = "D") + 
  theme_bw() + base_theme + theme(
    axis.title = element_blank(),
    legend.position=c(.8,.2),
    legend.background = element_rect(fill='transparent'),
    legend.key = element_rect(fill='transparent'),
    #axis.text.x = element_text(angle = 90, vjust = 0.5),  # Rotate x-axis labels and center them
    axis.text.y = element_text(angle = 90, hjust = 0.5),  # Rotate y-axis labels and center them
    axis.ticks.x = element_line(),  # Show ticks for x-axis
    axis.ticks.y = element_line()   # Show ticks for y-axis
  ) + 
  # Manually set breaks for every other tick mark on the x-axis
  scale_x_continuous(breaks = seq(-121.98, -121.88, by = 0.03)) +  # Adjust for your x-axis range
  scale_y_continuous(breaks = seq(36.515, 36.645, by = 0.03))  # Adjust for your y-axis range
g

#from https://stackoverflow.com/questions/70977700/creating-completely-customized-legends-in-ggplot2
dummy_guide <- function(
    labels = NULL,  
    ..., 
    title = NULL, 
    key   = draw_key_point,
    guide_args = list()
) {
  # Capture arguments
  aesthetics <- list(...)
  n <- max(lengths(aesthetics), 0)
  labels <- labels %||%  seq_len(n)
  
  # Overrule the alpha = 0 that we use to hide the points
  aesthetics$alpha <- aesthetics$alpha %||% rep(1, n)
  
  # Construct guide
  guide_args$override.aes <- guide_args$override.aes %||% aesthetics
  guide <- do.call(guide_legend, guide_args)
  
  # Allow dummy aesthetic
  update_geom_defaults("point", list(dummy = "x"))
  
  dummy_geom <- geom_point(
    data = data.frame(x = rep(Inf, n), y = rep(Inf, n), 
                      dummy = factor(labels)),
    aes(x, y, dummy = dummy), alpha = 0, key_glyph = key
  )
  dummy_scale <- discrete_scale(
    "dummy", "dummy_scale", palette = scales::identity_pal(), name = title,
    guide = guide
  )
  list(dummy_geom, dummy_scale)
}

g0 <- g + dummy_guide(
  labels = c("Intertidal site", "Mussel foraging"),
  shape = c(24, 19),
  fill = c("orange", "black"),
  size=4
)

#g0

################################################################################
#plot pisaster 

# Separate data into three segments
segment1 <- combined_data %>% filter(category == "P. ochraceus" & year <= 2007)
segment2 <- combined_data %>% filter(category == "P. ochraceus" & year > 2007 & year <= 2013)
segment3 <- combined_data %>% filter(category == "P. ochraceus" & year >= 2013)

# Gompertz growth function
gompertz_growth <- function(year, A, B, C) {
  A * exp(-B * exp(-C * (year - min(segment1$year))))
}

# Adjust initial values
initial_values <- list(A = max(segment1$response, na.rm = TRUE), 
                       B = 1, 
                       C = 0.1)

# Fit Gompertz growth model with constraint
trend_model <- nlsLM(response ~ gompertz_growth(year, A, B, C), 
                     data = segment1, 
                     start = initial_values,
                     lower = c(0, -Inf, -Inf), 
                     upper = c(0.3, Inf, Inf),
                     control = list(maxiter = 200))

# Calculate residuals and standard deviation of residuals
predicted_segment1 <- predict(trend_model, newdata = segment1)
residuals <- segment1$response - predicted_segment1
residual_sd <- sd(residuals, na.rm = TRUE)

# Predict responses for the gap years (2007-2013)
gap_years <- data.frame(year = 2007:2013)
predicted_responses <- predict(trend_model, newdata = gap_years)

# Add noise to the predicted responses
set.seed(1985)  
predicted_responses_with_error <- predicted_responses + rnorm(length(predicted_responses), mean = 0, sd = residual_sd)

# Combine gap years data with predicted responses with error
gap_years <- cbind(gap_years, response = predicted_responses_with_error)

# Combine data for plotting

segment1 <- segment1 %>% mutate(segment = "segment1")
gap_years <- gap_years %>% mutate(segment = "segment2")
segment3 <- segment3 %>% mutate(segment = "segment3")

plot_data <- bind_rows(segment1, gap_years, segment3) %>%
  mutate(category = "P. ochraceus",
         response = ifelse(year == 2007 & segment == "segment2",0.236363636,response))


# Create the initial plot with geom_smooth
p <- ggplot(plot_data, aes(x = year, y = response)) +
  # Points for segments other than segment2
  geom_point(
    data = plot_data %>% filter(segment != "segment2"),
    alpha = 0.4,
    color = "#E377C2",
    size=1
  ) +
  # Points for segment2 in black, excluding year 2007
  geom_point(
    data = plot_data %>% filter(segment == "segment2", year != 2007),
    alpha = 0.4,
    color = "black",
    size=1
  ) +
  # Smooth line with high 'n' and no visible line or ribbon
  geom_smooth(
    data = plot_data %>% filter(segment == "segment1" | segment == "segment2"),
    method = "auto",
    color = NA,
    fill = NA,
    alpha = 0,
    n = 10000,
    lwd=0.7
  )

# Extract smoothed data including confidence intervals
smoothed_data <- ggplot_build(p)$data[[3]][, c("x", "y", "ymin", "ymax")]

# Add period indicator
smoothed_data <- smoothed_data %>%
  mutate(
    period = ifelse(x >= 2007 & x <= 2013, "2007-2013", "Other")
  )

# Add smoothed line and error ribbon
p1 <- p +
  # Add error ribbon with variable fill
  geom_ribbon(
    data = smoothed_data,
    aes(x = x, ymin = ymin, ymax = ymax, fill = period),
    alpha = 0.2,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  # Add smoothed line with conditional aesthetics
  geom_line(
    data = smoothed_data,
    aes(x = x, y = y, color = period, linetype = period),
    inherit.aes = FALSE,
    show.legend = FALSE,
    lwd=0.7
  ) +
  # Smooth line for segment3
  geom_smooth(
    data = plot_data %>% filter(segment == "segment3"),
    method = "lm",
    formula = y ~ poly(x, 1),
    color = "#E377C2",
    fill = "#E377C2",
    alpha = 0.2,
    show.legend = FALSE,
    lwd=0.7
  ) +
  # Vertical line and annotations
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  annotate(
    geom = "text",
    label = "SSW",
    x = 2010.5,
    y = 0.18,
    size = 2.5
  ) +
  annotate(
    "segment",
    x = 2010.8,
    y = 0.17,
    xend = 2012.7,
    yend = 0.1,
    arrow = arrow(type = "closed", length = unit(0.02, "npc"))
  ) +
  # Themes and scales
  theme_minimal() +
  theme_bw() +
  base_theme +
  scale_y_continuous(
    limits = c(0, NA),
    oob = scales::squish
  ) +
  scale_x_continuous(limits = c(2000, 2023)) +
  # Manual scales for color, fill, and linetype
  scale_color_manual(
    values = c("2007-2013" = "grey", "Other" = "#E377C2"),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c("2007-2013" = "grey", "Other" = "#E377C2"),
    guide = "none"
  ) +
  scale_linetype_manual(
    values = c("2007-2013" = "dashed", "Other" = "solid"),
    guide = "none"
  ) +
  labs(
    x = "Year",
    y = expression(paste("Density ",italic("(n "), "per m²)")),
    title = "Pisaster",
    tag = "C"
  ) +
  theme(
    plot.title = element_text(face = "bold.italic"),
    legend.title = element_blank()
  )

p1



# Filter data
plot_data <- combined_data %>% filter(category == "P. ochraceus")

# Perform t-test
t_test <- t.test(response ~ period, data = plot_data)

# Categorize p-value for display
p_label <- case_when(
  t_test$p.value < 0.001 ~ "p < 0.001",
  t_test$p.value < 0.01  ~ "p < 0.01",
  t_test$p.value < 0.05  ~ "p < 0.05",
  TRUE                   ~ "n.s."  # Not significant
)


# Apply the buffer to y-axis limits
p2 <- ggplot(plot_data, aes(x = period, y = response, fill = period)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, color = "black", size=0.7) +
  geom_signif(comparisons = list(c("Before", "After")), 
              test = "t.test", 
              map_signif_level = TRUE, 
              tip_length = 0.02, 
              textsize = 3,
              y_position = max(plot_data$response)*1.1) +   
  theme_minimal() +
  theme_bw() +
  base_theme +
  scale_fill_manual(values = c("Before" = "#E377C2", "After" = "#E377C2")) +
  labs(x = "Period", y = "", title = "") +
  #annotate("text", x = 2.8, y = max(plot_data$response)*1.2, 
   #        label = p_label, size = 2.5, hjust = 1) +  # Adds formatted p-value
  scale_y_continuous(limits = c(min(plot_data$response), max(plot_data$response)*1.25)) +  # Adjusted y-axis limits
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

#p2



p <- ggpubr::ggarrange(p1, p2, widths = c(1.9,1))
p



################################################################################
#plot mussels

m1 <- ggplot(combined_data %>% filter(category == "Mussels"), aes(x = year, y = response)) +
  geom_point(alpha = 0.4, color = "#FF7F0E", show.legend = FALSE, size=1) +
  stat_smooth(geom = "line", span = 0.6, color = "#FF7F0E", show.legend = FALSE, lwd=0.7) +
  stat_smooth(
    method = "loess",
    geom = "ribbon",
    alpha = 0.2,
    fill = "#FF7F0E",
    color = NA,
    span = 0.6,
    show.legend = FALSE,
    lwd=1
  ) +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  labs(x = "Year", y = "Percent cover") +
  ggtitle("") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme_bw() +
  base_theme +
  scale_y_continuous(limits = c(0, NA), oob = scales::squish) +
  scale_x_continuous(limits = c(2000, 2023)) +
  scale_color_manual(values = "#FF7F0E") + 
  scale_fill_manual(values = "#FF7F0E")+
  labs(x = "", y = "Cover (%)", title = "Mussels", tag = "B") +
  theme(axis.text.x = element_blank())

#m1


# Filter data
plot_data <- combined_data %>% filter(category == "Mussels")

# Perform t-test
t_test <- t.test(response ~ period, data = plot_data)

# Categorize p-value for display
p_label <- case_when(
  t_test$p.value < 0.001 ~ "p < 0.001",
  t_test$p.value < 0.01  ~ "p < 0.01",
  t_test$p.value < 0.05  ~ "p < 0.05",
  TRUE                   ~ "n.s."  # Not significant
)


# Boxplot with significance bracket
m2 <- ggplot(plot_data, aes(x = period, y = response, fill = period)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, color = "black", size = 0.7) +
  geom_signif(comparisons = list(c("Before", "After")), 
              test = "t.test", 
              map_signif_level = TRUE, 
              tip_length = 0.02, 
              textsize = 3,
              y_position = max(plot_data$response)*1.1) + 
  theme_minimal() +
  theme_bw() +
  base_theme +
  scale_fill_manual(values = c("Before" = "#FF7F0E", "After" = "#FF7F0E")) +
  labs(x = "Period", y = "", title = "") +
 # annotate("text", x = 2.8, y = max(plot_data$response)*1.2, 
  #         label = p_label, size = 2.5, hjust = 1) +  # Adds formatted p-value
  scale_y_continuous(limits = c(min(plot_data$response), max(plot_data$response)*1.2)) +  # Adjusted y-axis limits
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

m2

m <- ggpubr::ggarrange(m1, m2, widths = c(1.9,1))
m


################################################################################

#calculate no. otters by time period

n_otters <- combined_data %>% filter(category == "Total independent otters") %>%
            mutate(Period = ifelse(year < 2013, "Before", "After")) %>%
            group_by(Period) %>%
            summarize(mean = mean(response, na.rm=TRUE),
                      sd = sd(response, na.rm=TRUE))

#plot sea otters

o1 <- ggplot(combined_data %>% filter(category == "Total independent otters"), aes(x = year, y = response)) +
  geom_point(alpha = 0.4, color = "#2CA02C", show.legend = FALSE, size=1) +
  stat_smooth(geom = "line", span = 0.6, color = "#2CA02C", show.legend = FALSE, lwd=0.7) +
  stat_smooth(
    method = "loess",
    geom = "ribbon",
    alpha = 0.2,
    fill = "#2CA02C",
    color = NA,
    span = 0.6,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  ggtitle("") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme_bw() +
  base_theme +
  scale_y_continuous(limits = c(200, NA), oob = scales::squish) +
  scale_x_continuous(limits = c(2000, 2023)) +
  scale_color_manual(values = "#2CA02C") + 
  scale_fill_manual(values = "#2CA02C")  +    
  labs(x = "", y = expression(paste(italic("n "), "independents")),
       title = "Sea otters", tag = "A")+
  theme(axis.text.x = element_blank())

o1

# Filter data
plot_data <- combined_data %>% filter(category == "Total independent otters")

# Perform t-test
t_test <- t.test(response ~ period, data = plot_data)

# Categorize p-value for display
p_label <- case_when(
  t_test$p.value < 0.001 ~ "p < 0.001",
  t_test$p.value < 0.01  ~ "p < 0.01",
  t_test$p.value < 0.05  ~ "p < 0.05",
  TRUE                   ~ "n.s."  # Not significant
)


# Boxplot with significance bracket
o2 <- ggplot(plot_data, aes(x = period, y = response, fill = period)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, color = "black", size = 0.7) +
  geom_signif(comparisons = list(c("Before", "After")), 
              test = "t.test", 
              map_signif_level = TRUE, 
              tip_length = 0.02, 
              textsize = 3, 
              y_position = max(plot_data$response)*1.1) + 
  theme_minimal() +
  theme_bw() +
  base_theme +
  scale_fill_manual(values = c("Before" = "#2CA02C", "After" = "#2CA02C")) +
  #annotate("text", x = 2.8, y = max(plot_data$response)*1.2, 
   #        label = p_label, size = 2.5, hjust = 1) +  # Adds formatted p-value
  scale_y_continuous(limits = c(min(plot_data$response), max(plot_data$response)*1.2)) +  # Adjusted y-axis limits
  labs(x = "Period", y = "", title = "") +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

o2

o <- ggpubr::ggarrange(o1, o2, widths = c(1.9,1))
o



#combine plots

f <- gridExtra::grid.arrange(o,m,p, ncol=1) 
p_final <- gridExtra::grid.arrange(f, g0, ncol = 2, widths = c(1.15,1))

p_final



#save
ggsave(p_final, filename = file.path(figdir, "Fig2_temporal_trendsv3.png"), 
       width = 8.5, height = 5.5, units = "in", dpi = 600, bg = "white") #last write 26 Sept 2024








