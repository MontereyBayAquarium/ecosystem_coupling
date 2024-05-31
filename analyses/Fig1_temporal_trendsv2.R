
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

#required packages
librarian::shelf(tidyverse, sf, zoo)

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
                  category = factor(category, levels = c("Total independent otters","Mussels","P. ochraceus"))
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


# Theme
base_theme <-  theme(axis.text=element_text(size=8, color = "black"),
                     axis.title=element_text(size=8,color = "black"),
                     plot.tag=element_text(size=7,color = "black"),
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
  #add foraging bouts for mussels
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
  )+
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  #add rocky intertidal sites
  geom_point(
    data = rocky_sites,
    aes(x = longitude, y = ifelse(monitoring_site == "Point Lobos",latitude+.001,latitude)),
    shape = 24,  
    size = 3,
    fill = "orange",
    show.legend = FALSE
  )+
  #add scale bar
  ggsn::scalebar(x.min = -121.99, x.max = -121.88, 
                 y.min = 36.519, y.max = 36.645,
                 #anchor=c(x=-124.7,y=41),
                 location="bottomright",
                 dist = 2, dist_unit = "km",
                 transform=TRUE, 
                 model = "WGS84",
                 st.dist=0.02,
                 st.size=3,
                 border.size=.5,
                 height=.02
  )+
  #add north arrow
  ggsn::north(x.min = -121.99, x.max = -121.88, 
              y.min = 36.519, y.max = 36.65,
              location = "topright", 
              scale = 0.05, 
              symbol = 10)+
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.515, 36.645), crs = 4326)+
  labs(title = "", tag = "D")+
  theme_bw() + base_theme + theme(
    axis.title = element_blank(),
    legend.position=c(.8,.2),
    #legend.position = "top",  # Position the legend at the top
    #legend.justification = "right",  # Align the legend to the right
    #legend.box = "vertical",  # Box style for the legend
    #legend.margin = margin(t = -10, r = 10),  # Adjust the top and right margins for positioning
    axis.text = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.key = element_rect(fill='transparent')
  )
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

g0

#######OLD
p1 <- ggplot(combined_data %>% filter(category == "P. ochraceus"), aes(x = year, y = response)) +
  geom_point(alpha = 0.4, color = "#E377C2", show.legend = FALSE) +
  stat_smooth(geom = "line", size = 1, span = 0.6, color = "#E377C2", show.legend = FALSE) +
  stat_smooth(
    method = "loess",
    geom = "ribbon",
    alpha = 0.2,
    fill = "#E377C2",
    color = NA,
    span = 0.6,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  #SSW
  geom_vline(xintercept = 2013, linetype = "dotted", size=0.3)+
  annotate(geom="text", label="SSW", x=2011.5, y=0.25, size=2.5) +
  annotate("segment", x = 2011.8, y = 0.22, xend = 2012.7, yend = 0.2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme_bw() +
  base_theme +
  scale_y_continuous(limits = c(0, NA), oob = scales::squish) +
  scale_x_continuous(limits = c(2000, 2023)) +
  scale_color_manual(values = "#E377C2") + 
  scale_fill_manual(values = "#E377C2") +
  labs(x = "Year", y = "Density (no. per m²)", title = "Pisaster", tag = "C") +
  theme(plot.title = element_text(face = "bold.italic"))

#p1


# Separate data into three segments
segment1 <- combined_data %>% filter(category == "P. ochraceus" & year <= 2007)
segment2 <- combined_data %>% filter(category == "P. ochraceus" & year >= 2008 & year <= 2012)
segment3 <- combined_data %>% filter(category == "P. ochraceus" & year >= 2013)

# Transform data for asymptotic fitting
transformed_segment1 <- mutate(segment1, reciprocal_year = 1/year)

# Model trend using linear regression with reciprocal transformation
trend_model <- lm(response ~ poly(reciprocal_year, 2), data = transformed_segment1)

# Predict responses for the gap years (2008-2012)
gap_years <- data.frame(year = 2007:2013, reciprocal_year = 1/(2007:2013))
predicted_responses <- predict(trend_model, newdata = gap_years)

# Combine data for plotting
plot_data <- rbind(segment1, segment3)

# Plot with updated layers
p1 <- ggplot(plot_data, aes(x = year, y = response)) +
  geom_point(alpha = 0.4, color = "#E377C2", show.legend = FALSE) +
  geom_smooth(data = segment1, method = "lm", formula = y ~ poly(x, 2), size = 1, color = "#E377C2", fill = "#E377C2", alpha = 0.2, show.legend = FALSE) +
  geom_smooth(data = segment3, method = "lm", formula = y ~ poly(x, 1), size = 1, color = "#E377C2", fill = "#E377C2", alpha = 0.2, show.legend = FALSE) +
  geom_line(data = gap_years, aes(y = predicted_responses), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  annotate(geom="text", label="SSW", x=2010.5, y=0.28, size=2.5) +
  annotate("segment", x = 2010.8, y = 0.25, xend = 2012.7, yend = 0.1,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme_bw() +
  base_theme +
  scale_y_continuous(limits = c(0, NA), oob = scales::squish) +
  scale_x_continuous(limits = c(2000, 2023)) +
  scale_color_manual(values = "#E377C2") + 
  labs(x = "Year", y = "Density (no. per m²)", title = "Pisaster", tag = "C") +
  theme(plot.title = element_text(face = "bold.italic"))

#p1



p2 <- ggplot(combined_data %>% filter(category == "Mussels"), aes(x = year, y = response)) +
  geom_point(alpha = 0.4, color = "#FF7F0E", show.legend = FALSE) +
  stat_smooth(geom = "line", size = 1, span = 0.6, color = "#FF7F0E", show.legend = FALSE) +
  stat_smooth(
    method = "loess",
    geom = "ribbon",
    alpha = 0.2,
    fill = "#FF7F0E",
    color = NA,
    span = 0.6,
    show.legend = FALSE
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
  labs(x = "", y = "Percent cover", title = "Mussels", tag = "B") +
  theme(axis.text.x = element_blank())

#p2


p3 <- ggplot(combined_data %>% filter(category == "Total independent otters"), aes(x = year, y = response)) +
  geom_point(alpha = 0.4, color = "#2CA02C", show.legend = FALSE) +
  stat_smooth(geom = "line", size = 1, span = 0.6, color = "#2CA02C", show.legend = FALSE) +
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
  labs(x = "", y = "Number of independents",title = "Sea otters", tag = "A") +
  theme(axis.text.x = element_blank())

#p3


p <- gridExtra::grid.arrange(p3,p2,p1, ncol=1) 

p_final <- gridExtra::grid.arrange(p,g0, ncol=2)




#save
ggsave(p_final, filename = file.path(figdir, "Fig1_temporal_trendsv2.png"), 
       width =7, height = 5.5, units = "in", dpi = 600, bg = "white")


