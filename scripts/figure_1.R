###############################################################
# Figure1
# meta, cols and orden across panels
###############################################################

library(ape)
library(ggtree)
library(treeio)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(scatterpie)
library(scales)

setwd("/Users/hugo/Desktop/vp1_parachovirus")

###############################################################
# Shared objects
###############################################################

tree <- read.tree("vp1_headers_updated_v2.nwk")

orden <- paste0("HPeV",1:19)

meta <- tibble(label=tree$tip.label) %>%
  separate(label,
           into=c("Accession","Virus","Genotype","Country","Year"),
           sep="\\|",
           remove=FALSE,
           fill="right") %>%
  mutate(
    Year=as.numeric(Year),
    Genotype=factor(Genotype,levels=orden)
  )

cols <- c(
  "HPeV1"="#4E79A7","HPeV2"="#F28E2B","HPeV3"="#59A14F",
  "HPeV4"="#E15759","HPeV5"="#76B7B2","HPeV6"="#EDC948",
  "HPeV7"="#B07AA1","HPeV8"="#FF9DA7","HPeV9"="#9C755F",
  "HPeV10"="#BAB0AC","HPeV11"="#8CD17D","HPeV12"="#499894",
  "HPeV13"="#D37295","HPeV14"="#86BCB6","HPeV15"="#A0CBE8",
  "HPeV16"="#FFBE7D","HPeV17"="#8A89A6")

###############################################################
# PANEL A
###############################################################

p.tree <-
  ggtree(tree, size = .18) %<+% meta +
  
  geom_tippoint(
    aes(color = Genotype),
    size = 0.60,
    alpha = .96
  ) +
  scale_color_manual(
    values = cols,
    breaks = orden,
    drop = FALSE,
    name = "Genotype"
  ) +
  theme_tree() +
  guides(
    colour = guide_legend(
      ncol = 3,
      byrow = TRUE,
      keyheight = unit(0.35, "cm"),
      keywidth  = unit(0.35, "cm"),
      override.aes = list(size = 3)
    )
  ) +
  
  theme(
    legend.position = c(0.82, 0.42),
    legend.direction = "vertical",
    
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      colour = NA
    ),
    
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    
    legend.text = element_text(size = 9),
    legend.title = element_text(
      face = "bold",
      size = 11
    )
  )

p.tree <-
  p.tree +
  labs(tag = "A") +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    )
  )
p.tree
###############################################################
## PANEL B
## World map
###############################################################

library(sf)
library(rnaturalearth)
library(scatterpie)
library(ggforce)
library(dplyr)
library(tidyr)
library(stringr)

###############################################################
## Prepare metadata
###############################################################

meta2 <- meta %>%
  transmute(
    genotype = as.character(Genotype),
    country = Country
  ) %>%
  mutate(
    country = str_replace_all(country,"_"," "),
    country = str_trim(country),
    country = case_when(
      country=="USA" ~ "United States of America",
      country=="UK" ~ "United Kingdom",
      country=="Russian Federation" ~ "Russia",
      country=="Czech Republic" ~ "Czechia",
      country=="Bangladeshi" ~ "Bangladesh",
      country=="Viet Nam" ~ "Vietnam",
      country=="Cote d'Ivoire" ~ "Côte d'Ivoire",
      TRUE ~ country
    )
  ) %>%
  filter(country!="Unknown")

###############################################################
## Country x genotype matrix
###############################################################

country.matrix <-
  meta2 %>%
  count(country,genotype) %>%
  pivot_wider(
    names_from=genotype,
    values_from=n,
    values_fill=0
  ) %>%
  left_join(
    meta2 %>%
      count(country,name="Genomes"),
    by="country"
  )

###############################################################
## World map
###############################################################

world <-
  ne_countries(
    scale="medium",
    returnclass="sf"
  ) %>%
  dplyr::select(
    country=admin,
    geometry
  )

coord <-
  cbind(
    st_drop_geometry(
      st_point_on_surface(world)
    ),
    st_coordinates(
      st_point_on_surface(world)
    )
  ) %>%
  rename(
    X=X,
    Y=Y
  )

map.data <-
  left_join(
    country.matrix,
    coord,
    by="country"
  )

###############################################################
## Pie radius
###############################################################

map.data$radius <-
  scales::rescale(
    sqrt(map.data$Genomes),
    to=c(2,10)
  )

genotype.cols <-
  intersect(
    orden,
    colnames(map.data)
  )

###############################################################
## Base map
###############################################################

p.map <-
  ggplot() +
  geom_sf(
    data=world,
    fill="grey82",
    colour="white",
    linewidth=.20
  )+
  geom_scatterpie(
    data=map.data,
    aes(
      x=X,
      y=Y,
      r=radius
    ),
    cols=genotype.cols,
    colour="grey35",
    linewidth=.12)+
  scale_fill_manual(
    values=cols,
    guide="none")+
  coord_sf(
    xlim=c(-180,185),
    ylim=c(-60,90),
    expand=FALSE,
    clip="off") +
  theme_void()+
  theme(
    panel.background=
      element_rect(
        fill="white",
        colour=NA
      ),
    plot.background=
      element_rect(
        fill="white",
        colour=NA))

p.map
###############################################################
## Genome legend (inside map)
###############################################################

legend.df <-
  data.frame(
    genomes=c(5,10,25,50,100),
    radius=c(2.2,3.0,4.1,5.2,6.6),
    x=165,
    y=c(10,-2,-18,-40,-68))

###############################################################
## Draw legend
###############################################################

p.map <-
  p.map +
  geom_rect(
    xmin=150,
    xmax=184,
    ymin=-82,
    ymax=20,
    fill="white",
    colour="grey70",
    linewidth=.25
  )+
  geom_circle(
    data=legend.df,
    aes(
      x0=x,
      y0=y,
      r=radius
    ),
    fill="white",
    colour="black",
    linewidth=.30
  )+
  geom_text(
    data=legend.df,
    aes(
      x=x+11,
      y=y,
      label=genomes
    ),
    hjust=0,
    size=3
  )+
  annotate(
    "text",
    x=167,
    y=17,
    label="No. genomes",
    fontface="bold",
    size=3.8)

###############################################################
## Final panel
###############################################################

panel.B <- p.map

panel.B <-
  panel.B +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 18
      )
    )
  )

panel.B$patches$annotation$tag_levels <- NULL
panel.B$patches$annotation$tag <- "B"
panel.B
###############################################################
## PANEL C
## Sampling timeline
###############################################################

timeline <- meta %>%
  transmute(
    Country = Country,
    Year = as.numeric(Year)
  ) %>%
  filter(!is.na(Year))

###############################################################
## Clean country names
###############################################################

timeline <- timeline %>%
  mutate(
    Country = str_replace_all(Country, "_", " "),
    Country = case_when(
      Country == "USA" ~ "United States",
      Country == "UK" ~ "United Kingdom",
      Country == "Bangladeshi" ~ "Bangladesh",
      Country == "Czech Republic" ~ "Czechia",
      Country == "Viet Nam" ~ "Vietnam",
      TRUE ~ Country
    )
  )

###############################################################
## Order countries by total genomes
###############################################################

country.order <- timeline %>%
  count(Country, sort = TRUE) %>%
  pull(Country)

timeline$Country <- factor(
  timeline$Country,
  levels = rev(country.order)
)

###############################################################
## Total genomes per country and year
###############################################################

timeline.plot <- timeline %>%
  count(
    Country,
    Year,
    name = "Genomes"
  )

###############################################################
## Draw figure
###############################################################

p.time <-
  
  ggplot(
    timeline.plot,
    aes(
      x = Year,
      y = Country
    )
  ) +
  
  geom_point(
    aes(size = Genomes),
    colour = "black",
    alpha = 0.75
  ) +
  
  scale_size_continuous(
    name = "No. Genomes",
    breaks = c(1,5,10,20,40),
    range = c(2,8)
  ) +
  
  scale_x_continuous(
    breaks = seq(
      min(timeline.plot$Year),
      max(timeline.plot$Year),
      by = 1
    ),
    expand = c(0.01,0)
  ) +
  
  labs(
    x = "Collection year",
    y = NULL
  ) +
  
  theme_bw(base_size = 10) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7
    ),
    
    axis.text.y = element_text(size = 9),
    
    axis.title.x = element_text(
      face = "bold",
      size = 11
    ),
    
    panel.grid.major.x = element_line(
      colour = "grey88",
      linewidth = 0.25
    ),
    
    panel.grid.major.y = element_line(
      colour = "grey90",
      linewidth = 0.25
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "right"
    
  )

p.time

p.time <-
  p.time +
  labs(tag = "C") +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    )
  )

library(patchwork)

right <-
  panel.B /
  p.time +
  plot_layout(
    heights = c(1.6,0.75)
  )

fig1 <-
  p.tree | right +
  plot_layout(
    widths = c(1.4, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 18
      )
    )
  )

fig1
