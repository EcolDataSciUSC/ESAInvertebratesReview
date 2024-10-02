# This script visualizes threat maps to invertebrates across aquatic and
# terrestrial domains using data scraped from U.S. Endangered Species Act
# documentation and those species with Critical Habitat Designations.

#### SETUP ####
# Load libraries
library(colorspace)
library(cowplot)
library(data.table)
library(sf)
library(terra)
library(tidyverse)

sf::sf_use_s2(FALSE)
crsrobin <- "+proj=robin +lon_0=-180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
`%!in%` <- negate(`%in%`)

# Load in the basemap
basemap <- sf::st_read("../../../Resources/ShapefilesAndRasters/Countries/ne_10m_admin_0_countries.shp") %>%
  sf::st_crop(xmin = -180, xmax = -10, ymin = -15, ymax = 90)

basemap_GUCNMI <- sf::st_read("../../../Resources/ShapefilesAndRasters/Countries/ne_10m_admin_0_countries.shp") %>%
  dplyr::filter(NAME %in% c("Guam", "N. Mariana Is."))

basemap_AS <- sf::st_read("../../../Resources/ShapefilesAndRasters/Countries/ne_10m_admin_0_countries.shp") %>%
  dplyr::filter(NAME %in% c("American Samoa"))

basemap_union <- sf::st_read("../../../Resources/ShapefilesAndRasters/Countries/ne_10m_admin_0_countries.shp") %>%
  dplyr::filter(NAME %in% c("United States of America", "Canada", "Mexico",
                            "Guatemala", "El Salvador", "Belize", "Honduras",
                            "Nicaragua", "Costa Rica", "Panama", "Colombia", 
                            "Venezuela", "Guyana", "Suriname", "French Guiana",
                            "Cuba", "Cayman Islands", "Jamaica", "Haiti", "Dominican Rep.", "Bahamas",
                            "Puerto Rico", "British Virgin Islands", "Anguilla",
                            "Cayman Islands", "Aruba", "Trinidad and Tobago",
                            "Guadeloupe", "U.S. Virgin Is.", "Guam", "N. Mariana Is.", "American Samoa",
                            "Brazil", "Ecuador", "Peru", "Bolivia")) %>%
  dplyr::summarise(geometry = st_union(geometry)) %>%
  sf::st_crop(sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -15, ymax = 90),
                          crs = st_crs(.))) %>%
  sf::st_break_antimeridian(lon_0 = 170) %>%
  sf::st_transform(crsrobin)

# Create a hex grid over the basemap
my_hex <- sf::st_make_grid(basemap, cellsize = 2, flat_topped = TRUE, square = FALSE) %>%
  sf::st_as_sf() %>%
  dplyr::mutate(GID = row_number())
hex1_index <- nrow(my_hex)

my_hex_GUCNMI <- sf::st_make_grid(basemap_GUCNMI, cellsize = 2, flat_topped = TRUE, square = FALSE) %>%
  sf::st_as_sf() %>%
  dplyr::mutate(GID = row_number() + hex1_index)
hex2_index <- nrow(my_hex_GUCNMI)

my_hex_AS <- sf::st_make_grid(basemap_AS, cellsize = 2, flat_topped = TRUE, square = FALSE) %>%
  sf::st_as_sf() %>%
  dplyr::mutate(GID = row_number() + hex1_index + hex2_index)

my_hex <- do.call(rbind, list(my_hex, my_hex_GUCNMI, my_hex_AS))

# Load in the scraped threat assessment data
threat_dat <- data.table::fread("../data/FINAL_THREATS.csv") %>%
  dplyr::select(scientificName, 13:24)

threat_dat <- threat_dat %>%
  dplyr::group_by(scientificName) %>%
  dplyr::mutate(T1 = mean(t1_pred), T2 = mean(t2_pred),
                T3 = mean(t3_pred), T4 = mean(t4_pred),
                T5 = mean(t5_pred), T6 = mean(t6_pred),
                T7 = mean(t7_pred), T8 = mean(t8_pred),
                T9 = mean(t9_pred), T11 = mean(t11_pred), T12 = mean(t12_pred),
                T13 = mean(t13_pred)) %>%
  dplyr::ungroup() %>% 
  dplyr::select(scientificName, T1, T2, T4, T5,
                                     T6, T7, T8, T9, T11, T12) %>% unique()

threat_dat_binary <- threat_dat
threat_dat_binary[,2:11] <- lapply(threat_dat_binary[,2:11], 
                                   ceiling)

# Correlation plot of threats
corr_dat <- round(cor(threat_dat[,2:11]), 2)
p_dat <- ggcorrplot::cor_pmat(threat_dat[,2:11])
ggcorrplot::ggcorrplot(corr_dat, type = "lower", lab = TRUE, p.mat = p_dat,
                       colors = c("#E46726", "white", "#6D9EC1"),
                       insig = "blank",  show.diag = TRUE)+
  scale_x_discrete(labels = c("Res./Comm. Development", "Agriculture", "Energy Prod. & Mining",
                              "Bio. Res. Use", "Human Intrusion",
                              "Natural Sys. Mod.", "Problematic Biotic Factors", "Pollution",
                              "Climate Change", "Intrinsic Factors"))+
  scale_y_discrete(labels = c("Res./Comm. Development", "Agriculture", "Energy Prod. & Mining",
                              "Bio. Res. Use", "Human Intrusion",
                              "Natural Sys. Mod.", "Problematic Biotic Factors", "Pollution",
                              "Climate Change", "Intrinsic Factors"))+
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white", color = "white"))
cowplot::ggsave2("../figures/threat_correlations.png",
                 dpi = 400, height = 7, width = 7)

#### LOAD GBIF AND CRITICAL HABITAT INFORMATION ####
# Load in the Critical Habitat Designation data and get a count of species
crit_hab_sf <- sf::st_read("../data/ECOSRanges/crithab_poly.shp") %>%
  dplyr::inner_join(threat_dat_binary, by = c("sciname" = "scientificName")) %>%
  dplyr::mutate(value = 1) %>%
  sf::st_transform("wgs84")

# final_spp <- data.table::fread("../data/FINAL_THREATS_SPP.csv") %>%
#   dplyr::pull(scientificName) %>% rgbif::name_backbone_checklist() %>%
#   dplyr::filter(matchType != "NONE") %>%
#   dplyr::pull(usageKey)
# occ_query <- rgbif::occ_download(rgbif::pred_in("taxonKey", final_spp),
#                                  user = "vmshirey",
#                                  pwd = "18MacsNow",
#                                  email = "vmshirey@gmail.com")
# rgbif::occ_download_get("0032928-240906103802322", path = "../data/",
#                         overwrite = TRUE)

gbif_occ <- data.table::fread("../data/gbif.txt", header = TRUE, 
                              stringsAsFactors = FALSE, quote = "") %>%
  dplyr::filter(!is.na(decimalLongitude)) %>%
  dplyr::select(order, species, decimalLongitude, decimalLatitude) %>% unique() %>%
  sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude")) %>%
  dplyr::filter(species != "Danaus plexippus")
sf::st_crs(gbif_occ) <- "wgs84"

final_spp_domains <- data.table::fread("../data/FINAL_THREATS_SPP.csv")
aquatic_spp <- final_spp_domains %>%
  dplyr::filter(domain == "Aquatic")
terrestrial_spp <- final_spp_domains %>%
  dplyr::filter(domain == "Terrestrial")



##### AQUATICS MAP #####
# Parse out aquatic organisms from occurrence data
aquatic_occ <- gbif_occ %>%
  dplyr::filter(species %in% aquatic_spp$scientificName)

aquatic_occ <- sf::st_intersection(my_hex, aquatic_occ) %>%
  sf::st_drop_geometry() %>% unique()

aquatic_occ_spp <- aquatic_occ %>%
  dplyr::filter(!is.na(GID)) %>%
  dplyr::select(species, GID) %>%
  dplyr::group_by(GID) %>%
  dplyr::mutate(numSppGBIF = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(GID, numSppGBIF) %>% unique() %>%
  dplyr::mutate(GID = GID)

my_hex_aquatic <- my_hex %>%
  dplyr::left_join(aquatic_occ_spp, by = "GID")

# Now get a count of unique threats to species
my_hex_aquatic_final <- aquatic_occ %>%
  left_join(threat_dat_binary, by = c("species" = "scientificName")) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(across(T1:T12, max)) %>%
  dplyr::ungroup() %>% unique() %>%
  dplyr::right_join(my_hex_aquatic, by = "GID") %>%
  sf::st_as_sf() %>%
  dplyr::group_by(GID) %>%
  dplyr::summarise(across(T1:T12, sum))

my_hex_threatMapping <- my_hex_aquatic_final %>%
  sf::st_break_antimeridian(lon_0 = 170) %>%
  sf::st_transform(crsrobin)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T1)),
          mapping = aes(fill = T1, color = T1))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T1.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T2)),
          mapping = aes(fill = T2, color = T2))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T2.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T4)),
          mapping = aes(fill = T4, color = T4))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T4.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T5)),
          mapping = aes(fill = T5, color = T5))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T5.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T6)),
          mapping = aes(fill = T6, color = T6))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T6.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T7)),
          mapping = aes(fill = T7, color = T7))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T7.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T8)),
          mapping = aes(fill = T8, color = T8))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T8.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T9)),
          mapping = aes(fill = T9, color = T9))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T9.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T11)),
          mapping = aes(fill = T11, color = T11))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T11.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T11)),
          mapping = aes(fill = T11, color = T11))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T11.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T12)),
          mapping = aes(fill = T12, color = T12))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "Dark Mint",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "Dark Mint",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/aquatic_T12.png", dpi = 400,
                 heigh = 6, width = 10)

my_hex_aquatic_final_binary <- my_hex_aquatic_final %>%
  dplyr::mutate(across(T1:T12, ~ifelse(. >= 1, 1, 0))) %>%
  dplyr::mutate(numThreats = rowSums(across(T1:T12))) %>%
  dplyr::mutate(numThreats = ifelse(is.na(numThreats), 0, numThreats))

my_hex_aquatic$numThreats <- my_hex_aquatic_final_binary$numThreats

# Generate a bivariate map of number of species and unique threats
my_pal <- tibble("low-low" = "#f3f3f3",
                 "low-mid" = "#b3d3e0",
                 "low-high" = "#4f9dc1",
                 "mid-low" = "#f3e6b2",
                 "mid-mid" = "#b3b3b3",
                 "mid-high" = "#376386",
                 "high-low" = "#f4b301",
                 "high-mid" = "#b36500",
                 "high-high" = "#000000") %>%
  gather("bivariateCat", "fill")

my_hex_aquatic <- my_hex_aquatic %>%
  dplyr::mutate(numSppCat = cut(numSppGBIF, c(0, 20, 40, 60), include.lowest = TRUE,
                                labels = c("low", "mid", "high"))) %>%
  dplyr::mutate(numThreatsCat = cut(numThreats, c(0, 5, 10, 15), include.lowest = TRUE,
                                    labels = c("low", "mid", "high"))) %>%
  dplyr::mutate(bivariateCat = paste0(numSppCat, "-", numThreatsCat)) %>%
  dplyr::left_join(my_pal, by = "bivariateCat") %>%
  sf::st_break_antimeridian(lon_0 = 170) %>%
  sf::st_transform(crsrobin)

species_by_threat_map_lower48 <- ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(my_hex_aquatic,
          mapping = aes(fill = fill, color = fill))+
  geom_sf(dplyr::filter(my_hex_aquatic, !is.na(fill)),
          mapping = aes(), fill = NA, color = "white")+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  scale_fill_identity(na.value = NA)+
  scale_color_identity(na.value = NA)+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"))
cowplot::ggsave2("../figures/species_by_threat_mainland_aquatic.png",
                 species_by_threat_map_lower48, height = 6, width = 8, dpi = 400)

bivariate_legend <- my_pal %>%
  tidyr::separate(bivariateCat, into = c("spp", "thrt"), sep = "-") %>%
  dplyr::mutate(spp = ifelse(spp == "low", 1, ifelse(spp == "mid", 2, 3))) %>%
  dplyr::mutate(thrt = ifelse(thrt == "low", 1, ifelse(thrt == "mid", 2, 3)))

species_by_threat_map_lower48_legend <- ggplot()+
  geom_tile(bivariate_legend,
            mapping = aes(x = thrt, y = spp, fill = fill))+
  scale_fill_identity()+
  labs(x = "More Overlapping Threats →",
       y = "More Listed Species →")+
  theme_map()+
  theme(axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        plot.background = element_rect(fill = "white", color = "white"))+
  coord_fixed()
cowplot::ggsave2("../figures/species_by_threat_legend.png",
                 species_by_threat_map_lower48_legend, height = 6, width = 6, dpi = 400)

##### TERRESTRIALS MAP #####
# Parse out terrestrial organisms from occurrence data
terrestrial_occ <- gbif_occ %>%
  dplyr::filter(species %in% terrestrial_spp$scientificName)

terrestrial_occ <- sf::st_intersection(my_hex, terrestrial_occ) %>%
  sf::st_drop_geometry() %>% unique()

terrestrial_occ_spp <- terrestrial_occ %>%
  dplyr::filter(!is.na(GID)) %>%
  dplyr::select(species, GID) %>%
  sf::st_drop_geometry() %>% unique() %>%
  dplyr::group_by(GID) %>%
  dplyr::mutate(numSppGBIF = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(GID, numSppGBIF) %>% unique() %>%
  dplyr::mutate(GID = GID)

my_hex_terrestrial <- my_hex %>%
  dplyr::left_join(terrestrial_occ_spp, by = "GID")

# Now get a count of unique threats to species
my_hex_terrestrial_final <- terrestrial_occ %>%
  left_join(threat_dat_binary, by = c("species" = "scientificName")) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(across(T1:T12, max)) %>%
  dplyr::ungroup() %>% unique() %>%
  dplyr::right_join(my_hex_terrestrial, by = "GID") %>%
  sf::st_as_sf() %>%
  dplyr::group_by(GID) %>%
  dplyr::summarise(across(T1:T12, sum))

my_hex_threatMapping <- my_hex_terrestrial_final %>%
  sf::st_break_antimeridian(lon_0 = 170) %>%
  sf::st_transform(crsrobin)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T1)),
          mapping = aes(fill = T1, color = T1))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "heat",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "heat",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/terrestrial_T1.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T2)),
          mapping = aes(fill = T2, color = T2))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "heat",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "heat",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/terrestrial_T2.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T4)),
          mapping = aes(fill = T4, color = T4))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "heat",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "heat",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/terrestrial_T4.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T5)),
          mapping = aes(fill = T5, color = T5))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "heat",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "heat",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/terrestrial_T5.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T6)),
          mapping = aes(fill = T6, color = T6))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "heat",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "heat",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/terrestrial_T6.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T7)),
          mapping = aes(fill = T7, color = T7))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "heat",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "heat",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/terrestrial_T7.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T8)),
          mapping = aes(fill = T8, color = T8))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "heat",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "heat",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/terrestrial_T8.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T9)),
          mapping = aes(fill = T9, color = T9))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "heat",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "heat",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/terrestrial_T9.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T11)),
          mapping = aes(fill = T11, color = T11))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "heat",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "heat",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/terrestrial_T11.png", dpi = 400,
                 heigh = 6, width = 10)

ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(dplyr::filter(my_hex_threatMapping, !is.na(T12)),
          mapping = aes(fill = T12, color = T12))+
  colorspace::scale_fill_continuous_sequential(name = "Frequency",
                                               palette = "heat",
                                               guide = guide_colorbar(position = "inside"))+
  colorspace::scale_color_continuous_sequential(name = "Frequency",
                                                palette = "heat",
                                                guide = guide_colorbar(position = "inside"))+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position.inside = c(0.025, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(c(5,5,5,5)),
        legend.key.width = unit(1, "cm"))
cowplot::ggsave2("../figures/supp/terrestrial_T12.png", dpi = 400,
                 heigh = 6, width = 10)


my_hex_terrestrial_final_binary <- my_hex_terrestrial_final %>%
  dplyr::mutate(across(T1:T12, ~ifelse(. >= 1, 1, 0))) %>%
  dplyr::mutate(numThreats = rowSums(across(T1:T12))) %>%
  dplyr::mutate(numThreats = ifelse(is.na(numThreats), 0, numThreats))

my_hex_terrestrial$numThreats <- my_hex_terrestrial_final_binary$numThreats

# Generate a bivariate map of number of species and unique threats
my_pal <- tibble("low-low" = "#f3f3f3",
                 "low-mid" = "#b3d3e0",
                 "low-high" = "#4f9dc1",
                 "mid-low" = "#f3e6b2",
                 "mid-mid" = "#b3b3b3",
                 "mid-high" = "#376386",
                 "high-low" = "#f4b301",
                 "high-mid" = "#b36500",
                 "high-high" = "#000000")  %>%
  gather("bivariateCat", "fill")

my_hex_terrestrial <- my_hex_terrestrial %>%
  dplyr::mutate(numSppCat = cut(numSppGBIF, c(0, 5, 10, 15), include.lowest = TRUE,
                                labels = c("low", "mid", "high"))) %>%
  dplyr::mutate(numThreatsCat = cut(numThreats, c(0, 5, 10, 15), include.lowest = TRUE,
                                    labels = c("low", "mid", "high"))) %>%
  dplyr::mutate(bivariateCat = paste0(numSppCat, "-", numThreatsCat)) %>%
  dplyr::left_join(my_pal, by = "bivariateCat") %>%
  sf::st_break_antimeridian(lon_0 = 170) %>%
  sf::st_transform(crsrobin)

species_by_threat_map_lower48 <- ggplot()+
  geom_sf(basemap_union,
          mapping = aes(), fill = "white", color = "black")+
  geom_sf(my_hex_terrestrial,
          mapping = aes(fill = fill, color = fill))+
  geom_sf(dplyr::filter(my_hex_terrestrial, !is.na(fill)),
          mapping = aes(), fill = NA, color = "white")+
  geom_sf(basemap_union,
          mapping = aes(), fill = NA, color = "black")+
  scale_fill_identity(na.value = NA)+
  scale_color_identity(na.value = NA)+
  cowplot::theme_half_open()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"))
cowplot::ggsave2("../figures/species_by_threat_mainland_terrestrial.png",
                 species_by_threat_map_lower48, height = 6, width = 8, dpi = 400)

bivariate_legend <- my_pal %>%
  tidyr::separate(bivariateCat, into = c("spp", "thrt"), sep = "-") %>%
  dplyr::mutate(spp = ifelse(spp == "low", 1, ifelse(spp == "mid", 2, 3))) %>%
  dplyr::mutate(thrt = ifelse(thrt == "low", 1, ifelse(thrt == "mid", 2, 3)))

species_by_threat_map_lower48_legend <- ggplot()+
  geom_tile(bivariate_legend,
            mapping = aes(x = thrt, y = spp, fill = fill))+
  scale_fill_identity()+
  labs(x = "More Overlapping Threats →",
       y = "More Listed Species →")+
  theme_map()+
  theme(axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, angle = 90),
        plot.background = element_rect(fill = "white", color = "white"))+
  coord_fixed()
cowplot::ggsave2("../figures/species_by_threat_legend.png",
                 species_by_threat_map_lower48_legend, height = 6, width = 6, dpi = 400)
