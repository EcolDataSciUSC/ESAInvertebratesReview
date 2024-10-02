# This script visualizes various aspects of the ESA as it applies to
# invertebrate species.

# Load libraries
library(colorspace)
library(cowplot)
library(data.table)
library(ggalluvial)
library(ggnewscale)
library(googlesheets4)
library(networkD3)
library(taxizedb)
library(tidyverse)
library(streamgraph)

`%!in%` <- negate(`%in%`)

#### FIGURE ONE ####
esa_dat <- fread("../data/summary_data.csv", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::mutate(itisGenus = word(ScientificName, start = 1))

# Calculate estimated number of species
gbif_list1 <- data.table::fread("../data/gbif_list1.csv") %>%
  dplyr::filter(taxonRank == "SPECIES",
                species != "",
                !is.na(species),
                species != " ",
                numberOfOccurrences >= 100) %>%
  dplyr::select(species) %>% unique() %>% arrange(species) 

gbif_list2 <- data.table::fread("../data/gbif_list2.csv") %>%
  dplyr::filter(taxonRank == "SPECIES",
                species != "",
                !is.na(species),
                species != " ",
                numberOfOccurrences >= 100) %>%
  dplyr::select(species) %>% unique() %>% arrange(species)

# Retrieve parent taxonomy from ITIS (note - need to clean up, especially at the subspecies
# level where ITIS doesn't have many names).
taxizedb::db_download_itis()
src_itis <- taxizedb::src_itis()

my_tax_ids <- taxizedb::name2taxid(esa_dat$ScientificName, db = "itis", out_type = "summary")
# parent_tax <- taxizedb::classification(my_tax_ids$id, db = "itis")
parent_tax <- readRDS("../data/parent_taxonomy.rds")

kingdom_tax <- do.call(rbind, (lapply(parent_tax, function(x) x[1,]))) %>%
  dplyr::mutate(tax_id = as.numeric(rownames(.))) %>%
  dplyr::left_join(my_tax_ids, by = c("tax_id" = "id"))
colnames(kingdom_tax) <- c("kingdom", "rank", "kingdom_id", "id", "scientificName")
kingdom_tax <- kingdom_tax %>%
  left_join(esa_dat, by = c("scientificName" = "ScientificName")) %>%
  dplyr::filter(kingdom == "Animalia")

class_tax <- do.call(rbind, (lapply(parent_tax, 
                                    function(x){
                                      dplyr::filter(x, x$rank == "class")
                             }))) %>%
  dplyr::mutate(tax_id = as.numeric(rownames(.))) %>%
  dplyr::left_join(my_tax_ids, by = c("tax_id" = "id"))
colnames(class_tax) <- c("class", "rank", "class_id", "id", "scientificName")
class_tax <- class_tax %>%
  left_join(esa_dat, by = c("scientificName" = "ScientificName")) %>%
  dplyr::filter(scientificName %in% kingdom_tax$scientificName) %>%
  dplyr::select(class, class_id, id, scientificName)

chordata_phylum_tax <- do.call(rbind, (lapply(parent_tax, function(x) x[4,]))) %>%
  dplyr::filter(rank == "phylum") %>%
  dplyr::mutate(tax_id = as.numeric(rownames(.)))

final_tax <- kingdom_tax %>% left_join(class_tax, by = "scientificName") %>%
  dplyr::mutate(id = id.x)

invert_tax <- final_tax %>%
  dplyr::filter(id %!in% chordata_phylum_tax$tax_id) %>%
  dplyr::mutate(higherTax = "Invertebrate")

# Randomly pull a subset of invertebrate documentation to be use to validate the use of
# ChatGPT for analysis and summary
invert_docs_subset <- invert_tax %>%
  dplyr::select(TaxonID, class, scientificName,
                SpeciesStatusAssessmentURLs, RecoveryPlanURLs,
                ReviewURLs, CriticalHabitatDesignationURLs) %>%
  tidyr::pivot_longer(ends_with("URLs"),
                      names_to = "documentType", values_to = "urls") %>%
  dplyr::filter(urls != "") %>%
  dplyr::mutate(urls = strsplit(urls, "; ")) %>%
  unnest(urls)
write.csv(invert_docs_subset, "../data/pdf_all_August2024.csv")

invert_docs <- invert_tax %>%
  dplyr::select(TaxonID, class, scientificName,
                SpeciesStatusAssessmentURLs, RecoveryPlanURLs,
                ReviewURLs, CriticalHabitatDesignationURLs) %>%
  tidyr::pivot_longer(ends_with("URLs"),
                      names_to = "documentType", values_to = "urls") %>%
  dplyr::filter(urls != "") %>%
  dplyr::mutate(urls = strsplit(urls, "; ")) %>%
  unnest(urls)

invert_dates <- invert_tax %>%
  dplyr::select(TaxonID, SpeciesStatusAssessmentDates, RecoveryPlanDates,
                ReviewDates, CriticalHabitatDesignationDates) %>%
  tidyr::pivot_longer(ends_with("Dates"),
                      names_to = "documentType", values_to = "dates") %>%
  dplyr::filter(dates != "") %>%
  dplyr::mutate(dates = strsplit(dates, "; ")) %>%
  unnest(dates) %>%
  dplyr::mutate(dates = lubridate::mdy(dates)) %>%
  dplyr::mutate(docYear = lubridate::year(dates)) %>%
  dplyr::group_by(docYear, documentType) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(docYear, documentType, n) %>% unique()

invert_dates_listing <- final_tax %>%
  dplyr::filter(ListingStatus != "Not Listed") %>%
  dplyr::filter(id %!in% chordata_phylum_tax$tax_id) %>%
  dplyr::filter(!is.na(ListingDates), ListingDates != "", ListingDates !=  " ") %>%
  dplyr::mutate(firstListingDate = word(ListingDates, start = 1, end = 1, sep = ";")) %>%
  dplyr::mutate(firstListingDate = lubridate::mdy(firstListingDate)) %>%
  dplyr::mutate(firstListingYear = lubridate::year(firstListingDate)) %>%
  dplyr::mutate(domain = ifelse(scientificName %in% aquatic_spp$scientificName,
                                "Aquatic", "Terrestrial")) %>%
  dplyr::mutate(counter = 1) %>%
  dplyr::group_by(domain) %>%
  dplyr::arrange(firstListingYear) %>%
  dplyr::mutate(cumsum = cumsum(counter)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(domain, firstListingYear) %>%
  dplyr::slice_tail(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(domain, firstListingYear, cumsum) %>% 
  unique()

ggplot()+
  geom_bar(invert_dates_listing,
            mapping = aes(x = firstListingYear, y = cumsum, 
                          color = domain, fill = domain), linewidth = 1,
           stat = "identity")+
  scale_color_manual(name = "Domain",
                     values = c("Aquatic" = "#579c97",
                                "Terrestrial" = "#d98557"))+
  scale_fill_manual(name = "Domain",
                     values = c("Aquatic" = "#579c97",
                                "Terrestrial" = "#d98557"))+
  labs(x = "Initial Listing Year", y = "Cumulative Number of\nListed Species")+
  cowplot::theme_cowplot()+
  cowplot::background_grid()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "top", legend.direction = "horizontal",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave2("../figures/supp/initialListingOverTime.png", dpi = 400,
        height = 4, width = 8)
  
# Filter the data into two groups (vertebrates and invertebrates)
invert_dat <- final_tax %>%
  dplyr::filter(id %!in% chordata_phylum_tax$tax_id) %>%
  dplyr::mutate(higherTax = "Invertebrate") %>%
  dplyr::mutate(hasReview = ifelse(ReviewDates != "" & ReviewDates != " ", TRUE, FALSE),
                hasPlan = ifelse(RecoveryPlanDates != "" & RecoveryPlanDates != " ", TRUE, FALSE),
                hasSSA = ifelse(SpeciesStatusAssessmentDates != "" & SpeciesStatusAssessmentDates != " ", TRUE, FALSE)) %>%
  dplyr::mutate(ListingStatus = ifelse(ListingStatus == "", "Not Listed", ListingStatus)) %>%
  dplyr::mutate(firstListingDate = word(ListingDates, start = 1, end = 1, sep = ";"),
                firstSSADate = word(SpeciesStatusAssessmentDates, start = 1, end = 1, sep = ";"),
                firstRPDate = word(RecoveryPlanDates, start = 1, end = 1, sep = ";"),
                firstRDate = word(ReviewDates, start = 1, end = 1, sep = ";"),
                firstCHDDate = word(CriticalHabitatDesignationDates, start = 1, end = 1, sep = ";")) %>%
  dplyr::mutate(firstListingDateYear = as.numeric(str_extract(firstListingDate, "\\d{4}")),
                firstSSADateYear = as.numeric(str_extract(firstSSADate, "\\d{4}")),
                firstRPDateYear = as.numeric(str_extract(firstRPDate, "\\d{4}")),
                firstRDateYear = as.numeric(str_extract(firstRDate, "\\d{4}")),
                firstCHDDateYear = as.numeric(str_extract(firstCHDDate, "\\d{4}"))) %>%
  dplyr::mutate(firstReviewDatePassed = (2024 - firstListingDateYear) > 5,
                SSADatePassed = ifelse(firstListingDateYear >= 2016 & (2024 - firstListingDateYear), TRUE, FALSE)) %>%
  dplyr::filter(ListingStatus != "Not Listed")
  
vert_dat <- final_tax %>%
  dplyr::filter(id %in% chordata_phylum_tax$tax_id) %>%
  dplyr::mutate(higherTax = "Vertebrate") %>%
  dplyr::mutate(hasReview = ifelse(ReviewDates != "" & ReviewDates != " ", TRUE, FALSE),
                hasPlan = ifelse(RecoveryPlanDates != "" & RecoveryPlanDates != " ", TRUE, FALSE),
                hasSSA = ifelse(SpeciesStatusAssessmentDates != "" & SpeciesStatusAssessmentDates != " ", TRUE, FALSE)) %>%
  dplyr::mutate(ListingStatus = ifelse(ListingStatus == "", "Not Listed", ListingStatus)) %>%
  dplyr::mutate(firstListingDate = word(ListingDates, start = 1, end = 1, sep = ";"),
                firstSSADate = word(SpeciesStatusAssessmentDates, start = 1, end = 1, sep = ";"),
                firstRPDate = word(RecoveryPlanDates, start = 1, end = 1, sep = ";"),
                firstRDate = word(ReviewDates, start = 1, end = 1, sep = ";"),
                firstCHDDate = word(CriticalHabitatDesignationDates, start = 1, end = 1, sep = ";")) %>%
  dplyr::mutate(firstListingDateYear = as.numeric(str_extract(firstListingDate, "\\d{4}")),
                firstSSADateYear = as.numeric(str_extract(firstSSADate, "\\d{4}")),
                firstRPDateYear = as.numeric(str_extract(firstRPDate, "\\d{4}")),
                firstRDateYear = as.numeric(str_extract(firstRDate, "\\d{4}")),
                firstCHDDateYear = as.numeric(str_extract(firstCHDDate, "\\d{4}"))) %>%
  dplyr::mutate(firstReviewDatePassed = (2024 - firstListingDateYear) > 5,
                SSADatePassed = ifelse(firstListingDateYear >= 2016 & (2024 - firstListingDateYear), TRUE, FALSE)) %>%
  dplyr::filter(ListingStatus != "Not Listed")

# Create an alluvial diagram by time
invert_time_alluvial <- invert_dat %>%
  dplyr::select(firstListingDateYear, firstSSADateYear, firstRPDateYear, firstRDateYear, firstCHDDateYear) %>%
  tidyr::pivot_longer(firstListingDateYear:firstCHDDateYear, 
                      names_to = "Document", values_to = "Year") %>%
  dplyr::group_by(Document, Year) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::select(Document, Year, n) %>% unique() %>%
  tidyr::pivot_wider(id_cols = Year, names_from = "Document", values_from = "n") %>%
  dplyr::arrange(Year) %>%
  tidyr::pivot_longer(!Year, names_to = "documentType", values_to = "count") %>%
  dplyr::filter(!is.na(Year), !is.na(count)) %>%
  dplyr::group_by(documentType) %>%
  dplyr::mutate(count_norm = count/sum(count)) %>%
  dplyr::ungroup()
invert_time_alluvial$documentType <- factor(invert_time_alluvial$documentType,
                                            levels = c("firstRDateYear", "firstRPDateYear",
                                                       "firstCHDDateYear", "firstSSADateYear",
                                                       "firstListingDateYear"))

vert_time_alluvial <- vert_dat %>%
  dplyr::select(firstListingDateYear, firstSSADateYear, firstRPDateYear, firstRDateYear, firstCHDDateYear) %>%
  tidyr::pivot_longer(firstListingDateYear:firstCHDDateYear, 
                      names_to = "Document", values_to = "Year") %>%
  dplyr::group_by(Document, Year) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::select(Document, Year, n) %>% unique() %>%
  tidyr::pivot_wider(id_cols = Year, names_from = "Document", values_from = "n") %>%
  dplyr::arrange(Year) %>%
  tidyr::pivot_longer(!Year, names_to = "documentType", values_to = "count") %>%
  dplyr::filter(!is.na(Year), !is.na(count)) %>%
  dplyr::group_by(documentType) %>%
  dplyr::mutate(count_norm = count/sum(count)) %>%
  dplyr::ungroup()
vert_time_alluvial$documentType <- factor(vert_time_alluvial$documentType,
                                            levels = c("firstRDateYear", "firstRPDateYear",
                                                       "firstCHDDateYear", "firstSSADateYear",
                                                       "firstListingDateYear"))

figure1a <- ggplot()+
  # Vertebrates
  geom_bar(vert_time_alluvial,
           mapping = aes(x = Year, y = count * -1, 
                         fill = documentType), stat = "identity", width = 1)+
  # Vertebrate Colorscheme
  scale_fill_manual(name = "Type of Document",
                    values = c("grey90", "grey80", "grey70", "grey60", "grey50"),
                    labels = c("Review",
                               "Reovery Plan",
                               "Critical Habitat Designation", 
                               "Species Status Assessment",
                               "Initial Listing"))+
  ggnewscale::new_scale_fill()+
  # Invertebrates
  geom_bar(invert_time_alluvial,
           mapping = aes(x = Year, y = count, 
                         fill = documentType), stat = "identity", width = 1)+
  # Invertebrate Colorscheme
  MoMAColors::scale_fill_moma_d(name = "Type of Document",
                                palette_name = "OKeeffe",
                                labels = c("Review",
                                           "Reovery Plan",
                                           "Critical Habitat Designation", 
                                           "Species Status Assessment",
                                           "Initial Listing"))+
  # Dividing Lines
  geom_vline(xintercept = 2016, linetype = 2)+
  geom_hline(yintercept = 0)+
  # Axes
  scale_x_continuous(limits = c(1970, 2022), breaks = seq(1970, 2022, by = 5))+
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 25),
                     labels = abs)+
  # Labels
  labs(x = "Year of Initial Documentation", y = "Number of Documents")+
  # Theme
  cowplot::theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "top", legend.direction = "horizontal",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave2("../figures/figure001a_bar.png", figure1a, dpi = 400, height = 4, width = 12)


figure1b <- ggplot()+
  # Vertebrates
  geom_bar(vert_time_alluvial,
           mapping = aes(x = Year, y = count/7378 * -1, 
                         fill = documentType), stat = "identity", width = 1)+
  # Vertebrate Colorscheme
  scale_fill_manual(name = "Type of Document",
                    values = c("grey90", "grey80", "grey70", "grey60", "grey50"),
                    labels = c("Review",
                               "Reovery Plan",
                               "Critical Habitat Designation", 
                               "Species Status Assessment",
                               "Initial Listing"))+
  ggnewscale::new_scale_fill()+
  # Invertebrates
  geom_bar(invert_time_alluvial,
           mapping = aes(x = Year, y = count/30284, 
                         fill = documentType), stat = "identity", width = 1)+
  # Invertebrate Colorscheme
  MoMAColors::scale_fill_moma_d(name = "Type of Document",
                                palette_name = "OKeeffe",
                                labels = c("Review",
                                           "Reovery Plan",
                                           "Critical Habitat Designation", 
                                           "Species Status Assessment",
                                           "Initial Listing"))+
  # Dividing Lines
  geom_vline(xintercept = 2016, linetype = 2)+
  geom_hline(yintercept = 0)+
  # Axes
  scale_x_continuous(limits = c(1970, 2022), breaks = seq(1970, 2022, by = 5))+
  scale_y_continuous(limits = c(-0.025, 0.025), breaks = seq(-0.025, 0.025, by = 0.01),
                     labels = abs)+
  # Labels
  labs(x = "Year of Initial Documentation", y = "Number of Documents/\nSpecies Diversity")+
  # Theme
  cowplot::theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "top", legend.direction = "horizontal",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave2("../figures/figure001b_bar.png", figure1b, dpi = 400, height = 4, width = 12)


figure1c <- ggplot()+
  # Vertebrates
  geom_bar(vert_time_alluvial,
            mapping = aes(x = Year, y = count_norm * -1, 
                          fill = documentType), stat = "identity", width = 1)+
  # Vertebrate Colorscheme
  scale_fill_manual(name = "Type of Document",
                    values = c("grey90", "grey80", "grey70", "grey60", "grey50"),
                    labels = c("Review",
                               "Reovery Plan",
                               "Critical Habitat Designation", 
                               "Species Status Assessment",
                               "Initial Listing"))+
  ggnewscale::new_scale_fill()+
  # Invertebrates
  geom_bar(invert_time_alluvial,
            mapping = aes(x = Year, y = count_norm, 
                          fill = documentType), stat = "identity", width = 1)+
  # Invertebrate Colorscheme
  MoMAColors::scale_fill_moma_d(name = "Type of Document",
                                palette_name = "OKeeffe",
                                labels = c("Review",
                                           "Reovery Plan",
                                           "Critical Habitat Designation", 
                                           "Species Status Assessment",
                                           "Initial Listing"))+
  # Dividing Lines
  geom_vline(xintercept = 2016, linetype = 2)+
  geom_hline(yintercept = 0)+
  # Axes
  scale_x_continuous(limits = c(1970, 2022), breaks = seq(1970, 2022, by = 5))+
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.1),
                     labels = abs)+
  # Labels
  labs(x = "Year of Initial Documentation", y = "Proportion of Total Documents")+
  # Theme
  cowplot::theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "top", legend.direction = "horizontal",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave2("../figures/figure001c_bar.png", figure1c, dpi = 400, height = 4, width = 12)


# Create an alluvial diagram by time just for invertebrate classes
gastropod_class_time_alluvial <- invert_dat %>%
  dplyr::filter(class %in% c("Gastropoda", "Malacostraca", "Bivalvia")) %>%
  dplyr::select(firstListingDateYear, firstSSADateYear, firstRPDateYear, firstRDateYear, firstCHDDateYear) %>%
  tidyr::pivot_longer(firstListingDateYear:firstCHDDateYear, 
                      names_to = "Document", values_to = "Year") %>%
  dplyr::group_by(Document, Year) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::select(Document, Year, n) %>% unique() %>%
  tidyr::pivot_wider(id_cols = Year, names_from = "Document", values_from = "n") %>%
  dplyr::arrange(Year) %>%
  tidyr::pivot_longer(!Year, names_to = "documentType", values_to = "count") %>%
  dplyr::filter(!is.na(Year), !is.na(count)) %>%
  dplyr::group_by(documentType) %>%
  dplyr::mutate(count_norm = count/sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = "Malacology")
gastropod_class_time_alluvial$documentType <- factor(gastropod_class_time_alluvial$documentType,
                                          levels = c("firstRDateYear", "firstRPDateYear",
                                                     "firstCHDDateYear", "firstSSADateYear",
                                                     "firstListingDateYear"))


insect_class_time_alluvial <- invert_dat %>%
  dplyr::filter(class == "Insecta") %>%
  dplyr::select(firstListingDateYear, firstSSADateYear, firstRPDateYear, firstRDateYear, firstCHDDateYear) %>%
  tidyr::pivot_longer(firstListingDateYear:firstCHDDateYear, 
                      names_to = "Document", values_to = "Year") %>%
  dplyr::group_by(Document, Year) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::select(Document, Year, n) %>% unique() %>%
  tidyr::pivot_wider(id_cols = Year, names_from = "Document", values_from = "n") %>%
  dplyr::arrange(Year) %>%
  tidyr::pivot_longer(!Year, names_to = "documentType", values_to = "count") %>%
  dplyr::filter(!is.na(Year), !is.na(count)) %>%
  dplyr::group_by(documentType) %>%
  dplyr::mutate(count_norm = count/sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = "Entomology")
insect_class_time_alluvial$documentType <- factor(insect_class_time_alluvial$documentType,
                                                     levels = c("firstRDateYear", "firstRPDateYear",
                                                                "firstCHDDateYear", "firstSSADateYear",
                                                                "firstListingDateYear"))

other_class_time_alluvial <- invert_dat %>%
  dplyr::filter(class %!in% c("Gastropoda", "Malacostraca", "Bivalvia", "Insecta")) %>%
  dplyr::select(firstListingDateYear, firstSSADateYear, firstRPDateYear, firstRDateYear, firstCHDDateYear) %>%
  tidyr::pivot_longer(firstListingDateYear:firstCHDDateYear, 
                      names_to = "Document", values_to = "Year") %>%
  dplyr::group_by(Document, Year) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::select(Document, Year, n) %>% unique() %>%
  tidyr::pivot_wider(id_cols = Year, names_from = "Document", values_from = "n") %>%
  dplyr::arrange(Year) %>%
  tidyr::pivot_longer(!Year, names_to = "documentType", values_to = "count") %>%
  dplyr::filter(!is.na(Year), !is.na(count)) %>%
  dplyr::group_by(documentType) %>%
  dplyr::mutate(count_norm = count/sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = "Other Invertebrates")
other_class_time_alluvial$documentType <- factor(other_class_time_alluvial$documentType,
                                                  levels = c("firstRDateYear", "firstRPDateYear",
                                                             "firstCHDDateYear", "firstSSADateYear",
                                                             "firstListingDateYear"))

final_class_stream <- rbindlist(list(gastropod_class_time_alluvial,
                                insect_class_time_alluvial,
                                other_class_time_alluvial))

ggplot(final_class_stream)+
  geom_bar(mapping = aes(x = Year, y = count_norm, 
                          fill = documentType),
            size = 0, bw = 0.45)+
  # Colorscheme
  MoMAColors::scale_fill_moma_d(name = "Type of Document",
                                palette_name = "OKeeffe",
                                labels = c("Review",
                                           "Reovery Plan",
                                           "Critical Habitat Designation", 
                                           "Species Status Assessment",
                                           "Initial Listing"))+
  # Dividing Lines
  geom_vline(xintercept = 2016, linetype = 2)+
  geom_hline(yintercept = 0)+
  # Axes
  scale_x_continuous(limits = c(1967, 2023), breaks = seq(1970, 2023, by = 5))+
  scale_y_continuous(limits = c(-00, 0.5), breaks = seq(0, 0.5, by = 0.1),
                     labels = abs)+
  # Labels
  labs(x = "Year of Initial Documentation", y = "Proportion of Total Documents")+
  # Facet Wrap
  facet_wrap(~group, nrow = 3,
             labeller = as_labeller(c(`Entomology` = "Insects",
                           `Malacology` = "Molluscs",
                           `Other Invertebrates` = "Other Invertebrates")))+
  # Theme
  cowplot::theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "top", legend.direction = "horizontal")
ggsave2("../figures/figure002.png", dpi = 400, height = 8, width = 12)

figure1b <- ggplot(final_class_stream)+
  geom_bar(mapping = aes(x = Year, y = count_norm, 
                         fill = documentType), stat = "identity", width = 1)+
  # Colorscheme
  MoMAColors::scale_fill_moma_d(name = "Type of Document",
                                palette_name = "OKeeffe",
                                labels = c("Review",
                                           "Reovery Plan",
                                           "Critical Habitat Designation", 
                                           "Species Status Assessment",
                                           "Initial Listing"))+
  # Dividing Lines
  geom_vline(xintercept = 2016, linetype = 2)+
  geom_hline(yintercept = 0)+
  # Axes
  scale_x_continuous(limits = c(1967, 2023), breaks = seq(1970, 2023, by = 5))+
  scale_y_continuous(limits = c(-00, 0.5), breaks = seq(0, 0.5, by = 0.1),
                     labels = abs)+
  # Labels
  labs(x = "Year of Initial Documentation", y = "Proportion of Total Documents")+
  # Facet Wrap
  facet_wrap(~group, nrow = 3,
             labeller = as_labeller(c(`Entomology` = "Insects",
                                      `Malacology` = "Molluscs",
                                      `Other Invertebrates` = "Other Invertebrates")))+
  # Theme
  cowplot::theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "none", legend.direction = "horizontal",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave2("../figures/figure002_bar.png", figure1b, dpi = 400, height = 8, width = 12)

# Combined Figure
cowplot::plot_grid(figure1a, figure1b,
                   nrow = 2)
ggsave2("../figures/FINAL_FIGURE001.png", dpi = 400, height = 12, width = 12)

###### FIGURE 2 ######
##### NB SCORES ENSEMBLED #####
# Read in the annotation data
my_anns <- data.table::fread("../data/preds/final_anns.csv", header = TRUE)

my_refs <- invert_docs %>%
  dplyr::mutate(documentKey = paste0(TaxonID, "-", urls))

my_anns_ref <- my_anns %>%
  dplyr::mutate(documentKey = paste0(TaxonID, "-", urls))

# Read in the scores
t1_pred <- data.table::fread("../data/preds/t1_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/112)

t2_pred <- data.table::fread("../data/preds/t2_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/51)

t3_pred <- data.table::fread("../data/preds/t3_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/127)

t4_pred <- data.table::fread("../data/preds/t4_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/4)

t5_pred <- data.table::fread("../data/preds/t5_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/53)

t6_pred <- data.table::fread("../data/preds/t6_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/32)

t7_pred <- data.table::fread("../data/preds/t7_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/62)

t8_pred <- data.table::fread("../data/preds/t8_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/367)

t9_pred <- data.table::fread("../data/preds/t9_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/290)

t10_pred <- data.table::fread("../data/preds/t10_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/1)

t11_pred <- data.table::fread("../data/preds/t11_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/112)

t12_pred <- data.table::fread("../data/preds/t12_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/191)

t13_pred <- data.table::fread("../data/preds/t13_preds_all.csv", header = TRUE) %>%
  dplyr::rowwise(V1) %>%
  dplyr::mutate(totalTrue = sum(c_across(where(is.logical)))/32)

pred_df <- data.frame(t1_pred = t1_pred$totalTrue,
                      t2_pred = t2_pred$totalTrue,
                      t3_pred = t3_pred$totalTrue,
                      t4_pred = t4_pred$totalTrue,
                      t5_pred = t5_pred$totalTrue,
                      t6_pred = t6_pred$totalTrue,
                      t7_pred = t7_pred$totalTrue,
                      t8_pred = t8_pred$totalTrue,
                      t9_pred = t9_pred$totalTrue,
                      t11_pred = t11_pred$totalTrue,
                      t12_pred = t12_pred$totalTrue,
                      t13_pred = t13_pred$totalTrue) %>%
  dplyr::mutate(idx = row_number())
pred_df <- pred_df[which(my_anns$`Threats Listed?` == TRUE | is.na(my_anns$`Threats Listed?`)),]

pred_df_full <- cbind(dplyr::filter(my_anns,
                                    `Threats Listed?` == TRUE | is.na(`Threats Listed?`)), 
                      pred_df)
write.csv(pred_df_full, "../data/FINAL_THREATS.csv")
final_spp_domains <- data.table::fread("../data/FINAL_THREATS_SPP.csv")

# Split into terrestrial and aquatic
gbif_occ <- data.table::fread("../data/gbif.txt", header = TRUE, stringsAsFactors = FALSE,
                              quote = "") %>%
  dplyr::filter(!is.na(decimalLongitude)) %>%
  dplyr::select(order, species, decimalLongitude, decimalLatitude) %>% unique() %>%
  sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"))
sf::st_crs(gbif_occ) <- "wgs84"

#### AQUATIC THREATS ####
aquatic_spp <- final_spp_domains %>%
  dplyr::filter(domain == "Aquatic")

aquatic_pred_df <- pred_df_full %>%
  dplyr::filter(scientificName %in% aquatic_spp$scientificName)

pred_df_aquatic_long <- aquatic_pred_df %>%
  dplyr::group_by(scientificName) %>%
  dplyr::mutate(housing_spMean = mean(t1_pred),
                agricult_spMean = mean(t2_pred),
                mining_spMean = mean(t3_pred),
                bioresuse_spMean = mean(t5_pred),
                humanintrus_spMean = mean(t6_pred),
                natsysmod_spMean = mean(t7_pred),
                problembio_spMean = mean(t8_pred),
                pollution_spMean = mean(t9_pred),
                climatechange_spMean = mean(t11_pred),
                intrinsic_spMean = mean(t12_pred)) %>%
  dplyr::ungroup() %>%
  dplyr::select(scientificName, housing_spMean:intrinsic_spMean) %>% unique() %>%
  pivot_longer(-scientificName, names_to = "cat", values_to = "prev")

pred_df_aquatic_summ <- pred_df_aquatic_long %>%
  dplyr::group_by(cat) %>%
  dplyr::mutate(meanPrev = median(prev),
                lowerPrev = quantile(prev, probs = 0.25),
                upperPrev = quantile(prev, probs = 0.75)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cat, meanPrev, lowerPrev, upperPrev) %>% unique() %>%
  dplyr::arrange(meanPrev) %>%
  dplyr::mutate(ordering = row_number())

pred_df_aquatic_cat_fact <- factor(pred_df_aquatic_summ$cat, 
                                   levels = pred_df_aquatic_summ$cat)

Figure2_a <- ggplot()+
  geom_jitter(pred_df_aquatic_long,
              mapping = aes(x = cat, y = prev, color = cat),
              alpha = 0.1)+
  geom_pointrange(pred_df_aquatic_summ,
                  mapping = aes(x = pred_df_aquatic_cat_fact, y = meanPrev, ymin = lowerPrev,
                                ymax = upperPrev, color = pred_df_aquatic_cat_fact), size = 2.5,
                  linewidth = 1)+
  geom_text(pred_df_aquatic_summ,
            mapping = aes(x = pred_df_aquatic_cat_fact, y = meanPrev,
                          label = round(meanPrev, 2)), color = "white",
            fontface = "bold")+
  colorspace::scale_color_discrete_sequential(palette = "DarkMint")+
  scale_x_discrete(labels = c("Bio. Res. Use", "Human Intrusion",
                              "Energy Prod. & Mining", "Agriculture",
                              "Res./Comm. Development", "Climate Change",
                              "Problematic Biotic Factors", "Intrinsic Factors",
                              "Natural Sys. Mod.", "Pollution"), 
                   name = "Threat Category")+
  scale_y_continuous(limits = c(-0.05,1.05), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent_format(),
                     name = "Predicted Prevelance\nAcross Species")+
  labs(title = "(a) Aquatic Invertebrates")+
  cowplot::theme_cowplot()+
  cowplot::background_grid()+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 13),
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "none")
cowplot::ggsave2("../figures/threatPrevalence_aquatic.png", Figure2_a,
                 dpi = 400, height = 5, width = 10)

#### TERRESTRIAL THREATS ####
terrestrial_spp <- final_spp_domains %>%
  dplyr::filter(domain == "Terrestrial")

terrestrial_pred_df <- pred_df_full %>%
  dplyr::filter(scientificName %in% terrestrial_spp$scientificName)

pred_df_terrestrial_long <- terrestrial_pred_df %>%
  dplyr::group_by(scientificName) %>%
  dplyr::mutate(housing_spMean = mean(t1_pred),
                agricult_spMean = mean(t2_pred),
                mining_spMean = mean(t3_pred),
                bioresuse_spMean = mean(t5_pred),
                humanintrus_spMean = mean(t6_pred),
                natsysmod_spMean = mean(t7_pred),
                problembio_spMean = mean(t8_pred),
                pollution_spMean = mean(t9_pred),
                climatechange_spMean = mean(t11_pred),
                intrinsic_spMean = mean(t12_pred)) %>%
  dplyr::ungroup() %>%
  dplyr::select(scientificName, housing_spMean:intrinsic_spMean) %>% unique() %>%
  pivot_longer(-scientificName, names_to = "cat", values_to = "prev")

pred_df_terrestrial_summ <- pred_df_terrestrial_long %>%
  dplyr::group_by(cat) %>%
  dplyr::mutate(meanPrev = median(prev),
                lowerPrev = quantile(prev, probs = 0.25),
                upperPrev = quantile(prev, probs = 0.75)) %>%
  dplyr::ungroup() %>%
  dplyr::select(cat, meanPrev, lowerPrev, upperPrev) %>% unique() %>%
  dplyr::arrange(meanPrev) %>%
  dplyr::mutate(ordering = row_number())

pred_df_terrestrial_cat_fact <- factor(pred_df_terrestrial_summ$cat, 
                                   levels = pred_df_terrestrial_summ$cat)

Figure2_b <- ggplot()+
  geom_jitter(pred_df_terrestrial_long,
              mapping = aes(x = cat, y = prev, color = cat),
              alpha = 0.1)+
  geom_pointrange(pred_df_terrestrial_summ,
                  mapping = aes(x = pred_df_terrestrial_cat_fact, y = meanPrev, ymin = lowerPrev,
                                ymax = upperPrev, color = pred_df_terrestrial_cat_fact), size = 2.5,
                  linewidth = 1)+
  geom_text(pred_df_terrestrial_summ,
            mapping = aes(x = pred_df_terrestrial_cat_fact, y = meanPrev,
                          label = round(meanPrev, 2)), color = "white",
            fontface = "bold")+
  colorspace::scale_color_discrete_sequential(palette = "Heat")+
  scale_x_discrete(labels = c("Energy Prod. & Mining", "Bio. Res. Use", 
                              "Agriculture", "Natural Sys. Mod.", "Human Intrusion",
                              "Intrinsic Factors", "Pollution", "Res./Comm. Development",
                              "Climate Change", "Problematic Biotic Factors"), 
                   name = "Threat Category")+
  scale_y_continuous(limits = c(-0.05,1.05), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent_format(),
                     name = "Predicted Prevelance\nAcross Species")+
  labs(title = "(b) Terrestrial Invertebrates")+
  cowplot::theme_cowplot()+
  cowplot::background_grid()+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 13),
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "none")
cowplot::ggsave2("../figures/threatPrevalence_terrestrial.png", Figure2_b,
                 dpi = 400, height = 5, width = 10)

##### FIGURE 2 #####
Figure2 <- cowplot::plot_grid(Figure2_a, Figure2_b, nrow = 2)
cowplot::ggsave2("../figures/Figure2.png", Figure2,
                 dpi = 400, height = 8, width = 10)