setwd("C:/Users/P70094486/OneDrive - Maastricht University/Desktop/rare epilepsy/EAGs/solid result")
getwd()


library(readxl)
library(openxlsx)
library(tidyr)
library(dplyr)

wp <- read.xlsx("wp_edge_v2.xlsx")

## delete the prefix in the columns

wp$sources <- gsub("https://identifiers.org/", "", wp$sources)
wp$targets <- gsub("https://identifiers.org/", "", wp$targets)
### Separate the column to source and identifier
wp <- separate(wp, "targets", into = c("targetdb", "targetid"), sep = "/", fill = "right")
wp$interaction <- sub(".*Pathway/([^_]+)_.*", "\\1", wp$interaction)
wp <- wp %>%
  mutate(database = "WikiPathways") %>%
  rename(pathwayid = interaction)

### unifying the lables by identifiers
wp_cleaned <- wp %>%
  group_by(sourceid) %>%
  mutate(sourceLabel = first(sourceLabel)) %>%
  ungroup()

wp_cleaned <- wp %>%
  group_by(targetid) %>%
  mutate(targetLabel = first(targetLabel)) %>%
  ungroup()

wp_cleaned <- wp_cleaned[!duplicated(wp_cleaned), ]

### How many bio molecules in total
length(unique(wp_cleaned$targetid, wp_cleaned$sourceid))

### Load the normal metabolites 
removemetabolite <- read.xlsx("metabolites to delete.xlsx")
to_remove <- tolower(c(removemetabolite$wikiid, removemetabolite$chemicalLabel, removemetabolite$chembl_id, removemetabolite$pubchem_id))

wp_cleaned_v2 <- wp_cleaned %>%
  filter(
    !tolower(sourceid) %in% to_remove,
    !tolower(sourceLabel) %in% to_remove,
    !tolower(targetid) %in% to_remove,
    !tolower(targetLabel) %in% to_remove,
  )

### Load the string and intact
intact <- read.xlsx("intact_edge.xlsx")
string <- read.xlsx("string_edges.xlsx")
wp_cleaned_v2 <- wp_cleaned_v2 %>%
  rename( "pathwayid/interactionweight" = "pathwayid")



intact <- mutate(intact, across(everything(), as.character))
string <- mutate(string, across(everything(), as.character))
wp_cleaned_v2 <- mutate(wp_cleaned_v2, across(everything(), as.character))
merged <- bind_rows(intact, string, wp_cleaned_v2)

### Put the weight from sources
weight <- read.xlsx("genename_source.xlsx")
counts_df <- weight %>%
  group_by(genename) %>%
  summarise(count = n())
length(unique(merged$targetLabel, merged$sourceLabel))
write.xlsx(merged, file = "merged.xlsx")
write.xlsx(counts_df, file = "nodeweight.xlsx")
### Laod the pathway in Cytoscape
library(RCy3)
library(igraph)
cytoscapePing()

colnames(merged)

browseVignettes("RCy3")
### how to add the heat in the end 
