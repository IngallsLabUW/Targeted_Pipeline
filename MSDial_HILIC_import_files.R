source("functions.R")

# Import all MSDial files --------------------------------------------------

filenames <- RemoveCsv(list.files(path = 'data_raw', pattern = '*.csv'))

for (i in filenames) {
  filepath <- file.path('data_raw', paste(i,".csv", sep = ""))
  assign(i, read.csv(filepath, stringsAsFactors = FALSE))
}


columns.to.drop <- c('Average.Rt.min.', 'Formula', 'Ontology', 'INCHIKEY', 'SMILES', 'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
'MS1.isotopic.spectrum', 'MS.MS.spectrum', 'Average.Mz', 'Post.curation.result', 'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
'm.z.matched', 'MS.MS.matched', 'Manually.modified', 'Total.score', 'RT.similarity', 'Dot.product', 'Reverse.dot.product', 'Fragment.presence..')

# Set header, filter unknowns ---------------------------------------
# Known issues: choosing which columns remain by hand.
# Also manually changing classes by choosing "10" rather than a filter.

runs <- grep("hilic", names(.GlobalEnv), value = TRUE, ignore.case = TRUE)
runlist <- do.call("list", mget(runs))

headers.set <- lapply(names(runlist), function(x) SetHeader(runlist[[x]]))
names(headers.set) <- runs

for (df in seq_along(headers.set)) { 
  headers.set[[df]] <- headers.set[[df]] %>% filter(!Metabolite.name == "Unknown")
  headers.set[[df]] <- headers.set[[df]] %>% select(-one_of(columns.to.drop))
  headers.set[[df]] <- headers.set[[df]] %>% rename(Metabolite.Name = Metabolite.name)
}

# Change variable classes -------------------------------------------------
classes.changed <- lapply(names(headers.set), function(x) ChangeClasses(headers.set[[x]]))
names(classes.changed) <- runs

list2env(classes.changed, globalenv())


# Change dataset orientation from wide to long -------------------------------------------------

# Positive
Area.pos <- RearrangeDatasets(Area_HILICPos_Example, parameter = "Area.Value")
Mz.pos   <- RearrangeDatasets(Mz_HILICPos_Example, parameter = "Mz.Value")
RT.pos   <- RearrangeDatasets(RT_HILICPos_Example, parameter = "RT.Value")
SN.pos   <- RearrangeDatasets(SN_HILICPos_Example, parameter = "SN.Value")


# Negative
Area.neg <- RearrangeDatasets(Area_HILICNeg_Example, parameter = "Area.Value")
Mz.neg   <- RearrangeDatasets(Mz_HILICNeg_Example, parameter = "Mz.Value")
RT.neg   <- RearrangeDatasets(RT_HILICNeg_Example, parameter = "RT.Value")
SN.neg   <- RearrangeDatasets(SN_HILICNeg_Example, parameter = "SN.Value")



# Combine to one dataset --------------------------------------------------
combined.pos <- Area.pos %>%
  left_join(Mz.pos) %>%
  left_join(SN.pos) %>%
  left_join(RT.pos) %>%
  mutate(Column = "HILICPos") %>%
  select(Replicate.Name, Column, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

combined.neg <- Area.neg %>%
  left_join(Mz.neg) %>%
  left_join(SN.neg) %>%
  left_join(RT.neg) %>%
  mutate(Column = "HILICNeg") %>%
  select(Replicate.Name, Column, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

# Standardize dataset --------------------------------------------------

combined.final <- combined.neg %>%
  bind_rows(combined.pos)

combined.final <- StandardizeMetabolites(combined.final) 

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/MSDial_combined_", currentDate, ".csv", sep = "")

write.csv(combined.final, csvFileName, row.names = FALSE)

rm(list = ls())