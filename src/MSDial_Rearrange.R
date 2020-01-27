# Set header, filter unknowns ---------------------------------------
columns.to.drop <- c('Average.Rt.min.', 'Formula', 'Ontology', 'INCHIKEY', 'SMILES', 
                     'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
                      'MS1.isotopic.spectrum', 'MS.MS.spectrum', 'Average.Mz', 'Post.curation.result', 
                      'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
                      'm.z.matched', 'MS.MS.matched', 'Manually.modified', 'Total.score', 'RT.similarity', 
                     'Dot.product', 'Reverse.dot.product', 'Fragment.presence..')

runs <- grep(matching.pattern, names(.GlobalEnv), value = TRUE, ignore.case = TRUE)
runlist <- do.call("list", mget(runs))

headers.set <- lapply(names(runlist), function(x) SetHeader(runlist[[x]]))
names(headers.set) <- runs

for (df in seq_along(headers.set)) {
  headers.set[[df]] <- headers.set[[df]] %>% filter(!Metabolite.name == "Unknown")
  headers.set[[df]] <- headers.set[[df]] %>% select(-one_of(columns.to.drop))
  headers.set[[df]] <- headers.set[[df]] %>% rename(Metabolite.Name = Metabolite.name)
}

# Change variable classes ---------------------------------------------------------------------
classes.changed <- lapply(names(headers.set), function(x) ChangeClasses(headers.set[[x]]))
names(classes.changed) <- runs

list2env(classes.changed, globalenv())


# Rearrange data and combine to one dataframe -------------------------------------------------
if (TRUE %in% grepl("positive|negative", names(.GlobalEnv), ignore.case = TRUE)) {

  # HILIC Positive
  Area.positive <- RearrangeDatasets(Area.positive, parameter = "Area.Value")
  Mz.positive   <- RearrangeDatasets(Mz.positive, parameter = "Mz.Value")
  RT.positive   <- RearrangeDatasets(RT.positive, parameter = "RT.Value")
  SN.positive   <- RearrangeDatasets(SN.positive, parameter = "SN.Value")

  # HILIC Negative
  Area.negative <- RearrangeDatasets(Area.negative, parameter = "Area.Value")
  Mz.negative   <- RearrangeDatasets(Mz.negative, parameter = "Mz.Value")
  RT.negative   <- RearrangeDatasets(RT.negative, parameter = "RT.Value")
  SN.negative   <- RearrangeDatasets(SN.negative, parameter = "SN.Value")


  # Combine to one dataset
  combined.pos <- Area.positive %>%
    left_join(Mz.positive) %>%
    left_join(SN.positive) %>%
    left_join(RT.positive) %>%
    mutate(Column = "HILICPos") %>%
    select(Replicate.Name, Column, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

  combined.neg <- Area.negative %>%
    left_join(Mz.negative) %>%
    left_join(SN.negative) %>%
    left_join(RT.negative) %>%
    mutate(Column = "HILICNeg") %>%
    select(Replicate.Name, Column, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

  combined.final <- combined.neg %>%
    bind_rows(combined.pos)

  } else {

  # Cyano
  Area <- RearrangeDatasets(Area.RP.Cyano, parameter = "Area.Value")
  Mz   <- RearrangeDatasets(Mz.RP.Cyano, parameter = "Mz.Value")
  RT   <- RearrangeDatasets(RT.RP.Cyano, parameter = "RT.Value")
  SN   <- RearrangeDatasets(SN.RP.Cyano, parameter = "SN.Value")

  # Combine to one dataset
  combined.final <- Area %>%
    left_join(Mz) %>%
    left_join(SN) %>%
    left_join(RT) %>%
    select(Replicate.Name, Area.Value, Mz.Value, RT.Value, SN.Value, everything())
}


# Standardize dataset --------------------------------------------------
combined.final <- StandardizeMetabolites(combined.final)

currentDate <- Sys.Date()
csvFileName <- paste("data_intermediate/MSDial_combined_", file.pattern, "_", currentDate, ".csv", sep = "")

write.csv(combined.final, csvFileName, row.names = FALSE)

rm(list = setdiff(ls()[!ls() %in% c("file.pattern")], lsf.str()))