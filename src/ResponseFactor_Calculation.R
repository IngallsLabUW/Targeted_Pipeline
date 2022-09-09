# Handle on-the-fly response factor calculation, and stop the program for those compounds that do not have response factors.


# Get response factors -----------------------------------------------------
Full.data.RF <- Full.stds.data %>%
  mutate(RF = Area.with.QC/Concentration_uM) %>%
  filter(!Compound_Type == "Internal Standard") %>%
  mutate(Replicate.Name = substr(Replicate.Name, 1, nchar(Replicate.Name)-2))

# In HILIC compounds, filter mixes.
if ("Column" %in% colnames(Full.data.RF)) {
  Full.data.RF <- Full.data.RF %>%
    filter(str_detect(Replicate.Name, as.character(HILIC_Mix)) | str_detect(Replicate.Name, "H2OInMatrix"))
}

# Calculate RF max and min using only standards in water.
Full.data.RF.dimensions <- Full.data.RF %>%
  filter(Type == "Standards_Water") %>%
  group_by(Metabolite.Name) %>%
  mutate(RF.max = max(RF, na.rm = TRUE),
         RF.min = min(RF, na.rm = TRUE))

Full.data.RF.dimensions$RF.max[is.infinite(Full.data.RF.dimensions$RF.max) | is.nan(Full.data.RF.dimensions$RF.max) ] <- NA
Full.data.RF.dimensions$RF.min[is.infinite(Full.data.RF.dimensions$RF.min) | is.nan(Full.data.RF.dimensions$RF.min) ] <- NA

Full.data.RF.dimensions <- Full.data.RF.dimensions %>%
  mutate(RF.diff = RF.max/RF.min) %>%
  unique()

print(paste("NAs or NaNs in the calculated response factors:", TRUE %in% is.na(Full.data.RF.dimensions)))
no.rfs <- Full.data.RF.dimensions[is.na(Full.data.RF.dimensions$RF),]
print(unique(no.rfs$Metabolite.Name))

# Calculate response factor ratios ----------------------------------------
# Calculate the response factor ratios using (Standards in Matrix - Water in Matrix) / (Standards in water) for each replicate.
temp.RF.ratios <- Full.data.RF %>%
  group_by(Metabolite.Name, Type) %>%
  mutate(RF.mean.per_sampleID = mean(RF, na.rm = TRUE)) %>%
  select(Replicate.Name, Metabolite.Name, Type, RF.mean.per_sampleID) %>%
  unique()

print(paste("NAs or NaNs in the calculated response factor ratios:", TRUE %in% is.na(temp.RF.ratios)))
metabolite.issues <- temp.RF.ratios[is.nan(temp.RF.ratios$RF.mean.per_sampleID),]
print(unique(metabolite.issues$Metabolite.Name))


Full.data.RF.ratios <- temp.RF.ratios %>%
  filter(!(is.infinite(RF.mean.per_sampleID) | is.nan(RF.mean.per_sampleID))) %>%
  group_by(Metabolite.Name) %>% filter(n() >= 3) %>%
  mutate(RF.ratio = ((RF.mean.per_sampleID[Type == "Standards_Matrix"] - RF.mean.per_sampleID[Type == "Water_Matrix"]) 
                     / RF.mean.per_sampleID[Type == "Standards_Water"])) %>%
  select(Metabolite.Name, RF.ratio) %>%
  unique()
