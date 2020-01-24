# This code retrieves mol/L from peak areas of targeted compounds.

# Get response factors for transect compounds ----------------------------------
Full.data.RF <- Full.data %>%
  mutate(RF = Area.with.QC/Conc..uM) %>%
  filter(!Compound.Type == "Internal Standard") %>%
  mutate(Replicate.Name = substr(Replicate.Name, 1, nchar(Replicate.Name)-2))

# In HILIC compounds, filter mixes.
if ("Column" %in% colnames(Full.data.RF)) {
  Full.data.RF <- Full.data.RF %>%
    filter(str_detect(Replicate.Name, as.character(HILICMix)) | str_detect(Replicate.Name, "H2OInMatrix"))
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

# If applicable, supplement data with information from calculated Ingalls QE.RF ratios.
missing.RF <- setdiff(unique(temp.RF.ratios$Metabolite.Name), Full.data.RF.ratios$Metabolite.Name)

test.standards <- Ingalls.Standards %>%
  filter(Metabolite.Name %in% missing.RF) %>%
  rename(RF.ratio = QE.RF.ratio) %>%
  select(Metabolite.Name, RF.ratio)

Full.data.RF.ratios <- Full.data.RF.ratios %>%
  as.data.frame() %>%
  rbind(test.standards) %>%
  filter(!is.na(RF.ratio)) %>%
  mutate(RF.ratio = as.numeric(RF.ratio))

currentDate = Sys.Date()
write.csv(Full.data.RF.ratios, paste("data_intermediate/MSDIAL_ResponseFactorRatios_", currentDate, ".csv", sep = ""))

# Quantify samples for the BMIS'd dataset ---------------------------------
BMISd.data.filtered <- BMISd.data %>%
  separate(Run.Cmpd, c("Sample.Name"), extra = "drop", fill = "right") %>%
  mutate(Metabolite.Name = Mass.Feature) %>%
  filter(Metabolite.Name %in% Full.data.RF.ratios$Metabolite.Name) %>%
  left_join(Full.data.RF.ratios) %>%
  left_join(Full.data.RF.dimensions %>% select(Metabolite.Name, RF.max, RF.min) %>% unique(), by = "Metabolite.Name") %>%
  select(Metabolite.Name, FinalBMIS, Sample.Name, Adjusted.Area, everything())


# Calculate umol/vial for compounds without an internal standard ----------
Quantitative.data <- BMISd.data.filtered %>%
  mutate(RF.ave = as.numeric(rowMeans(BMISd.data.filtered[, c("RF.min", "RF.max")]))) %>%
  mutate(umol.in.vial.ave = Adjusted.Area/RF.ave/RF.ratio,
         umol.in.vial.max = Adjusted.Area/RF.min/RF.ratio,
         umol.in.vial.min = Adjusted.Area/RF.max/RF.ratio) %>%
  select(Metabolite.Name:Adjusted.Area, everything())

# Pull out data for matched internal standards ----------------------------
IS.key <- BMISd.data.filtered %>%
  select(FinalBMIS, Metabolite.Name) %>%
  unique() %>%
  left_join(original.IS.key %>% select(FinalBMIS, Concentration_nM)) %>%
  filter(str_detect(FinalBMIS, Metabolite.Name))


# Calculate umol/vial for compounds with matched internal standards -----------------
IS.data <- Full.data %>%
  filter(Metabolite.Name %in% IS.key$FinalBMIS) %>%
  mutate(IS_Area = Area.with.QC,
         FinalBMIS = Metabolite.Name) %>%
  select(IS_Area, FinalBMIS, Replicate.Name) %>%
  left_join(IS.key %>% select(FinalBMIS, Metabolite.Name, Concentration_nM))

matched.IS.compounds <- data.frame(Compounds = c(IS.key[ ,"FinalBMIS"], as.character(IS.key[ ,"Metabolite.Name"])))

IS.sample.data <- QCd.data %>%
  left_join(IS.data %>% select(FinalBMIS, Metabolite.Name, Concentration_nM), by = "Metabolite.Name") %>%
  unique() %>%
  filter(Metabolite.Name %in% matched.IS.compounds$Compounds) %>%
  filter(!str_detect(Replicate.Name, "Std")) %>%
  mutate(Std.Type = ifelse(str_detect(Metabolite.Name, ","), "Internal_std", "Standard")) %>%
  mutate(testcol1 = ifelse(str_detect(Metabolite.Name, ","), sapply(strsplit(Metabolite.Name, ","), `[`, 1), Metabolite.Name)) %>%
  mutate(Names = ifelse(str_detect(testcol1, "-"), sapply(strsplit(testcol1, "-"), `[`, 2), testcol1)) %>%
  mutate(Pairs = ifelse(!str_detect(Metabolite.Name, ","), Metabolite.Name, paste(Names, "IS", sep = "_"))) %>%
  select(-c("Pairs", "testcol1", "Run.Type")) %>%
  arrange(Replicate.Name) %>%
  group_by(Names) %>%
  group_split()

IS.mid_frame <- lapply(IS.sample.data, function(x) group_by(x, Replicate.Name))

IS.mid_frame2 <- lapply(IS.mid_frame,
                        function(x) mutate(x,
                        umol.in.vial_IS = (Area.with.QC[Std.Type == "Standard"] / Area.with.QC[Std.Type == "Internal_std"]) * (Concentration_nM[Std.Type == "Standard"]/1000)))

IS.sample.data <- do.call(rbind, IS.mid_frame2) %>%
  filter(!str_detect(Metabolite.Name, ",")) %>%
  rename(Sample.Name = Replicate.Name) %>%
  select(Sample.Name:Area.with.QC, Concentration_nM, umol.in.vial_IS)

rm(list = c("matched.IS.compounds", "QCd.data", "IS.mid_frame", "IS.mid_frame2"))


# Add matched IS_smp info back into main frame ------------------------------------------------
All.Info <- Quantitative.data %>%
  select(Metabolite.Name, runDate:replicate, Adjusted.Area, Area.with.QC, RF.ratio:umol.in.vial.min) %>%
  unite(Sample.Name, c("runDate", "type", "SampID", "replicate"), remove = FALSE) %>%
  left_join(IS.sample.data %>% select(Sample.Name, Metabolite.Name, umol.in.vial_IS)) %>%
  mutate(umol.in.vial.ave = ifelse(is.na(umol.in.vial_IS), umol.in.vial.ave, umol.in.vial_IS),
         umol.in.vial.max = ifelse(is.na(umol.in.vial_IS), umol.in.vial.max, NA),
         umol.in.vial.min = ifelse(is.na(umol.in.vial_IS), umol.in.vial.min, NA)) %>%
  rename(Replicate.Name = Sample.Name) %>%
  select(-runDate, -type, -replicate) # %>%
  #filter(!str_detect(Replicate.Name, "DDA"))

# Add in dilution factor and filtered volume --------------------------------------------------
All.Info.Quantitative <- All.Info %>%
  mutate(nmol.in.Enviro.ave = (umol.in.vial.ave*10^-6*Injection.Volume/Volume.Filtered*1000*Dilution.Factor)) %>%
  left_join(Full.data %>% select(Metabolite.Name, Emperical.Formula)) %>%
  unique()

# Get molecules of carbon and nitrogen ------------------------------------
All.Info.Molecules <- All.Info.Quantitative  %>%
  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
                    str_extract(Emperical.Formula, "^C\\d"),
                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Emperical.Formula, "N\\D"),
                    1,
                    str_extract(Emperical.Formula, "N\\d"))) %>%
  mutate(N = as.numeric(str_replace_all(N, "N", ""))) %>%
  mutate(nmol.C.ave = nmol.in.Enviro.ave*C,
         nmol.N.ave = nmol.in.Enviro.ave*N ) %>%
  select(Metabolite.Name, SampID, Replicate.Name, everything())

# Summarize for each metabolite ------------------------------------
All.Info.Summed <- All.Info.Molecules %>%
  group_by(Metabolite.Name) %>%
  summarise(nmol.Enviro.med = median(nmol.in.Enviro.ave, na.rm  = T),
            nmol.Enviro.min = min(nmol.in.Enviro.ave, na.rm  = T),
            nmol.Enviro.max = max(nmol.in.Enviro.ave, na.rm  = T),
            nmol.C.med = median(nmol.C.ave, na.rm  = T),
            nmol.C.min = min(nmol.C.ave, na.rm  = T),
            nmol.C.max = max(nmol.C.ave, na.rm  = T)) %>%
  arrange(desc(nmol.Enviro.med))

# Summarize total carbon and nitrogen for each compound ------------------------------------
Final.All.perSampID <- All.Info.Molecules %>%
  select(SampID, nmol.C.ave, nmol.N.ave) %>%
  group_by(SampID) %>%
  summarise(totalCmeasured_nM_perID = sum(as.numeric(nmol.C.ave), na.rm = TRUE),
            totalNmeasured_nM_perID = sum(as.numeric(nmol.N.ave), na.rm = TRUE))


# Calculate mole fractions of each compound ------------------------------------
Final.Quantitative <- All.Info.Molecules %>%
  unique() %>%
  left_join(Final.All.perSampID) %>%
  mutate(ratioCN = totalCmeasured_nM_perID / totalNmeasured_nM_perID) %>%
  mutate(molFractionC = nmol.C.ave/totalCmeasured_nM_perID,
         molFractionN = nmol.N.ave/totalNmeasured_nM_perID) %>%
  select(Metabolite.Name, Replicate.Name, Adjusted.Area, Area.with.QC, RF.ratio:molFractionN) %>%
  unique()

Final.Quantitative.Summed <- Final.Quantitative %>%
  group_by(Metabolite.Name) %>%
  summarise(nmol.Enviro.med = median(nmol.in.Enviro.ave, na.rm  = T),
            nmol.Enviro.min = min(nmol.in.Enviro.ave, na.rm  = T),
            nmol.Enviro.max = max(nmol.in.Enviro.ave, na.rm  = T),
            nmol.C.med = median(nmol.C.ave, na.rm  = T),
            nmol.C.min = min(nmol.C.ave, na.rm  = T),
            nmol.C.max = max(nmol.C.ave, na.rm  = T),
            mol.C.Fraction.med = median(molFractionC, na.rm = T),
            mol.C.Fraction.min = min(molFractionC, na.rm = T),
            mol.C.Fraction.max = max(molFractionC, na.rm = T)) %>%
  arrange(desc(Metabolite.Name))
