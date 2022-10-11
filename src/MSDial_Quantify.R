
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
  select(Metabolite.Name:Adjusted.Area, everything()) %>%
  unique()

# Pull out data for matched internal standards ----------------------------
IS.key <- BMISd.data.filtered %>%
  select(FinalBMIS, Metabolite.Name) %>%
  unique() %>%
  left_join(Ingalls.Standards %>% select(Metabolite.Name, Concentration_uM) %>% rename(FinalBMIS = Metabolite.Name)) %>%
  filter(str_detect(FinalBMIS, Metabolite.Name))

# Calculate umol/vial for compounds with matched internal standards -----------------
IS.data <- Full.stds.data %>%
  filter(Metabolite.Name %in% IS.key$FinalBMIS) %>%
  mutate(IS_Area = Area.with.QC,
         FinalBMIS = Metabolite.Name) %>%
  select(IS_Area, FinalBMIS, Replicate.Name) %>%
  left_join(IS.key %>% select(FinalBMIS, Metabolite.Name, Concentration_uM))

matched.IS.compounds <- data.frame(Compounds = c(IS.key[ ,"FinalBMIS"], as.character(IS.key[ ,"Metabolite.Name"])))

IS.sample.data <- QCd.data %>%
  left_join(IS.data %>% select(FinalBMIS, Metabolite.Name, Concentration_uM), by = "Metabolite.Name") %>%
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
                        umol.in.vial_IS = (Area.with.QC[Std.Type == "Standard"] / Area.with.QC[Std.Type == "Internal_std"]) * (Concentration_uM[Std.Type == "Standard"]/1000)))

IS.sample.data <- do.call(rbind, IS.mid_frame2) %>%
  filter(!str_detect(Metabolite.Name, ",")) %>%
  rename(Sample.Name = Replicate.Name) %>%
  select(Sample.Name:Area.with.QC, Concentration_uM, umol.in.vial_IS)

rm(list = c("matched.IS.compounds", "IS.mid_frame", "IS.mid_frame2"))


# Add matched IS_smp info back into main frame ------------------------------------------------
All.Info <- Quantitative.data %>%
  select(Metabolite.Name, runDate:replicate, Adjusted.Area, Area.with.QC, RF.ratio:umol.in.vial.min) %>%
  unite(Sample.Name, c("runDate", "type", "SampID", "replicate"), remove = FALSE) %>%
  left_join(IS.sample.data %>% select(Sample.Name, Metabolite.Name, umol.in.vial_IS)) %>%
  mutate(umol.in.vial.ave = ifelse(is.na(umol.in.vial_IS), umol.in.vial.ave, umol.in.vial_IS),
         umol.in.vial.max = ifelse(is.na(umol.in.vial_IS), umol.in.vial.max, NA),
         umol.in.vial.min = ifelse(is.na(umol.in.vial_IS), umol.in.vial.min, NA)) %>%
  rename(Replicate.Name = Sample.Name) %>%
  select(-runDate, -type, -replicate) 

# Add in dilution factor and filtered volume --------------------------------------------------
All.Info.Quantitative <- All.Info %>%
  mutate(nmol.in.Enviro.ave = (umol.in.vial.ave*10^-6*Reconstitution.Volume/Volume.Filtered*1000*Dilution.Factor)) %>%
  left_join(Full.stds.data %>% select(Metabolite.Name, Empirical_Formula)) %>%
  unique()

# Get molecules of carbon and nitrogen ------------------------------------
All.Info.Molecules <- All.Info.Quantitative  %>%
  mutate(C = ifelse(is.na(str_extract(Empirical_Formula, "^C\\d\\d")),
                    str_extract(Empirical_Formula, "^C\\d"),
                    str_extract(Empirical_Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Empirical_Formula, "N\\D"),
                    1,
                    str_extract(Empirical_Formula, "N\\d"))) %>%
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
