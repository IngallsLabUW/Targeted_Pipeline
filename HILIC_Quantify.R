source("Functions.R")

# This code retrieves mol/L from peak areas of targeted compounds.
# This run is for the MESOSCOPE Eddy Transect HILIC run.

# Volume of seawater filtered for all TRANSECT samples was 5 L.

# User data ---------------------------------------------------------------
standards.pattern = "Ingalls"
BMIS.pattern = "BMIS"
QC.pattern = "QC"
names.pattern = "Names"
Column.Type = "HILIC"


Dilution.Factor = 2
Injection.Volume = 400 # nanomoles
Volume.Filtered = 10 # liters


# Import standards and filter NAs ---------------------------------------------------------------
filename <- RemoveCsv(list.files(path = 'data_extras/', pattern = standards.pattern))
filepath <- file.path('data_extras', paste(filename, ".csv", sep = ""))

Ingalls.Standards <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  filter(Column == Column.Type) %>%
  rename(Metabolite.Name = Compound.Name) %>%
  select(Metabolite.Name,Compound.Type, QE.RF.ratio, Conc..uM, HILICMix, Emperical.Formula) %>%
  filter(!is.na(Conc..uM)) 
Ingalls.Standards$Metabolite.Name <- TrimWhitespace(Ingalls.Standards$Metabolite.Name)


# Import BMIS'd sample file ---------------------------------------------------------------
filename <- RemoveCsv(list.files(path = 'data_processed/', pattern = BMIS.pattern))
filepath <- file.path('data_processed', paste(filename, ".csv", sep = ""))
BMISd.data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE))

# Import QC'd files and remove parameter data ------------------------------
filename <- RemoveCsv(list.files(path = 'data_processed/', pattern = QC.pattern))
filepath <- file.path('data_processed', paste(filename, ".csv", sep = ""))

QCd.data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  select(Replicate.Name, Metabolite.Name, Column, Area.with.QC, Area.Value, Run.Type)

# Import Internal standards key ---------------------------------------------------------------
filename <- RemoveCsv(list.files(path = 'data_extras/', pattern = names.pattern))
filepath <- file.path('data_extras', paste(filename, ".csv", sep = ""))

original.IS.key <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  rename(FinalBMIS = Internal_Standards)


# Identify duplicate compounds that are detected in both HILIC Pos and Neg runs ------------------------------
HILICS.duplicates <- IdentifyDuplicates(QCd.data)

duplicates.testing <- QCd.data %>%
  filter(Metabolite.Name %in% HILICS.duplicates$Metabolite.Name) %>%
  group_by(Metabolite.Name, Column) %>%
  summarize(mean(Area.with.QC, na.rm = TRUE))

QCd.data <- QCd.data %>%
  filter(!(Metabolite.Name %in% HILICS.duplicates$Metabolite.Name & Column == "HILICNeg"))


# Apply appropriate filters and isolate standards ---------------------------------------------------------------
Full.data <- QCd.data %>%
  filter(Metabolite.Name %in% Ingalls.Standards$Metabolite.Name) %>%
  filter(str_detect(Replicate.Name, "Std")) %>%
  left_join(Ingalls.Standards, by = "Metabolite.Name") %>%
  select(Replicate.Name, Metabolite.Name, Compound.Type, everything()) %>%
  unique()

# Remove bad compounds --------------------------------------------------------

# Check standard run types ----------------------------------------------------
Full.data <- CheckStandards(Full.data)

# Get response factors for transect compounds ----------------------------------
Full.data.RF <- Full.data %>%
  mutate(RF = Area.with.QC/Conc..uM) %>%
  filter(!Compound.Type == "Internal Standard") %>%
  mutate(Replicate.Name = substr(Replicate.Name, 1, nchar(Replicate.Name)-2)) %>%
  filter(str_detect(Replicate.Name, as.character(HILICMix)) | str_detect(Replicate.Name, "H2OInMatrix")) 
  

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
Full.data.RF.ratios <- Full.data.RF %>%
  group_by(Metabolite.Name, Type) %>%
  mutate(RF.mean.per_sampleID = mean(RF, na.rm = TRUE)) %>%
  select(Replicate.Name, Metabolite.Name, Type, RF.mean.per_sampleID) %>%
  unique() 

print(paste("NAs or NaNs in the calculatd response factor ratios:", TRUE %in% is.na(Full.data.RF.ratios)))
metabolite.issues <- Full.data.RF.ratios[is.nan(Full.data.RF.ratios$RF.mean.per_sampleID),]
print(unique(metabolite.issues$Metabolite.Name))

Full.data.RF.ratios <- Full.data.RF.ratios %>%
  filter(!(is.infinite(RF.mean.per_sampleID) | is.nan(RF.mean.per_sampleID))) %>%
  group_by(Metabolite.Name) %>% filter(n() >= 3) %>%
  mutate(RF.ratio = ((RF.mean.per_sampleID[Type == "Standards_Matrix"] - RF.mean.per_sampleID[Type == "Water_Matrix"]) / RF.mean.per_sampleID[Type == "Standards_Water"])) %>%
  select(Metabolite.Name, RF.ratio) %>%
  unique()

# If applicable, supplement data with information from calculated Ingalls QE.RF ratios. 

# test.standards <- Ingalls.Standards %>%l
#   filter(Metabolite.Name %in% Full.data.RF.ratios$Metabolite.Name) %>%
#   rename(RF.ratio = QE.RF.ratio) %>%
#   select(Metabolite.Name, RF.ratio)
# 
# Full.data.RF.ratios <- Full.data.RF.ratios %>%
#   as.data.frame() %>%
#   filter(!is.na(RF.ratio)) %>%
#   rbind(test.standards)


# Quantify samples for the BMIS'd dataset ---------------------------------
BMISd.data <- BMISd.data %>%
  #separate(Run.Cmpd, sep = " ", into = c("Sample.Name")) %>%
  separate(Run.Cmpd, c("Sample.Name"), extra = "drop", fill = "right") %>%
  mutate(Metabolite.Name = Mass.Feature) %>%
  filter(Metabolite.Name %in% Full.data.RF.ratios$Metabolite.Name) %>%
  left_join(Full.data.RF.ratios) %>%
  left_join(Full.data.RF.dimensions %>% select(Metabolite.Name, RF.max, RF.min) %>% unique(), by = "Metabolite.Name") %>%
  select(Metabolite.Name, FinalBMIS, Sample.Name, Adjusted.Area, everything())  


## SOMETHING WEIRD IS HAPPENING HERE with the actual concentrations. Check it out. Probably to do with high Taurine RF.Ratio! ###
# Calculate umol/vial for compounds without an internal standard ----------
Quantitative.data <- BMISd.data %>%
  mutate(RF.ave = as.numeric(rowMeans(BMISd.data[, c("RF.min", "RF.max")]))) %>%
  mutate(umol.in.vial.ave = Adjusted.Area/RF.ave/RF.ratio,
         umol.in.vial.max = Adjusted.Area/RF.min/RF.ratio,
         umol.in.vial.min = Adjusted.Area/RF.max/RF.ratio) %>%
  select(Metabolite.Name:Adjusted.Area, everything())

# Pull out data for matched internal standards ----------------------------
IS.key <- BMISd.data %>%
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

Matched.IS.Compounds <- data.frame(Compounds = c(IS.key[ ,"FinalBMIS"], as.character(IS.key[ ,"Metabolite.Name"])))

IS.sample.data <- QCd.data %>%
  left_join(IS.data %>% select(FinalBMIS, Metabolite.Name, Concentration_nM), by = "Metabolite.Name") %>%
  unique() %>%
  filter(Metabolite.Name %in% Matched.IS.Compounds$Compounds) %>%
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

rm(list = c("Matched.IS.Compounds", "QCd.data", "IS.mid_frame", "IS.mid_frame2"))


# Add matched IS_smp info back into main frame ------------------------------------------------
all.info <- Quantitative.data %>% 
  select(Metabolite.Name, runDate:replicate, Adjusted.Area, Area.with.QC, RF.ratio:umol.in.vial.min) %>%
  unite(Sample.Name, c("runDate", "type", "SampID", "replicate"), remove = FALSE) %>%
  left_join(IS.sample.data %>% select(Sample.Name, Metabolite.Name, umol.in.vial_IS)) %>%
  mutate(umol.in.vial.ave = ifelse(is.na(umol.in.vial_IS), umol.in.vial.ave, umol.in.vial_IS),
         umol.in.vial.max = ifelse(is.na(umol.in.vial_IS), umol.in.vial.max, NA),
         umol.in.vial.min = ifelse(is.na(umol.in.vial_IS), umol.in.vial.min, NA)) %>%
  rename(Replicate.Name = Sample.Name) %>%
  select(-runDate, -type, -replicate) %>%
  filter(!str_detect(Replicate.Name, "DDA"))

# Add in dilution factor and filtered volume --------------------------------------------------
all.info.quant <- all.info %>%
  mutate(nmol.in.Enviro.ave = (umol.in.vial.ave*10^-6*Injection.Volume/Volume.Filtered*1000*Dilution.Factor)) %>% 
  left_join(Full.data %>% select(Metabolite.Name, Emperical.Formula)) %>%
  unique()

# Get molecules of carbon and nitrogen ------------------------------------
all.info.molecules <- all.info.quant  %>%
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
all.info.summed <- all.info.molecules %>%
  group_by(Metabolite.Name) %>%
  summarise(nmol.Enviro.med = median(nmol.in.Enviro.ave, na.rm  = T),
            nmol.Enviro.min = min(nmol.in.Enviro.ave, na.rm  = T),
            nmol.Enviro.max = max(nmol.in.Enviro.ave, na.rm  = T),
            nmol.C.med = median(nmol.C.ave, na.rm  = T),
            nmol.C.min = min(nmol.C.ave, na.rm  = T),
            nmol.C.max = max(nmol.C.ave, na.rm  = T)) %>%
  arrange(desc(nmol.Enviro.med))

# Summarize total carbon and nitrogen for each compound ------------------------------------
quantitative.perID <- all.info.molecules %>%
  select(SampID, nmol.C.ave, nmol.N.ave) %>%
  group_by(SampID) %>%
  summarise(totalCmeasured_nM_perID = sum(as.numeric(nmol.C.ave), na.rm = TRUE),
            totalNmeasured_nM_perID = sum(as.numeric(nmol.N.ave), na.rm = TRUE))


# Calculate mole fractions of each compound ------------------------------------
quantitative.final <- all.info.molecules %>%
  unique() %>%
  left_join(quantitative.perID) %>%
  mutate(ratioCN = totalCmeasured_nM_perID / totalNmeasured_nM_perID) %>%
  mutate(molFractionC = nmol.C.ave/totalCmeasured_nM_perID,
         molFractionN = nmol.N.ave/totalNmeasured_nM_perID) %>%
  select(Metabolite.Name, Replicate.Name, Adjusted.Area, Area.with.QC, RF.ratio:molFractionN) %>%
  unique()

quantitative.summed <- quantitative.final %>%
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


# Save and export files ------------------------------------

currentDate <- Sys.Date()
csvFileName.summed <- paste("data_processed/Quantified_Summary_", currentDate, ".csv", sep = "")
csvFileName.final <- paste("data_processed/Quantified_Measurements_", currentDate, ".csv", sep = "")
csvFileName.perID <- paste("data_processed/Quantified_perSampID_", currentDate, ".csv", sep = "")


write.csv(quantitative.summed, csvFileName.summed, row.names = FALSE)
write.csv(quantitative.final, csvFileName.final, row.names = FALSE)
write.csv(quantitative.perID, csvFileName.perID, row.names = FALSE)



