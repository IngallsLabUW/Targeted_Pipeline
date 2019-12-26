# Things to Return --------------------------------------------------------
source("Functions.R")

## OUTPUT IS MESSED UP WITH SOME NAS, CHECK WHYYYYYYYYY


# IS_inspectPlot (plot to make sure there aren't any internal standards we should kick out)
# QuickReport (% that picked BMIS, with cut off values)
# ISTest_plot (plot to evaluate if you cut off is appropriate)
# BMIS_normalizedData (tibble with the info you actually want!)

# Set parameters -----------------------------------------------------------------
cut.off <- 0.3 # 30% decrease in RSD of pooled injections, aka improvement cutoff
cut.off2 <- 0.1 # RSD minimum
Column.Type = "HILIC"
pattern = "QC"

# Imports -----------------------------------------------------------------
# Sample Key
SampKey.all <- read.csv("data_extras/Sample.Key.EddyTransect.csv") %>%
  rename(Replicate.Name = Sample.Name) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-","."))

# Internal Standards
Internal.Standards <- read.csv("data_extras/Ingalls_Lab_Standards.csv") %>%
  filter(Column == Column.Type) %>%
  filter(Compound.Type == "Internal Standard")

Internal.Standards$Compound.Name <- TrimWhitespace(Internal.Standards$Compound.Name)

# QC'd HILIC output
filename <- RemoveCsv(list.files(path = 'data_processed/', pattern = pattern))
filepath <- file.path('data_processed', paste(filename, ".csv", sep = ""))

HILIC.QC <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, "Blk|Std")) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-",".")) 


HILICS.duplicates <- IdentifyDuplicates(HILIC.QC)

HILIC.QC <- HILIC.QC %>%
  filter(!(Metabolite.Name %in% HILICS.duplicates$Metabolite.Name & Column == "HILICNeg"))


# Match QC'd HILIC data with Internal Standards list -----------------------------------------------------------------
HILIC.withIS <- HILIC.QC %>%
  filter(Metabolite.Name %in% Internal.Standards$Compound.Name)

HILIC.NoIS <- HILIC.QC %>%
  filter(!Metabolite.Name %in% Internal.Standards$Compound.Name)

# Create HILICs Internal Standard data -----------------------------------------------------------------
HILIC.IS.data <- HILIC.withIS %>%
  select(Replicate.Name, Metabolite.Name, Area.with.QC) %>%
  mutate(Mass.Feature = Metabolite.Name) %>%
  select(-Metabolite.Name) 
#filter(!MassFeature == "Guanosine Monophosphate, 15N5")

# Add injection volume -----------------------------------------------------------------
SampKey <- SampKey.all %>%
  filter(Replicate.Name %in% HILIC.IS.data$Replicate.Name) %>% 
  select(Replicate.Name, Bio.Normalization) %>%
  mutate(Mass.Feature = "Inj_vol",
         Area.with.QC = Bio.Normalization) %>%
  select(Replicate.Name, Area.with.QC, Mass.Feature)

# Create Internal standard data to identify problematic compounds/replicates-----------------------------------------------------------------
HILIC.IS.data <- rbind(HILIC.IS.data, SampKey) %>%
  filter(!str_detect(Replicate.Name, "dda"))
# HILIC.IS.data[] <- lapply(HILIC.IS.data, gsub, pattern = 'Neg|Pos', replacement = '')

# Here is where we would hypothetically remove troublesome compounds.


# Identify internal standards without an Area, i.e. any NA values.
IS.Issues <- HILIC.IS.data[is.na(HILIC.IS.data$Area.with.QC), ]

# Visualize raw areas of Internal Standards -----------------------------------------------------------------
IS.Raw.Area.Plot <- ggplot(HILIC.IS.data, aes(x = Replicate.Name, y = Area.with.QC)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~Mass.Feature, scales = "free_y") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Internal Standard Raw Areas")
print(IS.Raw.Area.Plot)


# Edit data so names match-----------------------------------------------------------------
HILIC.long  <- HILIC.NoIS %>%
  rename(Mass.Feature = Metabolite.Name) %>%
  select(Replicate.Name, Mass.Feature, Area.with.QC) %>%
  filter(!str_detect(Replicate.Name, "dda")) %>%
  arrange(Replicate.Name)

# Test that names are equal across sample sets-----------------------------------------------------------------
test_isdata <- as.data.frame(sort(unique(HILIC.IS.data$Replicate.Name)), stringsAsFactors = FALSE)
test_long <- as.data.frame(sort(unique(HILIC.long$Replicate.Name)), stringsAsFactors = FALSE)
identical(test_isdata[[1]], test_long[[1]])

# Caluclate mean values for each Internal Standard----------------------------------------------------------------
HILIC.IS.means <- HILIC.IS.data %>%
  filter(!grepl("_Blk_", Replicate.Name)) %>%
  mutate(Mass.Feature = as.factor(Mass.Feature)) %>%
  group_by(Mass.Feature) %>%
  summarise(Average.Area = mean(as.numeric(Area.with.QC), na.rm = TRUE)) %>%
  mutate(Mass.Feature = as.character(Mass.Feature))

HILIC.IS.means[is.na(HILIC.IS.means)] <- NA


# Normalize to each internal Standard----------------------------------------------------------------
HILIC.binded <- rbind(HILIC.IS.data, HILIC.long) %>%
  arrange(Mass.Feature)

Split_Dat <- list()

for (i in 1:length(unique(HILIC.IS.data$Mass.Feature))) {
  Split_Dat[[i]] <- HILIC.binded %>%
    mutate(MIS = unique(HILIC.IS.data$Mass.Feature)[i]) %>%
    left_join(HILIC.IS.data %>%
                rename(MIS = Mass.Feature, IS_Area = Area.with.QC) %>%
                select(MIS, Replicate.Name, IS_Area), by = c("Replicate.Name", "MIS")) %>%
    left_join(HILIC.IS.means %>%
                rename(MIS = Mass.Feature), by = "MIS") %>%
    mutate(Adjusted.Area = Area.with.QC/IS_Area*Average.Area)
}

HILIC.area.norm <- do.call(rbind, Split_Dat) %>%
  select(-IS_Area, -Average.Area)

# Standardize name structure to: Date_type_ID_replicate_anythingextra) ----------------------------------------------------------------
HILIC.mydata.new <- HILIC.area.norm %>%
  separate(Replicate.Name, c("runDate", "type", "SampID", "replicate"), "_") %>%
  mutate(Run.Cmpd = paste(HILIC.area.norm$Replicate.Name, HILIC.area.norm$Mass.Feature))


# Find the B-MIS for each MassFeature----------------------------------------------------------------

# Look only at the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo),
# then choose which IS reduces the RSD the most (Poo.Picked.IS)
HILIC.poodat <- HILIC.mydata.new %>%
  filter(type == "Poo") %>%
  group_by(SampID, Mass.Feature, MIS) %>%
  summarise(RSD_ofPoo_IND = sd(Adjusted.Area, na.rm = TRUE) / mean(Adjusted.Area, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
  group_by(Mass.Feature, MIS) %>%
  summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo = ifelse(RSD_ofPoo == "NaN", NA, RSD_ofPoo)) 


HILIC.poodat <- HILIC.poodat %>%
  left_join(HILIC.poodat %>% group_by(Mass.Feature) %>%
              summarise(Poo.Picked.IS = unique(MIS)[which.min(RSD_ofPoo)] [1]))


# Get the original RSD, calculate RSD change, decide if MIS is acceptable----------------------------------------------------------------
HILIC.poodat <- left_join(HILIC.poodat, HILIC.poodat %>%
                            filter(MIS == "Inj_vol" ) %>%
                            mutate(Orig_RSD = RSD_ofPoo) %>%
                            select(-RSD_ofPoo, -MIS)) %>%
  mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
  mutate(percent.Change = del_RSD/Orig_RSD) %>%
  mutate(accept_MIS = (percent.Change > cut.off & Orig_RSD > cut.off2))


# Change the BMIS to "Inj_vol" if the BMIS is not an acceptable----------------------------------------------------------------

# Adds a column that has the BMIS, not just Poo.Picked.IS
# Changes the FinalBMIS to inject_volume if its no good
HILIC.fixedpoodat <- HILIC.poodat %>%
  filter(MIS == Poo.Picked.IS) %>% 
  mutate(FinalBMIS = ifelse(accept_MIS == "FALSE", "Inj_vol", Poo.Picked.IS)) %>%
  mutate(FinalRSD = RSD_ofPoo)

HILIC.newpoodat <- HILIC.poodat %>%
  left_join(HILIC.fixedpoodat %>% select(Mass.Feature, FinalBMIS)) %>%
  filter(MIS == FinalBMIS) %>%
  mutate(FinalRSD = RSD_ofPoo)

Try <- HILIC.newpoodat %>%
  filter(FinalBMIS != "Inj_vol")

QuickReport <- print(paste("% of MFs that picked a BMIS",
                           length(Try$Mass.Feature) / length(HILIC.newpoodat$Mass.Feature),
                           "RSD improvement cutoff", cut.off,
                           "RSD minimum cutoff", cut.off2,
                           sep = " "))


# Evaluate and visualize the results of your BMIS cutoff----------------------------------------------------------------
IS_toISdat <- HILIC.mydata.new %>%
  filter(Mass.Feature %in% HILIC.IS.data$Mass.Feature) %>%
  select(Mass.Feature, MIS, Adjusted.Area, type) %>%
  filter(type == "Smp") %>%
  group_by(Mass.Feature, MIS) %>%
  summarise(RSD_ofSmp = sd(Adjusted.Area, na.rm = TRUE)/mean(Adjusted.Area, na.rm = TRUE)) %>%
  left_join(HILIC.poodat %>% select(Mass.Feature, MIS, RSD_ofPoo, accept_MIS))

injectONlY_toPlot <- IS_toISdat %>%
  filter(MIS == "Inj_vol")


ISTest_plot <- ggplot() +
  geom_point(dat = IS_toISdat, shape = 21, color = "black", size = 2,aes(x = RSD_ofPoo, y = RSD_ofSmp, fill = accept_MIS)) +
  scale_fill_manual(values=c("white","dark gray")) +
  geom_point(dat = injectONlY_toPlot, aes(x = RSD_ofPoo, y = RSD_ofSmp), size = 3) +
  facet_wrap(~ Mass.Feature)
print(ISTest_plot)


# Return data that is normalized via BMIS----------------------------------------------------------------

## original
HILIC.BMIS_normalizedData <- HILIC.newpoodat %>% select(Mass.Feature, FinalBMIS, Orig_RSD, FinalRSD) %>%
  left_join(HILIC.mydata.new, by = "Mass.Feature") %>%
  filter(MIS == FinalBMIS) %>%
  unique()

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/BMIS_Output_", currentDate, ".csv", sep = "")


write.csv(HILIC.BMIS_normalizedData, csvFileName, row.names = FALSE)

rm(list = ls())






