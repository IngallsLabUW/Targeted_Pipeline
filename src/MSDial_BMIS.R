# B-MIS 

# Match QC'd data with Internal Standards list -----------------------------------------------------------------
Data.withIS <- QCd.data %>%
  filter(Metabolite.Name %in% Internal.Standards$Compound.Name)

Data.NoIS <- QCd.data %>%
  filter(!Metabolite.Name %in% Internal.Standards$Compound.Name)

# Create Internal Standard data -----------------------------------------------------------------
Int.Stds.data <- Data.withIS %>%
  select(Replicate.Name, Metabolite.Name, Area.with.QC) %>%
  mutate(Mass.Feature = Metabolite.Name) %>%
  select(-Metabolite.Name) #%>%
#filter(!Mass.Feature == "Guanosine Monophosphate, 15N5")

# Add injection volume -----------------------------------------------------------------
SampKey <- SampKey.all %>%
  filter(Replicate.Name %in% Int.Stds.data$Replicate.Name) %>%
  select(Replicate.Name, Bio.Normalization) %>%
  mutate(Mass.Feature = "Inj_vol",
         Area.with.QC = Bio.Normalization) %>%
  select(Replicate.Name, Area.with.QC, Mass.Feature)

# Create Internal standard data to identify problematic compounds/replicates------------------------------------
Int.Stds.data <- rbind(Int.Stds.data, SampKey) 

# Identify internal standards without an Area, i.e. any NA values.
IS.Issues <- Int.Stds.data[is.na(Int.Stds.data$Area.with.QC), ]
write.csv(IS.Issues, paste("data_intermediate/MSDial_InternalStdIssues_", currentDate, ".csv", sep = ""))


# Visualize raw areas of Internal Standards -----------------------------------------------------------------
IS.Raw.Area.Plot <- ggplot(Int.Stds.data, aes(x = Replicate.Name, y = Area.with.QC)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~Mass.Feature, scales = "free_y") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Internal Standard Raw Areas")

currentDate <- Sys.Date()
plotFileName <- paste("figures/IS.Raw.Areas_", currentDate, ".png", sep = "")

ggsave(file = plotFileName, dpi = 600, width = 8, height = 6, units = "in")
print(IS.Raw.Area.Plot)

# Edit data so names match, test that names are equal across sample sets------------------------------------------
Data.long  <- Data.NoIS %>%
  rename(Mass.Feature = Metabolite.Name) %>%
  select(Replicate.Name, Mass.Feature, Area.with.QC) %>%
  # Uncomment the next line to test the stop error message below.
  # filter(!str_detect(Replicate.Name, regex("dda", ignore_case = TRUE))) %>%
  arrange(Replicate.Name)

test_isdata <- as.data.frame(sort(unique(Int.Stds.data$Replicate.Name)), stringsAsFactors = FALSE)
test_long <- as.data.frame(sort(unique(Data.long$Replicate.Name)), stringsAsFactors = FALSE)
test <- identical(test_isdata[[1]], test_long[[1]])
print(paste("Your replicate names are identical:", test))

if(test == FALSE) 
  stop("Error: Your replicate names are not matched across datasets!")

# Caluclate mean values for each Internal Standard----------------------------------------------------------------
Int.Stds.means <- Int.Stds.data %>%
  mutate(Mass.Feature = as.factor(Mass.Feature)) %>%
  group_by(Mass.Feature) %>%
  summarise(Average.Area = mean(as.numeric(Area.with.QC), na.rm = TRUE)) %>%
  mutate(Mass.Feature = as.character(Mass.Feature))


# Normalize to each internal Standard----------------------------------------------------------------
Data.binded <- rbind(Int.Stds.data, Data.long) %>%
  arrange(Mass.Feature)

Split_Dat <- list()

for (i in 1:length(unique(Int.Stds.data$Mass.Feature))) {
  Split_Dat[[i]] <- Data.binded %>%
    mutate(MIS = unique(Int.Stds.data$Mass.Feature)[i]) %>%
    left_join(Int.Stds.data %>%
                rename(MIS = Mass.Feature, IS_Area = Area.with.QC) %>%
                select(MIS, Replicate.Name, IS_Area), by = c("Replicate.Name", "MIS")) %>%
    left_join(Int.Stds.means %>%
                rename(MIS = Mass.Feature), by = "MIS") %>%
    mutate(Adjusted.Area = Area.with.QC/IS_Area*Average.Area)
}

Data.area.norm <- do.call(rbind, Split_Dat) %>%
  select(-IS_Area, -Average.Area)

# Standardize name structure to: Date_type_ID_replicate_anythingextra ----------------------------------------------------
Mydata.new <- Data.area.norm %>%
  separate(Replicate.Name, c("runDate", "type", "SampID", "replicate"), "_") %>%
  mutate(Run.Cmpd = paste(Data.area.norm$Replicate.Name, Data.area.norm$Mass.Feature))


# Find the B-MIS for each MassFeature----------------------------------------------------------------

# Look only at the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo),
# then choose which IS reduces the RSD the most (Poo.Picked.IS)
Poodata <- Mydata.new %>%
  filter(type == "Poo") %>%
  group_by(SampID, Mass.Feature, MIS) %>%
  summarise(RSD_ofPoo_IND = sd(Adjusted.Area, na.rm = TRUE) / mean(Adjusted.Area, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
  group_by(Mass.Feature, MIS) %>%
  summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo = ifelse(RSD_ofPoo == "NaN", NA, RSD_ofPoo))


Poodata <- Poodata %>%
  left_join(Poodata %>% group_by(Mass.Feature) %>%
              summarise(Poo.Picked.IS = unique(MIS)[which.min(RSD_ofPoo)] [1]))


# Get the original RSD, calculate RSD change, decide if MIS is acceptable ----------------------------------------------------
Poodata <- left_join(Poodata, Poodata %>%
                       filter(MIS == "Inj_vol" ) %>%
                       mutate(Orig_RSD = RSD_ofPoo) %>%
                       select(-RSD_ofPoo, -MIS)) %>%
  mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
  mutate(percent.Change = del_RSD/Orig_RSD) %>%
  mutate(accept_MIS = (percent.Change > cut.off & Orig_RSD > cut.off2))

# Change the BMIS to "Inj_vol" if the BMIS is not an acceptable --------------------------------------------------------------

# Adds a column that has the BMIS, not just Poo.Picked.IS
# Changes the FinalBMIS to inject_volume if it's no good
Fixed.poodata <- Poodata %>%
  filter(MIS == Poo.Picked.IS) %>%
  mutate(FinalBMIS = ifelse(accept_MIS == "FALSE", "Inj_vol", Poo.Picked.IS)) %>%
  mutate(FinalRSD = RSD_ofPoo)

New.poodata <- Poodata %>%
  left_join(Fixed.poodata %>% select(Mass.Feature, FinalBMIS)) %>%
  filter(MIS == FinalBMIS) %>%
  mutate(FinalRSD = RSD_ofPoo)

Try <- New.poodata %>%
  filter(FinalBMIS != "Inj_vol")

QuickReport <- print(paste("Percent of Mass Features that picked a BMIS:",
                           length(Try$Mass.Feature) / length(New.poodata$Mass.Feature), "|",
                           "RSD improvement cutoff", cut.off, "|",
                           "RSD minimum cutoff", cut.off2,
                           sep = " "))

reportFileName = paste("data_intermediate/MSDial_QuickReport", file.pattern, "_", currentDate, ".txt", sep = "")
cat(QuickReport, file = reportFileName)

# Evaluate and visualize the results of your BMIS cutoff----------------------------------------------------------------
IS_toISdat <- Mydata.new %>%
  filter(Mass.Feature %in% Int.Stds.data$Mass.Feature) %>%
  select(Mass.Feature, MIS, Adjusted.Area, type) %>%
  filter(type == "Smp") %>%
  group_by(Mass.Feature, MIS) %>%
  summarize(RSD_of_Smp = sd(Adjusted.Area, na.rm = TRUE) / mean(Adjusted.Area, na.rm = TRUE)) %>%
  left_join(Poodata %>% select(Mass.Feature, MIS, RSD_ofPoo, accept_MIS))

injectONlY_toPlot <- IS_toISdat %>%
  filter(MIS == "Inj_vol")

ISTest_plot <- ggplot() +
  geom_point(dat = IS_toISdat, shape = 21, color = "black", size = 2, aes(x = RSD_ofPoo, y = RSD_of_Smp, fill = accept_MIS)) +
  geom_point(dat = injectONlY_toPlot, aes(x = RSD_ofPoo, y = RSD_of_Smp), size = 3) +
  facet_wrap(~ Mass.Feature) +
  ggtitle(paste("Results of BMIS Cutoff:", cut.off, "RSD decrease,", cut.off2, "RSD minimum."))

plotFileName <- paste("figures/BMIS_Evaluation_", currentDate, ".png", sep = "")

ggsave(file = plotFileName, dpi = 600, width = 8, height = 6, units = "in")
print(ISTest_plot)


# Return data that is normalized via BMIS----------------------------------------------------------------
BMIS.normalized.data <- New.poodata %>% select(Mass.Feature, FinalBMIS, Orig_RSD, FinalRSD) %>%
  left_join(Mydata.new, by = "Mass.Feature") %>%
  filter(MIS == FinalBMIS)

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/MSDial_BMIS_Output_", Column.Type, "_", currentDate, ".csv", sep = "")


write.csv(BMIS.normalized.data, csvFileName, row.names = FALSE)

rm(list = setdiff(ls()[!ls() %in% c("file.pattern")], lsf.str()))