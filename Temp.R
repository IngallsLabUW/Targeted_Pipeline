library(tidyverse)

mytheme <- theme(panel.grid.major = element_line(colour="black", size = (0.03)),
                 panel.grid.minor = element_line(size = (0.2), colour="grey"))

QCd <- read.csv("data_processed/Skyline_QE_QC_Output_Vitamins_2020-08-19.csv", 
                stringsAsFactors = FALSE) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  select(Replicate.Name, Precursor.Ion.Name, Area, Area.with.QC) %>%
  filter(!str_detect(Replicate.Name, "Blk|Std")) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-",".")) %>%
  filter(!str_detect(Replicate.Name, 
                     "TruePoo|CultureMED4|Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m|Process|Std")) %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  mutate(Control.Status = ifelse(str_detect(Supergroup, "T0"),
                                 "Incubation", ifelse(str_detect(Supergroup, "DSW"), "DeepSeaWater", 
                                                      ifelse(str_detect(Supergroup, "Control"), "Control", "Treatments")))) %>%
  mutate(Treatment.Status = ifelse(Control.Status == "Control", "Control",
                                   ifelse(Control.Status == "DeepSeaWater", "DeepSeaWater",
                                          ifelse(Control.Status == "Incubation", "TimeZero",
                                                 ifelse(str_detect(Supergroup, "DMBnoBT"), "DMBnoB12",
                                                        ifelse(str_detect(Supergroup, "WBT"), "B12",
                                                               ifelse(str_detect(Supergroup, "DMB"), "DMB", "noB12"))))))) %>%
  select(Replicate.Name, Supergroup, Precursor.Ion.Name, Area.with.QC, Treatment.Status)



all.plot <- ggplot(QCd, aes(x=Replicate.Name, y=Area, fill = Nutrient.Addition)) +
  facet_wrap(~Precursor.Ion.Name, scales = "free") +
  geom_bar(stat = "identity") +
  ggtitle("Vitamins: B12 Incubation, all peaks") +
  theme(axis.text.x = element_blank())
print(all.plot + mytheme)

QCd.plot <- ggplot(QCd, aes(x=Replicate.Name, y=Area.with.QC, fill = Nutrient.Addition)) +
  facet_wrap(~Precursor.Ion.Name, scales = "free") +
  geom_bar(stat = "identity") +
  ggtitle("Vitamins: B12 Incubation, QC'd") +
  theme(axis.text.x = element_blank())
print(QCd.plot + mytheme)

NonQCd <- QCd %>%
  filter(is.na(Area.with.QC)) 

NonQCd.plot <- ggplot(NonQCd, aes(x=Replicate.Name, y=Area, fill=Nutrient.Addition)) +
  facet_wrap(~Precursor.Ion.Name) +
  geom_bar(stat = "identity") +
  ggtitle("Vitamins: B12 Incubation, ONLY nonQC PEAKS") +
  theme(axis.text.x = element_blank())
print(NonQCd.plot + mytheme)



all.treatments.plot <- ggplot(all.treatments, aes(x=Replicate.Name, y=Area.with.QC, fill=Treatment.Status)) +
  facet_wrap(~Precursor.Ion.Name, scales = "free") +
  geom_bar(stat = "identity") +
  ggtitle("Vitamins: B12 Incubation. All Treatments.") +
  theme(axis.text.x = element_blank())
print(all.treatments.plot + mytheme)

print(paste("Total number of peaks:", length(QCd$Area)))
print(paste("Number of NAs thrown out due to area minimum (less than 5000):", sum(is.na(QCd$Area.with.QC))))


