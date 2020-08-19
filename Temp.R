library(tidyverse)

QCd <- read.csv("data_processed/Skyline_QE_QC_Output_Vitamins_2020-08-18.csv", 
                stringsAsFactors = FALSE) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  select(Replicate.Name, Precursor.Ion.Name, Area, Area.with.QC) %>%
  filter(!str_detect(Replicate.Name, "Blk|Std")) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-",".")) %>%
  # filter(!Precursor.Ion.Name == "B2-IS") %>%
  mutate(Nutrient.Addition = ifelse(str_detect(Replicate.Name, "DMB"), "DMB_Added",
                            ifelse(str_detect(Replicate.Name, "WBT"), "B12_Added", NA)))

mytheme <- theme(panel.grid.major = element_line(colour="black", size = (0.03)),
                 panel.grid.minor = element_line(size = (0.2), colour="grey"))

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

print(paste("Total number of peaks:", length(QCd$Area)))
print(paste("Number of NAs thrown out due to area minimum (less than 5000):", sum(is.na(QCd$Area.with.QC))))


