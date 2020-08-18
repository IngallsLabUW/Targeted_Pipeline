library(tidyverse)

QCd <- read.csv("data_processed/Skyline_QE_QC_Output_Vitamins_2020-08-17.csv", 
                stringsAsFactors = FALSE) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, "Blk|Std")) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-",".")) %>%
  filter(!Precursor.Ion.Name == "B2-IS")

mytheme <- theme(panel.grid.major = element_line(colour="black", size = (0.03)),
                 panel.grid.minor = element_line(size = (0.2), colour="grey"))


QCd.plot <- ggplot(QCd, aes(x=Replicate.Name, y=Area.with.QC)) +
  facet_wrap(~Precursor.Ion.Name) +
  geom_bar(stat = "identity") +
  ggtitle("Vitamins: B12 Incubation") +
  theme(axis.text.x = element_blank())
  
print(QCd.plot + mytheme)
