# Imports for BMIS

# Imports -----------------------------------------------------------------

# Internal Standards
Internal.Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
                               stringsAsFactors = FALSE, header = TRUE) %>%
  filter(Column == Column) %>%
  filter(Compound_Type == "Internal Standard")
Internal.Standards$Compound_Name <- TrimWhitespace(Internal.Standards$Compound_Name)

# QC'd output
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = software.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

QCd.data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, "Blk|Std")) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-",".")) ## This line is for an error in the example data!

