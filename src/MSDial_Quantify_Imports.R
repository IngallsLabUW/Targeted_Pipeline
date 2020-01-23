# Imports for the quantification step


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
  select(-c(Description, Value, Mz.Value, RT.Value, SN.Value, Alignment.ID, contains("Flag"))) %>%
  select(Replicate.Name, Metabolite.Name, Area.with.QC, Area.Value, Run.Type, everything())

# Import Internal standards key ---------------------------------------------------------------
filename <- RemoveCsv(list.files(path = 'data_extras/', pattern = names.pattern))
filepath <- file.path('data_extras', paste(filename, ".csv", sep = ""))

original.IS.key <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  rename(FinalBMIS = Internal_Standards)

# Apply appropriate filters and isolate standards ---------------------------------------------------------------
Full.data <- QCd.data %>%
  filter(Metabolite.Name %in% Ingalls.Standards$Metabolite.Name) %>%
  filter(str_detect(Replicate.Name, "Std")) %>%
  left_join(Ingalls.Standards, by = "Metabolite.Name") %>%
  select(Replicate.Name, Metabolite.Name, Compound.Type, everything()) %>%
  unique()
