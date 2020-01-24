# Quality control script

# Import files  ----------------------------
filename <- RemoveCsv(list.files(path = 'data_intermediate/', pattern = file.pattern))
filepath <- file.path('data_intermediate', paste(filename, ".csv", sep = ""))

combined <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  select(Replicate.Name:Alignment.ID, Metabolite.Name) %>%
  mutate(Run.Type = (tolower(str_extract(Replicate.Name, "(?<=_)[^_]+(?=_)")))) 

msdial.runtypes <- IdentifyRunTypes(combined)

RT.table <- combined %>%
  filter(Run.Type == "std") %>%
  mutate(RT.Value = na_if(RT.Value, 0)) %>%
  arrange(Metabolite.Name) %>%
  group_by(Metabolite.Name) %>%
  mutate(RT.min = min(RT.Value, na.rm = TRUE)) %>%
  mutate(RT.max = max(RT.Value, na.rm = TRUE)) %>%
  mutate(RT.diff = abs(RT.max - RT.min)) %>%
  select(Metabolite.Name:RT.diff) %>%
  unique()

blank.table <- combined %>%
  filter(Run.Type == "blk") %>%
  mutate(Blk.Area = Area.Value) %>%
  arrange(Metabolite.Name) %>%
  group_by(Metabolite.Name) %>%
  mutate(Blk.min = min(Area.Value)) %>%
  mutate(Blk.max = max(Area.Value)) %>%
  select(Metabolite.Name:Blk.max) %>%
  select(-Blk.Area) %>%
  unique()


# Create signal to noise (SN) and area minimum flags --------------------------------
SN.Area.Flags <- combined %>%
  arrange(Metabolite.Name) %>%
  mutate(SN.Flag       = ifelse(((SN.Value) < SN.min), "SN.Flag", NA)) %>%
  mutate(Area.Min.Flag = ifelse((Area.Value < area.min), "Area.Min.Flag", NA))

# Create retention time flags ---------------------------------------
add.RT.Flag <- SN.Area.Flags %>%
  left_join(RT.table %>% select(-Run.Type)) %>%
  mutate(RT.Flag = ifelse((RT.Value >= (RT.max + RT.flex) | RT.Value <= (RT.min - RT.flex)), "RT.Flag", NA)) %>%
  select(-c("RT.max", "RT.min", "RT.diff"))

# Create blank flags ---------------------------------------
add.blk.Flag <- add.RT.Flag %>%
  left_join(blank.table %>% select(-Run.Type)) %>%
  mutate(Blank.Flag = ifelse((Area.Value / Blk.max) < blk.thresh, "Blank.Flag", NA)) %>%
  select(-c("Blk.min", "Blk.max"))


# Combine all the flags ---------------------------------------------------
final.table <- add.blk.Flag %>%
  mutate(all.Flags      = paste(SN.Flag, Area.Min.Flag, RT.Flag, Blank.Flag, sep = ", ")) %>%
  mutate(all.Flags      =all.Flags %>% str_remove_all("NA, ") %>% str_remove_all("NA")) %>%
  mutate(all.Flags      = ifelse(all.Flags == "", NA, all.Flags)) %>%
  mutate(Area.with.QC   = ifelse(is.na(Area.Min.Flag), Area.Value, NA)) %>%
  select(Replicate.Name:Area.Value, Area.with.QC, everything())


# Print to file with comments and a new name ------------------------------
Description <- c(as.character(anydate(Sys.Date())),
                "Hello! Welcome to the world of MSDIAL QE Quality Control! ",
                 "Minimum area for a real peak: ",
                 "RT flexibility: ",
                 "Blank can be this fraction of a sample: ",
                 "S/N ratio: ")
Value <- as.character(c(NA, NA, area.min, RT.flex, blk.thresh, SN.min))

df <- data.frame(Description, Value)
final.table <- bind_rows(df, final.table)

rm(list = setdiff(ls()[!ls() %in% c("file.pattern", "final.table", "RT.table", "blank.table")], lsf.str()))