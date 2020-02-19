# Skyline TQS Quality Control

# Import datafiles and accompanying master files --------------------------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_raw", pattern = file.pattern))
filepath <- file.path('data_raw', paste(filenames, ".csv", sep = ""))
skyline.output <- assign(make.names(filenames), read.csv(filepath, stringsAsFactors = FALSE)) 

filenames <- RemoveCsv(list.files(path = "data_extras", pattern = "master", ignore.case = TRUE))
filepath <- file.path("data_extras", paste(filenames, ".csv", sep = ""))
master.file <- assign(make.names(filenames), read.csv(filepath, stringsAsFactors = FALSE)) %>%
  rename(Second.Trace = X2nd.trace)

# Sanity check for runtypes and fragments  ---------------------------------------------------------------------
# Stop program if this run has more or fewer runtypes than the normal std, blk, poo, and smp.
skyline.runtypes <- IdentifyRunTypes(skyline.output)

# Stop program if a precursor mz has more than two daughters.
fragments.checked <- CheckFragments(skyline.output, runtype = "Std") 

if (FALSE %in% fragments.checked$Two.Fragments) {
  stop("Some compounds have fewer than two fragments!")
}

# Create comparison tables for parameters  ---------------------------------------------------------------------

# Ion Ratios
# Find Ion Ratio by dividing the area of the quantitative trace by the area of the secondary trace. 
# Find the minimum and maximum IR to create reference table of IR ranges.
ion.ratio.table <- fragments.checked %>%
  group_by(Precursor.Ion.Name, Replicate.Name) %>%
  mutate(Std.Ion.Ratio = ifelse(Quan.Trace == TRUE, (Area[Quan.Trace == TRUE]) / (Area[Second.Trace == TRUE]), NA)) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(IR.min = min(Std.Ion.Ratio, na.rm = TRUE)) %>%
  mutate(IR.max = max(Std.Ion.Ratio, na.rm = TRUE)) %>%
  select(Precursor.Ion.Name, IR.min, IR.max) %>%
  unique()


# Retention Times 
# Find the minimum and maximum Retention Times and take the average.
# Use this as a reference table for acceptable Retention Times.
RT.table <- skyline.output %>%
  filter(str_detect(Replicate.Name, regex("Std", ignore_case = TRUE))) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(RT.min = min(Retention.Time, na.rm = TRUE)) %>%
  mutate(RT.max = max(Retention.Time, na.rm = TRUE)) %>%
  mutate(RT.Reference = mean(Retention.Time, na.rm = TRUE)) %>%
  select(Precursor.Ion.Name, RT.min, RT.max) %>%
  unique()


# Blank Table
# Isolate the blanks in the sample and add a column with maximum blank for each Precursor ion name.
blank.table <- skyline.output %>%
  merge(y = master.file,
        by.x = c("Precursor.Ion.Name", "Product.Mz"),
        by.y = c("Compound.Name", "Daughter"),
        all.x = TRUE) %>%
  mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE),
         Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
  filter(str_detect(Replicate.Name, "Blk"),
         Quan.Trace == TRUE) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Blank.max = max(Area, na.rm = TRUE)) %>%
  select(Precursor.Ion.Name, Blank.max) %>%
  unique()

# # Height
# # Isolate all pooled and sample Heights
# height.table <- skyline.output %>%
#   select(Replicate.Name, Precursor.Ion.Name, Precursor.Mz, Product.Mz, Height) %>%
#   filter(str_detect(Replicate.Name, "Smp|Poo")) 

# Area  
# Isolate all pooled and sample Areas.
# area.table <- skyline.output %>%
#   select(Replicate.Name, Precursor.Ion.Name, Area) %>%
#   filter(str_detect(Replicate.Name, "Smp|Poo"))  


# Signal to Noise 
# Isolate all pooled and sample runs. Find the Signal to Noise
# by dividing the Background of each run by its Area.
# SN.table <- skyline.output %>%
#   merge(y = master.file,
#         by.x = c("Precursor.Ion.Name", "Product.Mz"),
#         by.y = c("Compound.Name", "Daughter"),
#         all.x = TRUE) %>%
#   mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
#   mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
#   filter(str_detect(Replicate.Name, "Smp|Poo")) %>%
#   filter(Quan.Trace == TRUE) %>%
#   select(Replicate.Name, Precursor.Ion.Name, Area, Background) %>%
#   mutate(Signal.to.Noise = (Area / Background)) 

# Construct   ---------------------------------------
# Construct final output of sample and pooled runs.
all.standards <- CheckStdFragments(skyline.output) 

all.samples <- CheckSmpFragments(skyline.output)

all.samples <- all.standards %>%
  left_join(skyline.output %>% filter(str_detect(Replicate.Name, "Smp|Poo")))

# Ion Ratio Flags  ---------------------------------------
# If the Ion Ratio falls outside of the IR.Table range +/- the
# IR.flex value, add a flag.
IR.flags.added <- all.samples %>%
  mutate(IR.Ratio = ifelse(TRUE %in% Significant.Size, (Area[Quan.Trace == TRUE] / Area[Second.Trace == TRUE]), NA)) %>%
  left_join(ion.ratio.table, by = "Precursor.Ion.Name") %>%
  mutate(IR.Flag = ifelse((IR.Ratio < (IR.min - IR.flex) & IR.Ratio > (IR.max + IR.flex)), "IR.Flag", NA)) %>%
  select(Replicate.Name:Second.Trace, Protein.Name:Background, Height, IR.Flag)

# Retention Time Flags  ---------------------------------------
# If the Retention Time is "RT.flex" further away from the RT.Reference 
# Range from the RT.Range Table, add a flag. 
RT.flags.added <- IR.flags.added %>%
  left_join(RT.Range.Table) %>%
  mutate(RT.Flag = ifelse((Retention.Time >= (RT.max + RT.flex) | Retention.Time <= (RT.min - RT.flex)), "RT.Flag", NA)) %>%
  select(Replicate.Name:IR.Flag, RT.Flag)

# Blank Flags  ---------------------------------------
# If the Area divided by the Blank.Reference value is
# greater than the set blk.thresh value, add a flag.
Blank.flags.added <- RT.flags.added %>%
  left_join(Blank.Table) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Blank.Reference = Area * blk.thresh) %>%
  mutate(blank.Flag = ifelse(((Protein.Name != "Internal Std") & (Area * blk.thresh) < Blank.max), 
                             "blank.Flag", 
                             ifelse(((Protein.Name == "Internal Std") & (Area * blk.thresh < Blank.max)), "IS.blank.Flag", NA))) %>%
  select(Replicate.Name:RT.Flag, blank.Flag)


# Height Flags  ---------------------------------------
# Add a height.min.flag if the Height falls below the min.height
# value. Add an overloaded flag if the Height falls above the
# max.height value.
Height.flags.added <- Blank.flags.added %>%
  left_join(Height.Table) %>%
  mutate(height.min.Flag = ifelse((Height < min.height), "height.min.Flag", NA)) %>%
  mutate(overloaded.Flag = ifelse((Height > max.height), "overloaded.Flag", NA)) 

# Area Flags  ---------------------------------------
# If the Area is less than the area.min value, add a flag.
Area.flags.added <- Height.flags.added %>%
  mutate(area.min.Flag = ifelse((Area < area.min), "area.min.Flag", NA)) 


# Signal to Noise Flags  ---------------------------------------
# If the Signal to Noise ratio is less than the SN.min, add a flag.
SN.flags.added <- Area.flags.added %>%
  left_join(SN.Table) %>%
  mutate(SN.Flag = ifelse((Signal.to.Noise < SN.min), "SN.Flag", NA)) 

# All Flags  ---------------------------------------
# Add a column with all flags from the previous steps. 
final.table <- SN.flags.added %>%
  mutate(all.Flags = paste(IR.Flag, RT.Flag, blank.Flag, height.min.Flag, overloaded.Flag, area.min.Flag, SN.Flag, sep = ", ")) %>%
  mutate(all.Flags = as.character(all.Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA")))


# Remove Secondary trace ---------------------------------------
# Filter rows where Second.Trace == TRUE, keeping only Quan.Trace.
# Remove columns once finished.
final.table <- final.table %>%
  filter(Quan.Trace == TRUE) %>%
  select(Replicate.Name:Area, Retention.Time:all.Flags)


# Standards & blank addition  ---------------------------------------
# Test for standards and blanks in the run. Add those standards
# and blanks back into the final table.

Stds.test <- grepl("_Std_", skyline.output$Replicate.Name)
Blks.test <- grepl("_Blk_", skyline.output$Replicate.Name)

if (any(Stds.test == TRUE)) {
  print("There are standards in this run. Joining to the bottom of the dataset!", quote = FALSE)
  ##
  standards <- skyline.output %>%
    filter(str_detect(Replicate.Name, "Std")) %>%
    merge(y = master.file,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    filter(Quan.Trace == TRUE) 
  final.table <- rbind.fill(final.table, standards)
} else {
  print("No standards exist in this set.")
}

if (any(Blks.test == TRUE)) {
  print("There are blanks in this run. Joining to the bottom of the dataset!", quote = FALSE)
  ##
  blanks <- skyline.output %>%
    filter(str_detect(Replicate.Name, "Blk")) %>%
    merge(y = master.file,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    filter(Quan.Trace == TRUE) 
  final.table <- rbind.fill(final.table, blanks)
} else {
  print("No blanks exist in this set.")
}

# Rename and save  ---------------------------------------
# Add comments restating the given QC parameters. Save to 
# current working directory with a new name, 
# "TQSQC_<original file name>.csv

# Print to file with comments and a new name ------------------------------
Description <- c("Hello! Welcome to the world of Skyline TQS Quality Control! ",
                 "Maximum height for a peak: ",
                 "Minimum height for a peak: ",
                 "Maximum area for a real peak: ",
                 "RT flexibility: ",
                 "Blank can be this fraction of a sample: ",
                 "S/N ratio: " ,
                 "Ion ratio flexibility", 
                 "Processed on: ")
Value <- c(NA, max.height, min.height, area.min, RT.flex, blk.thresh, SN.min, IR.flex, Sys.time())

df <- data.frame(Description, Value)
final.table <- bind_rows(df, final.table)


rm(list = ls()[!ls() %in% c("final.table", lsf.str())])

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/QC_Cyano_Output_", currentDate, ".csv", sep = "")

write.csv(final.table, csvFileName, row.names = FALSE)


# 
# con <- file(paste("TQSQC_", basename(input.file), sep = ""), open = "wt")
# writeLines(paste("Hello! Welcome to the world of TQS Quality Control! ",
#                  "Minimum height for a real peak: ", min.height, ". ",
#                  "Minimum area for a real peak: ", area.min, ". ",
#                  "RT flexibility: ", RT.flex, ". ",
#                  "Ion ratio (IR) flexibility: ", IR.flex, ". ",
#                  "Blank can be this fraction of a sample: ", blk.thresh, ". ",
#                  "S/N ratio: " , SN.min, ". ",
#                  "Processed on: ", Sys.time(), ". ",
#                  sep = ""), con)
# write.csv(final.table, con, row.names = FALSE)
# close(con)








