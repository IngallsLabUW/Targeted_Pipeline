# Skyline TQS + QE Quality Control


# Import datafiles and accompanying master files --------------------------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_intermediate", pattern = file.pattern))
filepath <- file.path("data_intermediate", paste(filenames, ".csv", sep = ""))
skyline.output <- assign(make.names(filenames), read.csv(filepath, stringsAsFactors = FALSE)) 

if (instrument.pattern == "TQS") {
  filenames <- RemoveCsv(list.files(path = "data_extras", pattern = "master", ignore.case = TRUE))
  filepath <- file.path("data_extras", paste(filenames, ".csv", sep = ""))
  master.file <- assign(make.names(filenames), read.csv(filepath, stringsAsFactors = FALSE)) %>%
    dplyr::rename(Second.Trace = X2nd.trace)
}

# Sanity check for runtypes  ---------------------------------------------------------------------
# Stop program if this run has more or fewer runtypes than the normal std, blk, poo, and smp.
skyline.runtypes <- IdentifyRunTypes(skyline.output)

# # In HILIC compounds, filter mixes.
# if ("Column" %in% colnames(Full.data.RF)) {
#   Full.data.RF <- Full.data.RF %>%
#     filter(str_detect(Replicate.Name, as.character(HILICMix)) | str_detect(Replicate.Name, "H2OInMatrix"))
# }


# Depending on instrument.pattern, create comparison tables --------------------------------------

if (instrument.pattern == "TQS") {
  # Check for fragments in TQS data.
  fragments.checked <- CheckFragments(skyline.output, runtype = "Std") 
  
  # Stop program if a precursor mz has more than two daughters.
  if (FALSE %in% fragments.checked$Two.Fragments) {
    stop("Some compounds have fewer than two fragments!")
  }
  
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
  
  # Height
  # Isolate all pooled and sample Heights
  height.table <- skyline.output %>%
    select(Replicate.Name, Precursor.Ion.Name, Precursor.Mz, Product.Mz, Height) %>%
    filter(str_detect(Replicate.Name, "Smp|Poo")) 
  
  # Signal to Noise 
  # Isolate all pooled and sample runs. Find the Signal to Noise
  # by dividing the Background of each run by its Area.
  SN.table <- skyline.output %>%
    merge(y = master.file,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    filter(str_detect(Replicate.Name, "Smp|Poo")) %>%
    filter(Quan.Trace == TRUE) %>%
    select(Replicate.Name, Precursor.Ion.Name, Area, Background) %>%
    mutate(Signal.to.Noise = (Area / Background))
  
} else{
  
  print(paste("This is a", instrument.pattern, "run. No fragmentation check necessary."))
  
  # Blank Table
  blank.table <- skyline.output %>%
    filter(str_detect(Replicate.Name, regex("Blk", ignore_case = TRUE))) %>%
    select(Precursor.Ion.Name, Area) %>%
    group_by(Precursor.Ion.Name) %>% 
    filter(Area == max(Area)) %>%
    rename(Blank.max = Area) %>%
    unique()
  
  # Height
  # Isolate all pooled and sample heights
  height.table <- skyline.output %>%
    select(Replicate.Name, Precursor.Ion.Name, Height) %>%
    filter(str_detect(Replicate.Name, regex("Smp|Poo", ignore_case = TRUE))) 
}


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

# Area  
# Isolate all pooled and sample Areas.
area.table <- skyline.output %>%
  select(Replicate.Name, Precursor.Ion.Name, Area) %>%
  filter(str_detect(Replicate.Name, regex("Smp|Poo", ignore_case = TRUE)))

# Signal to Noise 
# Isolate all pooled and sample runs. Find the Signal to Noise
# by dividing the Background of each run by its Area.
SN.table <- skyline.output %>%
  filter(str_detect(Replicate.Name, regex("Smp|Poo", ignore_case = TRUE))) %>%
  select(Replicate.Name, Precursor.Ion.Name, Area, Background) %>%
  mutate(Signal.to.Noise = (Area / Background))


# Construct final comparative table ---------------------------------------

if (instrument.pattern == "TQS") {
  all.standards <- CheckFragments(skyline.output, runtype = "Std") 
  
  all.samples <- CheckFragments(skyline.output, runtype = "Smp")
  
  all.samples <- all.samples %>%
    left_join(skyline.output %>% filter(str_detect(Replicate.Name, "Smp|Poo")))
  
  # Ion Ratio Flags  ---------------------------------------
  # If the Ion Ratio falls outside of the IR.Table range +/- the
  # IR.flex value, add a flag.
  all.samples <- all.samples %>%
    group_by(Precursor.Ion.Name) %>%
    mutate(IR.Ratio = ifelse(TRUE %in% Significant.Size, (Area[Quan.Trace == TRUE] / Area[Second.Trace == TRUE]), NA)) %>%
    left_join(ion.ratio.table, by = "Precursor.Ion.Name") %>%
    mutate(IR.Flag = ifelse(((IR.Ratio < (IR.min - IR.flex)) | (IR.Ratio > (IR.max + IR.flex))), "IR.Flag", NA)) %>%
    select(Replicate.Name:Second.Trace, Protein.Name:Background, Height, IR.Flag)
  
} else {
  
  all.samples <- skyline.output
  
}

# Retention Time Flags  ---------------------------------------
# If the Retention Time is "RT.flex" further away from the RT.Reference 
# Range from the RT.Range Table, add a flag. 
RT.flags.added <- all.samples %>%
  left_join(RT.table) %>%
  mutate(RT.Flag = ifelse((Retention.Time >= (RT.max + RT.flex) | Retention.Time <= (RT.min - RT.flex)), 
                          "RT.Flag", NA))

# Blank Flags  ---------------------------------------
# If the Area divided by the Blank.Reference value is
# greater than the set blk.thresh value, add a flag.
Blank.flags.added <- RT.flags.added %>%
  left_join(blank.table) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Blank.Reference = Area / Blank.max) %>%
  ########################333
  mutate(Protein.Name = ifelse(str_detect(Replicate.Name, ","), "Internal Std", "Non IS")) %>%
  ##############################
  mutate(blank.Flag = ifelse(((Protein.Name != "Internal Std") & ((Area / Blank.max) < blk.thresh)), 
                             "blank.Flag", 
                             ifelse(((Protein.Name == "Internal Std") & ((Area / Blank.max) < blk.thresh)),
                                    "IS.blank.Flag", NA)))


# Height Flags  ---------------------------------------
# Add a height.min.flag if the Height falls below the min.height
# value. Add an overloaded flag if the Height falls above the
# max.height value.
Height.flags.added <- Blank.flags.added %>%
  left_join(height.table) %>%
  mutate(height.min.Flag = ifelse((Height < height.min), "height.min.Flag", NA)) %>%
  mutate(overloaded.Flag = ifelse((Height > height.max), "overloaded.Flag", NA))

# Area Flags  ---------------------------------------
# If the Area is less than the area.min value, add a flag.
Area.flags.added <- Height.flags.added %>%
  mutate(area.min.Flag = ifelse((Area < area.min), "area.min.Flag", NA)) %>%
  mutate(Area.with.QC   = ifelse(is.na(area.min.Flag), Area, NA)) %>%
  select(Replicate.Name:Area, Area.with.QC, everything())

# Signal to Noise Flags  ---------------------------------------
# If the Signal to Noise ratio is less than the SN.min, add a flag.
SN.flags.added <- Area.flags.added %>%
  left_join(SN.table) %>%
  mutate(SN.Flag = ifelse((Signal.to.Noise < SN.min), "SN.Flag", NA))

# All Flags  ---------------------------------------
# Add a column with all flags from the previous steps. 
semifinal.table <- SN.flags.added 

semifinal.table <- semifinal.table %>%
  unite(all.Flags, contains("Flag"), sep = ", ", remove = FALSE) %>%
  mutate(all.Flags = as.character(all.Flags %>% str_remove_all("NA,|NA") %>% trimws()))
semifinal.table$all.Flags <- gsub('^\\,|\\,$', '', semifinal.table$all.Flags)

final.table <- semifinal.table %>%
  select(Replicate.Name:Mass.Error.PPM, contains("Flag"))
final.table[final.table==""]<-NA


# Remove Secondary trace ---------------------------------------
# Filter rows where Second.Trace == TRUE, keeping only Quan.Trace.
# Remove columns once finished.
if (instrument.pattern == "TQS") {
  final.table <- final.table %>%
    filter(Quan.Trace == TRUE) %>%
    select(Replicate.Name:Area, Retention.Time:all.Flags)
}


# Standards & blank addition  ---------------------------------------
# Test for standards and blanks in the run. Add those standards
# and blanks back into the final table.
# 
# Stds.test <- grepl("_Std_", skyline.output$Replicate.Name)
# Blks.test <- grepl("_Blk_", skyline.output$Replicate.Name)
# 
# if (any(Stds.test == TRUE)) {
#   print("There are standards in this run. Joining to the bottom of the dataset!", quote = FALSE)
#   final.table <- JoinStandardsBlanks(skyline.output, machine = instrument.pattern, runtype = "Std")
# } else {
#   print("No standards exist in this set.")
# }
# 
# if (any(Blks.test == TRUE)) {
#   print("There are blanks in this run. Joining to the bottom of the dataset!", quote = FALSE)
#   final.table <- JoinStandardsBlanks(skyline.output, machine = instrument.pattern, runtype = "Blk")
# } else {
#   print("No blanks exist in this set.")
# }

# Rename and save  ---------------------------------------
# Add comments restating the given QC parameters. Save to 
# current working directory with a new name, 
# "TQSQC_<original file name>.csv

# Print to file with comments and a new name ------------------------------
if (instrument.pattern == "TQS") {
  Description <- c(as.character(anydate(Sys.Date())),
                   "Hello! Welcome to the world of Skyline TQS Quality Control! ",
                   "Maximum height for a real peak: ",
                   "Minimum height for a real peak: ",
                   "Maximum area for a real peak: ",
                   "RT flexibility: ",
                   "Blank can be this fraction of a sample: ",
                   "S/N ratio: " ,
                   "Ion ratio flexibility", 
                   "Processed on: ")
  
  Value <- c(NA, NA, height.max, height.min, area.min, RT.flex, blk.thresh, SN.min, IR.flex, Sys.time())
} else{
  Description <- c(as.character(anydate(Sys.Date())),
                   "Hello! Welcome to the world of Skyline QE Quality Control! ",
                   "Maximum height for a real peak: ",
                   "Minimum height for a real peak: ",
                   "Maximum area for a real peak: ",
                   "RT flexibility: ",
                   "Blank can be this fraction of a sample: ",
                   "S/N ratio: " ,
                   "Processed on: ")
  Value <- c(NA, NA, height.max, height.min, area.min, RT.flex, blk.thresh, SN.min, Sys.time())
  
}

df <- data.frame(Description, Value)
final.table <- bind_rows(df, final.table)


rm(list = setdiff(ls()[!ls() %in% c("software.pattern", "file.pattern", "instrument.pattern",
                                    "final.table", "ion.ratio.table", "RT.table", "blank.table",
                                    "height.table", "area.table", "SN.table")], lsf.str()))









