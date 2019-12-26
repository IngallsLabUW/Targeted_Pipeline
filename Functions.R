## Function definitions ##
library(rlist)
library(tidyverse)
library(tidyr)
options(scipen=999)


# STILL NEEDS SOME WORK
ChangeClasses <- function(df) {
  # Identifies columns starting with X and changes their class to numeric.
  #
  # Args
  #   df: MSDial dataframe reorganized to drop all empty rows at the top.
  #
  # Returns
  #   df: MSDial dataframed with modified sample column classes.
  #
  #col.test <- grepl("^X", names(df))
  for (i in c(10:ncol(df))) {
    df[, i] <- as.numeric(as.character(df[, i]))
  }
  return(df)
}

RearrangeDatasets <- function(df, parameter) {
  # Shortcut for altering multiple datasets using the tidyr::gather() function.
  #
  # Args
  #   df: MSDial dataframe with first n empty rows removed.
  #   parameter: Table value. This parameter will become the column name when 
  #              changed to long format.
  #
  # Returns
  #   df: MSDial dataframe, changed to long format and with a custom-named value column.
  df <- df %>%
    tidyr::gather(
      key = "Replicate.Name",
      value = "parameter",
      starts_with("X")) %>%
    select(Replicate.Name, parameter, everything())
  
  names(df)[2] <- parameter
  
  return(df)
}

RemoveCsv <- function(full.filepaths) {
  # Gathers all files in given directory and drops the csv extension.
  #
  # Args
  #   full.filepaths: list of files in a directory matching given patterns.
  #
  # Returns
  #   no.path: list of files, with filepath and csv extension removed.
  #
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}

SetHeader <- function(df) {
  # Test for blank rows in MSDial output, and filter them out. Replace column names with syntactically correct headers.
  #
  # Args
  #   df: MSDial dataframe. 
  #
  # Returns
  #   df: MSDial dataframe, first n blank rows removed and headers set.
  #
  df <- df[!(is.na(df[1]) | df[1] == ""), ]
  colnames(df) <- make.names(as.character(unlist(df[1,])))
  df <- df[-1, ]
  
  return(df)
}






CheckStandards <- function (df) {
  # Mutates a new column identifying standard run types, then prints number of unique run types.
  #
  # Args
  #   df: Dataset of containing a Replicate.Name column, pre-filtered to include only standard runs.
  #
  # Returns
  #   df.checked: Dataset with a new column describing run types, and a printed message stating how many 
  #               unique types there are.
  #
  df.checked <- df %>%
    mutate(Type = paste(Env = ifelse(str_detect(Replicate.Name, "StdsMix|InH2O"), "Standards", "Water"),
                        Matrix = ifelse(str_detect(Replicate.Name, "InMatrix"), "Matrix", "Water"), sep = "_"))
  
  print(paste("Number of standard run types:", length(unique(df.checked$Type))))
  print(unique(df.checked$Type))
  
  return(df.checked)
}


StandardizeVariables <- function(df) {
  if (c("ReplicateName", "AreaValue", "MZValue", "RTValue", "SNValue") %in% colnames(df))
  {
    df <- df %>%
      rename(Replicate.Name = ReplicateName) %>%
      rename(Area.Value = AreaValue) %>%
      rename(MZ.Value = MZValue) %>%
      rename(RT.Value = RTValue) %>%
      rename(SN.Value = SNValue)
  }
  return(df)
}

TrimWhitespace <- function (x) gsub("^\\s+|\\s+$", "", x)

CheckBlankMatcher <- function(blank.matcher) {
  # Takes a blank matcher file and separates any multi-value variable
  # columns into their own row.
  #
  # Args:
  #   blank.matcher: CSV entered by user to match samples with
  #   appropriate blanks.
  #
  # Returns:
  #   blank.matcher: new CSV with any duplicate values separated
  #   into their own rows.
  #
  blank.matcher <- do.call("rbind", Map(data.frame,
                                        Blank.Name = strsplit(as.character(blank.matcher$Blank.Name), ","),
                                        Replicate.Name = (blank.matcher$Replicate.Name))
  )
  blank.matcher <- blank.matcher[c(2, 1)]
  
  return(blank.matcher)
}

IdentifyRunTypes <- function(msdial.file) {
  # Identify run typfes and return each unique value present in the Skyline output.
  #
  # Args
  #   msdial.file: Raw output file from Skyline.
  #
  # Returns
  #   run.type: list of labels identifying the run types, isolated from Replicate.Name.
  #   Options conssist of samples (smp), pooled (poo), standards (std), and blanks (blk).
  #
  run.type <- tolower(str_extract(msdial.file$Replicate.Name, "(?<=_)[^_]+(?=_)"))
  print(paste("Your runtypes are:", toString(unique(run.type))))
  return(run.type)
}



StandardizeMetabolites <- function(df) {
  df.standardized <- df %>%
    mutate(Metabolite.Name = ifelse(str_detect(Metabolite.Name, "Ingalls_"), sapply(strsplit(Metabolite.Name, "_"), `[`, 2), Metabolite.Name)) 
  
  df.standardized$Replicate.Name <- gsub("^.{0,1}", "", df.standardized$Replicate.Name)
  
  return(df.standardized)
}

IdentifyDuplicates <- function(df) {
  # Determine which compounds are detected in both positive and negative HILIC runs.
  # 
  # Args
  #   df: MSDial dataframe, containing all required parameters (MZ, SN, Area, etc),
  #       and modified to long form instead of wide.
  # 
  # Returns
  #   duplicates: Simple dataframe of listed compounds that have been identified as duplicates.
  #
  duplicates <- df %>%
    group_by(Metabolite.Name, Replicate.Name) %>%
    mutate(number = 1) %>%
    mutate(ticker = cumsum(number)) %>%
    filter(ticker == 2) %>%
    ungroup() %>%
    select(Metabolite.Name) %>%
    unique()
  return(duplicates)
}



# Old Functions -----------------------------------------------------------

# FilterUnknowns <- function(df) {
#   df.filtered <- df %>%  
#     filter_(!Metabolite.name == "Unknown") #%>%
#     # select(-c(Average.Rt.min., Formula, Ontology, INCHIKEY, SMILES, Isotope.tracking.parent.ID, Isotope.tracking.weight.number, 
#     #           MS1.isotopic.spectrum, MS.MS.spectrum, Average.Mz, Post.curation.result, Fill.., Annotation.tag..VS1.0., RT.matched, 
#     #           m.z.matched, MS.MS.matched, Manually.modified, Total.score:Fragment.presence..))
#   return(df.filtered)
# }


# RearrangeDatasets <- function(df, parameter) {
#   df.gathered <- df %>%
#     tidyr::gather(
#       key = "Replicate.Name",
#       value = eval(parse(text = paste("Area.Value"))),
#       starts_with("X")) %>%
#     select(Replicate.Name, parameter, everything())
#   return(df.gathered)
# }


## TBD VARIABLE RENAMES
# asNumeric <- function(x) as.numeric(as.character(x))
# asCharacter <- function(x) as.character(x)
# 
# factorsNumeric <- function(df) ifelse(grepl('^X', df), 
#                                      modifyList(df, lapply(df[, sapply(df, is.factor)], asNumeric)), 
#                                      modifyList(df, lapply(df[, sapply(df, is.factor)], asCharacter)))

# mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name))

# Do we need this function?
# CheckBlankMatcher <- function(blank.matcher) {
#   # Takes a blank matcher file and separates any multi-value variable
#   # columns into their own row.
#   #
#   # Args:
#   #   blank.matcher: CSV entered by user to match samples with
#   #   appropriate blanks.
#   #
#   # Returns:
#   #   blank.matcher: new CSV with any duplicate values separated
#   #   into their own rows.
#   #
#   blank.matcher <- do.call("rbind", Map(data.frame,
#                                         Blank.Name = strsplit(as.character(blank.matcher$Blank.Name), ","),
#                                         Replicate.Name = (blank.matcher$Replicate.Name))
#   )
#   blank.matcher <- blank.matcher[c(2, 1)]
#   
#   return(blank.matcher)
#}


# IdentifyDuplicates <- function(df) {
#   test <- which(duplicated(df$Compound.Name))
#   duplicates <- as.data.frame(df$Compound.Name[test]) %>%
#     rename(Compound.Name = 1) %>%
#     arrange(Compound.Name)
#   
#   return(duplicates)
# }



# rm(list = ls()[!ls() %in% c("combined", lsf.str())])