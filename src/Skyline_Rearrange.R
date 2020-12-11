## Skyline Rearrange and Compound Name Check

csvFileName <- paste("data_intermediate/", software.pattern, "_combined_", file.pattern, "_", currentDate, ".csv", sep = "")

if (matching.pattern == "pos|neg") {
  
  replace_nonvalues <- function(x) (gsub("#N/A", NA, x))
  
  skyline.HILIC.pos <- skyline.HILIC.pos %>%
    mutate(Column = "HILICpos")
  skyline.HILIC.neg <- skyline.HILIC.neg %>%
    mutate(Column = "HILICneg")
  
  combined.skyline <- skyline.HILIC.pos %>%
    rbind(skyline.HILIC.neg) %>%
    select(Replicate.Name, Precursor.Ion.Name, Retention.Time, Area, Background, Height, Mass.Error.PPM, Column) %>%
    mutate_all(replace_nonvalues) 
  
  skyline.classes.changed <- ChangeClasses(combined.skyline, start.column = 3, end.column = 7) 

  #write.csv(skyline.classes.changed, csvFileName, row.names = FALSE)
  
} else {
  
  skyline.classes.changed <- ChangeClasses(skyline.RP.Cyano, start.column = 4, end.column = 8) 
  skyline.RP.Cyano <- skyline.classes.changed
  rm(skyline.classes.changed)
  
  #write.csv(skyline.RP.Cyano, csvFileName, row.names = FALSE)
  
}

# Check compatibility and replace names -----------------------------------

# Import standards sheet
original.file <- combined.skyline

Ingalls.Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                              stringsAsFactors = FALSE, header = TRUE) %>%
  select(Compound.Name, Compound.Name_old) %>%
  unique()


if ("Precursor.Ion.Name" %in% colnames(original.file)) {
  print("Checking for updated compound names.")
  
  compound.names <- original.file %>%
    select(Precursor.Ion.Name) %>%
    unique() %>%
    rename(Compound.Name_old = Precursor.Ion.Name)
  
  combined.names <- compound.names %>%
    left_join(Ingalls.Standards) %>%
    rename(Precursor.Ion.Name = Compound.Name_old)
  
  combined.skyline <- combined.names %>%
    left_join(original.file) %>%
    mutate(Compound.Name = ifelse(is.na(Compound.Name), Precursor.Ion.Name, Compound.Name)) %>%
    select(-Precursor.Ion.Name) %>%
    rename(Precursor.Ion.Name = Compound.Name)
  print("Compound names updated. Exporting file to data_intermediate.")
  
  #write.csv(combined.skyline, csvFileName, row.names = FALSE)
  
} else {
  
  stop("Your column is incorrectly named.") 
  
}


