## Skyline Rearrange and Compound Name Check

## THIS FUNCTION IS NOT FINAL. IT WILL NEED TO BE ADJUSTED
## TO PROPERLY HANDLE DIFFERENT INSTRUMENT INPUTS.

# CheckCompoundNames <- function (df) {
#   
#   Ingalls.Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
#                                 stringsAsFactors = FALSE, header = TRUE) %>%
#     select(Compound.Name, Compound.Name_old) %>%
#     unique()
#   
#   if ("Precursor.Ion.Name" %in% colnames(df)) {
#     print("Checking for updated compound names.")
#     
#     compound.names <- df %>%
#       select(Precursor.Ion.Name) %>%
#       unique() %>%
#       rename(Compound.Name_old = Precursor.Ion.Name)
#     
#     combined.names <- compound.names %>%
#       left_join(Ingalls.Standards) %>%
#       rename(Precursor.Ion.Name = Compound.Name_old)
#     
#     combined.final <- combined.names %>%
#       left_join(df) %>%
#       mutate(Compound.Name = ifelse(is.na(Compound.Name), Precursor.Ion.Name, Compound.Name)) %>%
#       select(-Precursor.Ion.Name) %>%
#       rename(Precursor.Ion.Name = Compound.Name)
#     print("Compound names updated. Exporting file to data_intermediate.")
#     
#   } else {
#     
#     stop("Your column is incorrectly named.") 
#     
#   }
#   
#   return(combined.final)
# }

# Remove any illegal values
replace_nonvalues <- function(x) (gsub("#N/A", NA, x))

csvFileName <- paste("data_intermediate/", software.pattern, "_combined_", file.pattern, "_", currentDate, ".csv", sep = "")

if (matching.pattern == "pos|neg") {
  
  skyline.HILIC.pos <- skyline.HILIC.pos %>%
    mutate(Column = "HILICpos")
  skyline.HILIC.neg <- skyline.HILIC.neg %>%
    mutate(Column = "HILICneg")
  
  combined.skyline <- skyline.HILIC.pos %>%
    rbind(skyline.HILIC.neg) %>%
    select(Replicate.Name, Precursor.Ion.Name, Retention.Time, Area, Background, Height, Mass.Error.PPM, Column) %>%
    mutate_all(replace_nonvalues) 
  
  # Change variable classes
  skyline.classes.changed <- ChangeClasses(combined.skyline, start.column = 3, end.column = 7) 
  
  # Check and update any old compound names
  #skyline.names.checked <- CheckCompoundNames(skyline.classes.changed)
  
  #write.csv(skyline.names.checked, csvFileName, row.names = FALSE)
  write.csv(skyline.classes.changed, csvFileName, row.names = FALSE)
  
} else {
  
  skyline.RP.Cyano2 <- skyline.RP.Cyano %>%
    mutate_all(replace_nonvalues)
  
  # Change variable classes
  skyline.classes.changed <- ChangeClasses(skyline.RP.Cyano2, start.column = 4, end.column = 8) 
  
  # Check and update any old compound names
  skyline.names.checked <- CheckCompoundNames(skyline.classes.changed)
  
  skyline.RP.Cyano <- skyline.names.checked
  rm(skyline.classes.changed, skyline.names.checked)
  
  write.csv(skyline.RP.Cyano, csvFileName, row.names = FALSE)
  
}







