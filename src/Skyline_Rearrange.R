## Skyline Rearrange and Compound Name Check

# Define custom file name for export
csvFileName <- paste("data_intermediate/", software.pattern, "_combined_", file.pattern, "_", currentDate, ".csv", sep = "")

# Function to remove syntactically incorrect values usually produced by Skyline
replace_nonvalues <- function(x) (gsub("#N/A", NA, x))

# Replace original compound names with updated Standards names
update_compound_names <- function(df) {
  names.changed <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
                                stringsAsFactors = FALSE, header = TRUE) %>%
    select(Compound_Name, Compound_Name_Original) %>%
    unique() %>%
    full_join(df %>% rename(Compound_Name_Original = Precursor.Ion.Name)) %>%
    filter(Compound_Name_Original %in% df$Precursor.Ion.Name) %>%
    select(Precursor.Ion.Name = Compound_Name, everything(), -Compound_Name_Original)
}

# Identify positive and negative HILIC runs
if (runtype.pattern == "pos|neg") {
  
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
  
  # Fix old compound names
  skyline.names.updated <- update_compound_names(skyline.classes.changed)
  
  # Export rearranged dataframe
  write.csv(skyline.names.updated, csvFileName, row.names = FALSE)
  
} else {
  
  skyline.reversephase.nonvalues <- skyline.reversephase %>%
    mutate_all(replace_nonvalues)
  
  # Change variable classes
  skyline.classes.changed <- ChangeClasses(skyline.reversephase.nonvalues, start.column = 4, end.column = 8) 
  
  # Fix old compound naames
  skyline.names.updated <- update_compound_names(skyline.classes.changed)
  
  # Export rearranged dataframe
  write.csv(skyline.names.updated, csvFileName, row.names = FALSE)
  
}
