## Skyline Rearrange

csvFileName <- paste("data_intermediate/", software.pattern, "_combined_", file.pattern, "_", currentDate, ".csv", sep = "")

if (file.pattern == "HILIC") {
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

  write.csv(skyline.classes.changed, csvFileName, row.names = FALSE)
  
} else {
  
  skyline.classes.changed <- ChangeClasses(skyline.RP.Cyano, start.column = 4, end.column = 8) 
  skyline.RP.Cyano <- skyline.classes.changed
  rm(skyline.classes.changed)
  
  write.csv(skyline.RP.Cyano, csvFileName, row.names = FALSE)
  
}


