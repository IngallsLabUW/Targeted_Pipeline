# Duplicates testing

HILICS.duplicates <- IdentifyDuplicates(QCd.data)

if ("Column" %in% colnames(QCd.data)) {
  duplicates.testing <- QCd.data %>%
    filter(Metabolite.Name %in% HILICS.duplicates["Metabolite.Name"]) %>%
    group_by(Metabolite.Name, Column) %>%
    summarize(mean(Area.with.QC, na.rm = TRUE))
  
} else {
  print("Non-HILIC data: no duplicates to detect.")
}

if (exists("duplicates.testing") == TRUE) {
  QCd.data <- QCd.data %>%
    filter(!(Metabolite.Name %in% HILICS.duplicates$Metabolite.Name & Column == "HILICNeg"))
} else {
  print("Non-HILIC data: no duplicates to remove.")
}