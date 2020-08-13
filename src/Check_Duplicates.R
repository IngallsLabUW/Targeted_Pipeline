# Duplicates testing

HILICS.duplicates <- IdentifyDuplicates(QCd.data)

if ("Column" %in% colnames(QCd.data)) {
  duplicates.testing <- QCd.data %>%
    filter(Metabolite.Name %in% HILICS.duplicates$Metabolite.Name) %>%
    group_by(Metabolite.Name, Column) %>%
    mutate(Means = mean(Area.with.QC, na.rm = TRUE)) %>%
    mutate(Std.Devs = sd(Area.with.QC, na.rm = TRUE)) %>%
    ungroup() %>%
    select(Metabolite.Name, Column, Means, Std.Devs) %>%
    unique()
  

} else {
  print("Non-HILIC data: no duplicates to detect.")
}

