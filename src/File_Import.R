
# Import all files --------------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_raw", pattern = file.pattern))

for (i in filenames) {
  filepath <- file.path("data_raw", paste(i,".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = TRUE))
}
