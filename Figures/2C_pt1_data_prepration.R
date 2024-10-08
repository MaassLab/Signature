#_____Load required packages_____________________________________________________
library(tidyr)
library(dplyr)
library(data.table)

#_____Read in arguments_________________________________________________________
args = commandArgs(trailingOnly = TRUE)

data_path <- args[1]
outpath <- args[2]

#____split column ID____________________________________________________________
data <- read.table(data_path), header= T)
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
data$ID <- sub("B", "\\.B", as.character(data$ID))
data <- data %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
head(data)
data<- data[,-c(1,4,7)]
data

#make 62 files for each data set
for (col in names(data)[5:ncol(data)]) {


  if (file.exists(paste0(cells_pathway, col, ".txt"))) {
    next
  }

  new_table <- data[c( "chrA", "st1", "chrB", "st2", col)]

  file_name <- paste0(col, ".txt")
  write.table(new_table, file = sprintf("%s/%s", outpath, file_name), sep = "\t", row.names = FALSE)
  
  
}

#______________________remove rows with NA in the last column____________________

na.omit_last <- function(data_list) {
  last_col <- ncol(data_list)
  na_rows <- is.na(data_list[, last_col])
  data_list[!na_rows, ]
}

for (filename in list.files(path = outpath ,pattern = "\\.txt$")) {
  
  data_list <- read.table(filename, header = TRUE, sep = "\t")
  
  data_list <- na.omit_last(data_list)
  
  write.table(data_list, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
}


