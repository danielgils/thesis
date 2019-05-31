# Script to download files with the number of downloads on DOE packages
# Files are stored in
# http://cran-logs.rstudio.com/

# Here's an easy way to get all the URLs in R
start <- as.Date('2018-01-01')
today <- as.Date('2018-10-22')

all_days <- seq(start, today, by = 'day')

year <- as.POSIXlt(all_days)$year + 1900
urls <- paste0('http://cran-logs.rstudio.com/', year, '/', all_days, '.csv.gz')
dest <- paste0('Data/', year, '/', all_days, '.csv.gz') #I had to create 2018 folder first
# You can then use download.file to download into a directory.

for (i in 1:length(urls)) {
  download.file(urls[i], dest[i], "internal", quiet = FALSE, mode = "w")
}

# If you only want to download the files you don't have, try:
missing_days <- setdiff(all_days, tools::file_path_sans_ext(dir(), TRUE))


              
