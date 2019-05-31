#' Analysis of first semester of downloaded packages
#' @author Daniel Gil
#' October 2018

#' Here's an easy way to get all the URLs in R
start <- as.Date('2018-01-01')
today <- as.Date('2018-10-22')

all_days <- seq(start, today, by = 'day')

year <- as.POSIXlt(all_days)$year + 1900
# urls <- paste0('http://cran-logs.rstudio.com/', year, '/', all_days, '.csv.gz')
dest <- paste0('Data/', year, '/', all_days, '.csv.gz') #I had to create 2018 folder

data <- read.table(gzfile(dest[1]))

# All files were corrupted, so i didnt continue
                   