\name{castaway_convert}
\alias{castaway convert}
\title{Castaway CTD Convert function}
\usage{
castaway_convert(directory)
}
\description{
converts comma separated value (.CSV) castaway CTD data into a tibble of the header information in key-value format, and a tibble of the cast data.
}
\examples{
castaway_convert(directory = "data/")
}
