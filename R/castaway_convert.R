#' Castaway CTD to netCDF conversion
#'
#' @param files A character vector of a directory in which castaway CTD CSV
#'   files are stored. Must be in "Raw" format from the Castaway CTD Software
#'   (only variable for \code{Time (Seconds)},	\code{Pressure
#'   (Decibar)},	\code{Temperature (Celsius)}, and \code{Conductivity
#'   (MicroSiemens per Centimeter)})
#' @param destination an existing folder to save netCDF files to
#'
#' @return \code{castaway_convert()} creates netCDF files from castaway CTD CSV files
#' @export
#'
#' @examples castaway_convert("data/")
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import readr
#' @import ncdf4
#' @import lubridate
#'

castaway_convert <- function(files, destination = "R:/Science/CESD/CESD_DataManagement/data_out/castawayCTD") {
  for (i in files) {
    # read data file
    data <- readLines(i) %>%
      as_tibble() %>%
      mutate(header = grepl("% ", value))


    cast_header <- data %>%
      filter(header == TRUE) %>%
      select(value) %>%
      separate(
        value,
        into = c("key", "value"),
        sep = ",",
        fill = 'right'
      ) %>%
      mutate(
        key = str_replace(
          string = .$key,
          pattern = "% ",
          replacement = ""
        ),
        key = str_trim(string = key, side = "both"),
        value = str_trim(string = value, side = "both")
      ) %>%
      filter(key != "")

    castaway_data_cols <- data %>%
      slice(sum(data$header) + 1) %>%
      select(value) %>%
      str_replace_all(pattern = " ", replacement = "_") %>%
      str_replace_all(pattern = "\\(", replacement = "") %>%
      str_replace_all(pattern = "\\)", replacement = "") %>%
      str_to_lower() %>%
      str_split(pattern = ",") %>%
      unlist()

    cast_data <- data %>%
      filter(header == FALSE) %>%
      select(value) %>%
      slice(-1) %>%
      separate(
        col = value,
        into = castaway_data_cols,
        sep = ",",
        convert = TRUE
      )

    # Dimensions
    # Castaway data will be treated as a single dimension profile datatype

    dimLon <-
      ncdim_def(
        name = 'lon',
        units = 'degrees_east',
        longname = 'Longitude',
        vals = as.numeric(cast_header$value[cast_header$key == 'Start longitude'])
      )

    dimLat <-
      ncdim_def(
        name = 'lat',
        units = 'degrees_north',
        longname = 'Latitude',
        vals = as.numeric(cast_header$value[cast_header$key == 'Start latitude'])
      )

    dimTime <-
      ncdim_def(
        name = 'time',
        units = 'seconds since 1970-01-01',
        longname = 'Time',
        vals = as.numeric(ymd_hms(cast_header$value[cast_header$key == 'Cast time (UTC)']) + cast_data$time_seconds)
      )

    # dimnchar <- ncdim_def("nchar", "", 1:25, create_dimvar = FALSE)

    # Variables
    varLongitude <- ncvar_def(
      name = 'lon',
      units = 'degrees_east',
      dim = list(dimLon, dimLat, dimTime),
      missval = NA,
      longname = "Longitude",
      prec = 'double'
    )

    varLatitude <- ncvar_def(
      name = 'lat',
      units = 'degrees_north',
      dim = list(dimLon, dimLat, dimTime),
      missval = NA,
      longname = "Latitude",
      prec = 'double'
    )

    varPressure <- ncvar_def(
      name = 'Pres_Z',
      units = 'dbar',
      dim = list(dimLon, dimLat, dimTime),
      missval = NA,
      longname = "Pressure (spatial coordinate) exerted by the water body by profiling pressure sensor and correction to read zero at sea level",
      prec = 'double'
    )

    varTemperature <- ncvar_def(
      name = 'WC_temp_CTD',
      units = 'degrees C',
      dim = list(dimLon, dimLat, dimTime),
      missval = NA,
      longname = 'Temperature of the water body by CTD or STD',
      prec = 'double'
    )

    varConductivity <- ncvar_def(
      name = 'conductivity',
      units = 'microseimens per cm',
      dim = list(dimLon, dimLat, dimTime),
      missval = NA,
      longname = 'Conductivity (\u00B5S/cm)',
      prec = 'double'
    )

    castaway_nc <-
      nc_create(
        paste0(destination, "/", cast_header$value[cast_header$key == 'File name'], ".nc"),
        vars = list(varPressure,
                    varTemperature,
                    varConductivity)
      )

    # Include all header information in the netCDF file as attribute data
    for (i in cast_header$key) {
      ncatt_put(
        nc = castaway_nc,
        varid = 0,
        attname = i,
        attval = cast_header$value[cast_header$key == i]
      )
    }

    ncvar_put(
      castaway_nc,
      varid = varPressure,
      vals = cast_data$pressure_decibar,
      count = c(1, 1,-1)
    )

    ncvar_put(
      castaway_nc,
      varid = varTemperature,
      vals = cast_data$temperature_celsius,
      count = c(1, 1, -1)
    )

    ncvar_put(
      castaway_nc,
      varid = varConductivity,
      vals = cast_data$conductivity_microsiemens_per_centimeter,
      count = c(1, 1, -1)
    )

    # close the nc file to ensure no data is lost
    nc_close(castaway_nc)
    message(paste0("netCDF file ", cast_header$value[cast_header$key == "File name"], ".nc created."))
  }
}
