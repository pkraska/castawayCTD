#' Castaway CTD to netCDF conversion
#'
#' @param files A character vector of a directory in which castaway CTD CSV
#'   files are stored.
#' @param destination an existing folder to save netCDF files to
#' @param level either \code{raw} or \code{processed} based on the level of
#'   processing completed using the castaway CTD software
#' @param project_lead an optional field to add the project lead's name to the file global attributes
#' @param project_name an optional field to add a project name to the file global attributes
#'
#' @return \code{castaway_convert()} creates netCDF files from castaway CTD CSV
#'   files
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

castaway_convert <-
  function(files,
           destination = "R:/Science/CESD/CESD_DataManagement/data_out/castawayCTD",
           level = "raw",
           project_lead,
           project_name) {
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
      if (cast_header$value[cast_header$key == 'Cast data'] == "Raw") {
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
        
        ncatt_put(
          nc = castaway_nc,
          varid = 0,
          attname = 'project_name',
          attval = project_name
        )
        
        ncatt_put(
          nc = castaway_nc,
          varid = 0,
          attname = 'project_lead',
          attval = project_lead
        )
        
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
        message(paste0(
          "netCDF file ",
          cast_header$value[cast_header$key == "File name"],
          ".nc created in folder ",
          destination,
          "."
        ))
      }
      else {
        if (cast_header$value[cast_header$key == 'Cast data'] == "Processed") {
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
          
          dimDepth <-
            ncdim_def(
              name = 'depth',
              units = 'm',
              longname = 'Depth in meters',
              vals = cast_data$depth_meter
            )
          
          dimTime <-
            ncdim_def(
              name = 'time',
              units = 'seconds since 1970-01-01',
              longname = 'Time',
              vals = as.numeric(ymd_hms(cast_header$value[cast_header$key == 'Cast time (UTC)']))
            )
          
          # Variables
          varLongitude <- ncvar_def(
            name = 'lon',
            units = 'degrees_east',
            dim = list(dimLon, dimLat, dimDepth, dimTime),
            missval = NA,
            longname = "Longitude",
            prec = 'double'
          )
          
          varLatitude <- ncvar_def(
            name = 'lat',
            units = 'degrees_north',
            dim = list(dimLon, dimLat, dimDepth, dimTime),
            missval = NA,
            longname = "Latitude",
            prec = 'double'
          )
          
          varPressure <- ncvar_def(
            name = 'Pres_Z',
            units = 'dbar',
            dim = list(dimLon, dimLat, dimDepth, dimTime),
            missval = NA,
            longname = "Pressure (spatial coordinate) exerted by the water body by profiling pressure sensor and correction to read zero at sea level",
            prec = 'double'
          )

          
          varTemperature <- ncvar_def(
            name = 'WC_temp_CTD',
            units = 'degrees C',
            dim = list(dimLon, dimLat, dimDepth, dimTime),
            missval = NA,
            longname = 'Temperature of the water body by CTD or STD',
            prec = 'double'
          )
          
          varConductivity <- ncvar_def(
            name = 'conductivity',
            units = 'microseimens per cm',
            dim = list(dimLon, dimLat, dimDepth, dimTime),
            missval = NA,
            longname = 'Conductivity (\u00B5S/cm)',
            prec = 'double'
          )
          
          varConductance <- ncvar_def(
            name = 'specific conductance',
            units = 'microseimens per cm',
            dim = list(dimLon, dimLat, dimDepth, dimTime),
            missval = NA,
            longname = 'Specific Conductance (\u00B5S/cm)',
            prec = 'double'
          )
          
          varSalinity <- ncvar_def(
            name = 'salinity',
            units = 'PSS',
            dim = list(dimLon, dimLat, dimDepth, dimTime),
            missval = NA,
            longname = 'Salinity (PSS)',
            prec = 'double'
          )
          
          varSound <- ncvar_def(
            name = 'sound velocity',
            units = 'm/s',
            dim = list(dimLon, dimLat, dimDepth, dimTime),
            missval = NA,
            longname = 'Calculated sound velocity in seawater (m/s)',
            prec = 'double'
          )
          
          varDensity <- ncvar_def(
            name = 'density',
            units = 'kg/m^3',
            dim = list(dimLon, dimLat, dimDepth, dimTime),
            missval = NA,
            longname = 'Calculated Seawater Density (kg/m^3)',
            prec = 'double'
          )
          
          castaway_nc <-
            nc_create(
              paste0(destination, "/", cast_header$value[cast_header$key == 'File name'], ".nc"),
              vars = list(
                varPressure,
                # varDepth,
                varTemperature,
                varConductivity,
                varConductance,
                varSalinity,
                varSound,
                varDensity
              )
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
          
          ncatt_put(
            nc = castaway_nc,
            varid = 0,
            attname = 'project_name',
            attval = project_name
          )
          
          ncatt_put(
            nc = castaway_nc,
            varid = 0,
            attname = 'project_lead',
            attval = project_lead
          )
          ncvar_put(
            castaway_nc,
            varid = varPressure,
            vals = cast_data$pressure_decibar,
            count = c(1, 1, -1, 1)
          )
          
          ncvar_put(
            castaway_nc,
            varid = varTemperature,
            vals = cast_data$temperature_celsius,
            count = c(1, 1, -1, 1)
          )
          
          ncvar_put(
            castaway_nc,
            varid = varConductivity,
            vals = cast_data$conductivity_microsiemens_per_centimeter,
            count = c(1, 1, -1, 1)
          )
          
          ncvar_put(
            castaway_nc,
            varid = varConductance,
            vals = cast_data$specific_conductance_microsiemens_per_centimeter,
            count = c(1, 1, -1, 1)
          )
          
          ncvar_put(
            castaway_nc,
            varid = varSalinity,
            vals = cast_data$salinity_practical_salinity_scale,
            count = c(1, 1, -1, 1)
          )
          
          ncvar_put(
            castaway_nc,
            varid = varSound,
            vals = cast_data$sound_velocity_meters_per_second,
            count = c(1, 1, -1, 1)
          )
          
          ncvar_put(
            castaway_nc,
            varid = varDensity,
            vals = cast_data$density_kilograms_per_cubic_meter,
            count = c(1, 1, -1, 1)
          )
          # close the nc file to ensure no data is lost
          nc_close(castaway_nc)
          message(paste0(
            "netCDF file ",
            cast_header$value[cast_header$key == "File name"],
            ".nc created in folder ",
            destination,
            "."
          ))
        }
      }
    }
  }
