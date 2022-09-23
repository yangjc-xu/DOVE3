#' Specify Exposure Variables
#'
#' This function is used in the model statement of dove3() to specify
#'  the vaccination time, vaccine type, infection time, and infection type.
#'
#' @param vaccine_time The variable for the time when the subject is vaccinated,
#'   with NA, Inf, or an arbitrary value that is larger than the study end time
#'   if the subject is never vaccinated during the study period.
#'
#' @param vaccine_type The variable for the vaccine type the subject receives,
#'   must be an integer between 1 and the total number of vaccine types,
#'   with NA or an arbitrary value within this range if the subject is never
#'   vaccinated during the study period.
#'
#' @param infection_time The variable for the time when the subject is infected,
#'   with NA, Inf, or an arbitrary value that is larger than the study end time if the
#'   subject is never infected during the study period.
#'
#' @param infection_type The variable for the dominant variant at the time when
#'   the subject is infected, must be an integer between 1 and the total
#'   number of variants under investigation, with NA or an arbitrary value within
#'   this range if the subject is never infected during the study period.
#'
#' @returns This function is intended to be used only in the model statement
#'  of dove3(). The result, a matrix, is used internally.
#'
#' @name exposure
#' @rdname exposure
#' @export

exposure = function(vaccine_time, vaccine_type, infection_time, infection_type) {

  ### vaccine_time

  # must be provided as a numeric vector.

  if (missing(x = vaccine_time)) {
    stop("must provide a vaccine_time argument",
         call. = FALSE)
  }

  if (!is.numeric(x = vaccine_time)) {
    stop ("vaccine_time is not numeric", call. = FALSE)
  }

  ind = which(is.na(vaccine_time))
  vaccine_time[ind] = 9999

  ### vaccine_type

  # must be provided as a vector of integers, with unique values 1,2,...K.

  if (missing(x = vaccine_type)) {
    stop("must provide a vaccine_type argument",
         call. = FALSE)
  }

  if (!is.numeric(x = vaccine_type)) {
    stop ("vaccine_type must be integer", call. = FALSE)
  }

  if(min(vaccine_type, na.rm = T) != 1 |
     max(vaccine_type, na.rm = T) != length(unique(na.omit(vaccine_type)))) {
    stop ("vaccine_type must range from 1 to total number of vaccine types", call. = FALSE)
  }

  ind = which(is.na(vaccine_type))
  vaccine_type[ind] = 1

  ### infection_time

  # must be provided as a numeric vector.

  if (missing(x = infection_time)) {
    stop("must provide a infection_time argument",
         call. = FALSE)
  }

  if (!is.numeric(x = infection_time)) {
    stop ("infection_time is not numeric", call. = FALSE)
  }

  ind = which(is.na(infection_time))
  infection_time[ind] = 9999

  ### infection_type

  # must be provided as a vector of integers, with unique values 1,2,....

  if (missing(x = infection_type)) {
    stop("must provide a infection_type argument",
         call. = FALSE)
  }

  if (!is.numeric(x = infection_type)) {
    stop ("infection_type must be integer", call. = FALSE)
  }

  if(min(infection_type, na.rm = T) != 1 |
     max(infection_type, na.rm = T) != length(unique(na.omit(infection_type)))) {
    stop ("infection_type must range from 1 to total number of infection types", call. = FALSE)
  }

  ind = which(is.na(infection_type))
  infection_type[ind] = 1

  ### check if all inputs are of the same length

  if (length(x = vaccine_time) != length(x = vaccine_type) |
      length(x = vaccine_time) != length(x = infection_time) |
      length(x = vaccine_time) != length(x = infection_type)) {
    stop("all inputs must be of same length", call. = FALSE)
  }

  dm <- cbind(vaccine_time, vaccine_type, infection_time, infection_type)

  cname <- c("vaccine_time", "vaccine_type", "infection_time", "infection_type")
  dimnames(x = dm) <- list(NULL, cname)

  return( dm )
}
