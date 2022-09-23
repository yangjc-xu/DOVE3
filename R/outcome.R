#' Specify Outcome Variables
#'
#' This function is used in the model statement of dove3() to specify
#'  the subject id, entry time, event time, and censoring time.
#'
#' @param subject_id The variable for the subject ID that is used to identify
#'   records from the same subject, must be an integer between 1 and the total
#'   number of subjects.
#'
#' @param entry_time The variable for the time when the subject enters the risk set.
#'   For example, entry_time is set to be 0 for all subjects if effectiveness of
#'   primary vaccine series is of interest, and is set to be the time of
#'   primary vaccine series if effectiveness of boosters is of interest.
#'
#' @param event_time The variable for the time when the subject has an event.
#'
#' @param censor_time The variable for the time when the subject is censored.
#'   For example, subjects are censored when they receive boosters if
#'   effectiveness of primary vaccine series is of interest.
#'   Likewise, subjects are censored when they receive second boosters if
#'   effectiveness of first boosters is of interest.
#'
#' @returns This function is intended to be used only in the model statement
#'  of dove3(). The result, a matrix, is used internally.
#'
#' @name outcome
#' @rdname outcome
#' @export

outcome = function(subject_id, entry_time, event_time, censor_time) {

  ### subject_id

  # must be provided as a vector of integers.

  if (missing(x = subject_id)) {
    stop("must provide a subject_id argument",
         call. = FALSE)
  }

  if (!is.integer(x = subject_id)) {
    stop ("subject_id must be integer", call. = FALSE)
  }

  if(min(subject_id, na.rm = T) != 1 |
     max(subject_id, na.rm = T) != length(unique(na.omit(subject_id)))) {
    stop ("subject_id must range from 1 to total number of subjects", call. = FALSE)
  }

  ### entry_time

  # must be provided as a numeric vector.

  if (missing(x = entry_time)) {
    stop("must provide a entry_time argument",
         call. = FALSE)
  }

  if (!is.numeric(x = entry_time)) {
    stop ("entry_time is not numeric", call. = FALSE)
  }

  ### event_time

  # must be provided as a numeric vector.

  if (missing(x = event_time)) {
    stop("must provide a event_time argument",
         call. = FALSE)
  }

  if (!is.numeric(x = event_time)) {
    stop ("event_time is not numeric", call. = FALSE)
  }

  # if NA, set to be censor_time + 1 day
  ind = which(is.na(event_time))
  event_time[ind] = censor_time[ind]+1

  ### censor_time

  # must be provided as a numeric vector.

  if (missing(x = censor_time)) {
    stop("must provide a censor_time argument",
         call. = FALSE)
  }

  if (!is.numeric(x = censor_time)) {
    stop ("censor_time is not numeric", call. = FALSE)
  }

  ### check if all inputs are of the same length

  if (length(x = subject_id) != length(x = entry_time) |
      length(x = subject_id) != length(x = event_time) |
      length(x = subject_id) != length(x = censor_time)) {
    stop("all inputs must be of same length", call. = FALSE)
  }

  ### check if entry_time <= censor_time

  tst <- {entry_time > censor_time}

  if (any(tst, na.rm = TRUE)) {
    message(sum(tst, na.rm = TRUE),
            " don't satisfy entry_time <= censor_time; cases set to NA")
    entry_time[tst] <- NA
  }

  dm <- cbind(subject_id, entry_time, event_time, censor_time)

  cname <- c("subject_id", "entry_time", "event_time", "censor_time")
  dimnames(x = dm) <- list(NULL, cname)

  return( dm )
}
