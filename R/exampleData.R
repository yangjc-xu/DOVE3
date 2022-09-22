#' Toy Dataset For Illustration
#'
#' This data set is provided for the purposes of illustrating the use of
#' the software. It was a simulated data set that mimics the surveillance data of
#' vaccination and infection.
#'
#'
#' @usage data(exampleData)
#'
#' @format idoveData is a data.frame containing 8,000 records.
#'   The data.frame contains 11 columns,
#'   \describe{
#'   \item{subject.id}{The subject ID of each record.}
#'   \item{event.time}{The event time in days.}
#'   \item{censor.time}{The censoring time for each subject in days.}
#'   \item{entry.time}{The entry time in days.}
#'   \item{Vtime}{The time of vaccination/booster in days.}
#'   \item{Vtype}{The type of vaccination/booster.}
#'   \item{infection.time}{The time of infection in days.}
#'   \item{infection.type}{The type of infection.}
#'   \item{age}{A categorical variable of age (<18, 18-34, 35-49, 50-64, >=65).}
#'   \item{gender}{A binary indicator of gender (male/female).}
#'   \item{priority}{A composite baseline risk score taking values 1-3.}
#'   }
#'
#' @name exampleData
#' @keywords datasets
NULL
