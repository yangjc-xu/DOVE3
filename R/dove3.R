#' Effectiveness of Vaccination and Prior Infection
#'
#' Estimates the potentially waning long-term effectiveness of vaccination and
#'  prior infection against COVID-19 outcomes in observational studies.
#'  Effects of all exposures (i.e., vaccination and prior infection) are estimated
#'  simultaneously under a single Cox regression model, allowing the outcome
#'  of interest to be a recurrent event.
#'
#' The log intensity or hazard ratio for an exposure is a piecewise
#'  linear function of time since the exposure.
#'  Specific constraint on the first and last piece is allowed.
#'  In addition, vaccine effects can be evaluated among previously uninfected
#'  and infected subjects separately.
#'
#' The information required for an analysis is
#'   \describe{
#'     \item{Subject ID:}{Number that identifies subjects.}
#'     \item{Entry Time:}{Calendar time when the subject enters the risk set.}
#'     \item{Event Time:}{Calendar time when the subject experiences
#'       the clinical outcome of interest (e.g., SARS-CoV-2 Infection,
#'       hospitalization, or death caused by infection), with NA, Inf, or an
#'       arbitrary value greater than the censoring time if the subject does not
#'       experience an event.}
#'     \item{Censoring Time:}{Calendar time when the subject moves out of the risk set.}
#'     \item{Vaccination Time:}{Calendar time when vaccination takes place,
#'       with NA, Inf, or an arbitrary value greater than the study end time
#'       if the subject is not vaccinated during the study period.}
#'     \item{Vaccine Type:}{Categorical variable indicating which vaccine type
#'       the subject receives. Must be an integer between 1 and the total number
#'       of vaccine types, with NA or an arbitrary value within this range if the
#'       subject is not vaccinated during the study period.}
#'     \item{Infection Time:}{Calendar time when the subject is infected,
#'       with NA, Inf, or an arbitrary value that is larger than the study end
#'       time if the subject is not infected during the study period.}
#'    \item{Infection Type:}{Categorical variable indicating which dominant variant
#'       the infection is associated with. Must be an integer between 1 and the total
#'       number of variants under investigation, with NA or an arbitrary value within
#'       this range if the subject is not infected during the study period.}
#'     \item{Covariates:}{Baseline covariates (e.g., age, gender, ethnicity).}
#'    }
#'
#' Note that a subject can have multiple rows and each row
#' corresponds to an infection.
#' The event_time is the time of the clinical outcome of interest caused by the
#' infection in the same row, with NA, Inf, or an
#' arbitrary value greater than the censor_time if that infection does not cause
#' the outcome.
#' If a subject never experiences the clinical outcome of interest during the
#' entire study period, this subject will have only one record, with the event_time
#' being NA, Inf, or an arbitrary value greater than the censor_time, and with the
#' infection_time being NA, Inf, or any arbitrary value.
#'
#' All the time variables are measured from the same time origin
#' and are specified in units of days.
#' For each individual, the entry_time and censor_time must satisfy
#' entry_time \eqn{\le} censor_time.
#' If entry_time > censor_time, the case will be
#' removed from the analysis and a message will be generated.
#' The definitions of entry_time and censor_time depend on specific analyses.
#' For example, entry_time is set to be 0 for all subjects if effectiveness of primary
#' vaccine series is of interest, and is set to be the completion time of the
#' primary vaccine series if effectiveness of boosters is of interest.
#' Thus, the entry_time and censor_time are not necessarily the beginning and
#' end of the follow-up period for a subject, and it is possible that
#' event_time is smaller than entry_time or larger than censor_time.
#'
#' All the categorical variables (i.e., subject_id, vaccine_type, infection_type)
#' must be labelled as 1, 2, \eqn{\dots}, \eqn{K}, where \eqn{K} is the total
#' number of categories.
#'
#'
#' The general structure of the formula input is
#'   \preformatted{
#'   outcome(subject_id, entry_time, event_time, censor_time) ~ covariates +
#'     exposure(vaccine_time, vaccine_type, infection_time, infection_type)
#'   }
#'
#' The response variable contains the information of event time and event status
#' for each infection and must be specified through function 'outcome()'.
#' All the variables related to vaccination and infection are specified
#' through function 'exposure()'.
#' The covariates can be either numerical or categorical.
#' If categorical covariates are provided, all other categories are compared to
#' the first category.
#'
#' @rdname dove3
#' @name dove3
#'
#' @param formula A formula object, with the response on the left-hand side of a
#'   '~' operator, and the covariates and exposure() function on the right-hand side.
#'   The response contains the outcome information and
#'   must be specified through function 'outcome()'.
#'   The exposure() function must be used to specify the information of each exposure.
#'   See ?exposure, ?outcome, and Details for further information.
#'
#' @param data A data.frame object. The data.frame in which to interpret the
#'   variable names in formula. Must contain the subject ID,
#'   the entry and censoring times, the event time,
#'   the vaccination time and vaccine type, the infection time and
#'   infection type (determined by the dominant variant),
#'   and the covariates. See Details.
#'
#' @param vaccine_infection_interaction A logical object. If TRUE, interaction term
#'   between vaccination and prior infection is included in the model. That is,
#'   vaccine effects are allowed to be different between previously uninfected subjects
#'   and previously infected subjects. If FALSE (default), there is no interaction
#'   term in the model, i.e., average vaccine effects among all subjects are estimated,
#'   regardless of the prior infection status.
#'
#' @param vaccine_knots A list object or NULL. The \eqn{k}th element specifies
#'   the knots of the piecewise linear function for the log rate/hazard ratio
#'   of the \eqn{k}th vaccine type. If NULL, knots will be placed
#'   at every month by default.
#'   This input is ignored when vaccine_infection_interaction = TRUE.
#'
#' @param vaccine_uninfected_knots A list object or NULL.
#'   The \eqn{k}th element specifies the knots of the piecewise linear function
#'   for the log rate/hazard ratio of the \eqn{k}th vaccine type
#'   given no prior infection before vaccination.
#'   If NULL, knots will be placed at every month by default.
#'   This input is ignored when vaccine_infection_interaction = FALSE.
#'
#' @param vaccine_infected_knots A list object or NULL.
#'   The \eqn{k}th element specifies the knots of the piecewise linear function
#'   for the log rate/hazard ratio of the \eqn{k}th vaccine type
#'   given at least one prior infection before vaccination.
#'   If NULL, knots will be placed at every other month by default.
#'   This input is ignored when vaccine_infection_interaction = FALSE.
#'
#' @param prior_infection_knots A list object or NULL.
#'   The \eqn{k}th element specifies the knots of the piecewise linear function
#'   for the log rate/hazard ratio of the \eqn{k}th infection type.
#'   The first knot should be placed at reinfection_cutoff.
#'   If NULL, a default set of knots will be used, with the first knot placed
#'   at reinfection_cutoff and the subsequent knots placed at every three months.
#'
#' @param related_vaccine_types A list object or NULL.
#'   Each element of the list takes the form c(i, j) and imposes a constraint that the
#'   slope of the first piece of the piecewise linear function is the same between
#'   the ith and jth vaccine types.
#'   This input is intended to account for the fact that the 1-dose and 2-dose
#'   regimens of an mRNA vaccine should have the same effectiveness from
#'   the receipt of the first dose until the receipt of the second dose.
#'   If NULL (default), there is no constraint on the piecewise linear functions
#'   for different vaccine types during the estimation.
#'
#' @param last_piece_constant A logical object. If FALSE (default), effectiveness
#'   is allowed to vary after the last knot. If TRUE, effectiveness is assumed to be
#'   constant after the last knot.
#'
#' @param reinfection_cutoff A positive scalar object. Positive results at two time
#'   points separated by a gap larger than this number (in days) are considered as
#'   two different infections.
#'
#' @param plots A logical object. If TRUE (default), plots of the estimated
#'   effectiveness of each vaccine type and prior infection type against the
#'   clinical outcome of interest and their 95\% confidence intervals will be
#'   automatically generated. If FALSE, plots will not be generated.
#'
#'
#' @returns An list object with elements
#'
#'   \item{covariates}{A matrix containing the estimated (log) hazard ratio of
#'     each covariate, together with the estimated standard error, the 95\%
#'     confidence interval, and the two-sided p-value for testing no covariate
#'     effect.}
#'
#'   \item{effectiveness}{A list of matrices, one for each type of exposure 
#'     (vaccination comes first and prior infection next, both in order of type).
#'     Each matrix contains the daily effectiveness estimates in reducing the rate
#'     or hazard of the clinical outcome of interest, together with the
#'     standard errors and the 95\% confidence intervals.}
#'
#'   \item{plots}{A list of plot objects returned by ggplot(), one for each type 
#'     of exposure (vaccination comes first and prior infection next, 
#'     both in order of type).}
#'
#' @references Lin D, Gu Y, Xu Y, et al. Association of Primary and Booster
#'   Vaccination and Prior Infection With SARS-CoV-2 Infection and Severe
#'   COVID-19 Outcomes. JAMA. Published online September 26, 2022.
#'   doi:10.1001/jama.2022.17876
#'
#' @export
#' @import methods Rcpp
#'
#' @useDynLib DOVE3
#'
#' @include CoxReg.R exposure.R outcome.R
#' @importFrom stats model.response update.formula complete.cases
#'
#' @examples
#' data(exampleData)
#'
#' # specify the knots for each exposure
#' vaccine.knots = list("vac.type1" = c(30),
#'                      "vac.type2" = c(30,60),
#'                      "vac.type3" = c(30,60))
#' prior.infection.knots = list("inf.type1" = c(14, 120),
#'                              "inf.type2" = c(14,120))
#'
#' # Fit the simple model without interaction or related vaccine types
#' fit1 <- dove3(formula = outcome(subject.id, entry.time, event.time, censor.time) ~
#'                 age + gender + priority +
#'                 exposure(Vtime, Vtype, infection.time, infection.type),
#'       data = exampleData,
#'       vaccine_knots = vaccine.knots,
#'       prior_infection_knots = prior.infection.knots)
#'
#' # Specify the knots for vaccination without and with prior infection
#' vaccine.uninfected.knots = list("vac.noinf.type1" = c(30),
#'                                 "vac.noinf.type2" = c(30,60),
#'                                 "vac.noinf.type3" = c(30,60))
#'
#' vaccine.infected.knots = list("vac.noinf.type1" = c(30),
#'                                 "vac.noinf.type2" = c(60),
#'                                 "vac.noinf.type3" = c(60))
#'
#' # Fit the model with interaction between vaccination and prior infection status,
#' # and impose a constraint on the first pieces of the first two vaccine types.
#' fit2 <- dove3(formula = outcome(subject.id, entry.time, event.time, censor.time) ~
#'                 age + gender + priority +
#'                 exposure(Vtime, Vtype, infection.time, infection.type),
#'       data = exampleData,
#'       vaccine_infection_interaction = TRUE,
#'       vaccine_uninfected_knots = vaccine.uninfected.knots,
#'       vaccine_infected_knots = vaccine.infected.knots,
#'       prior_infection_knots = prior.infection.knots,
#'       related_vaccine_types = list(c(1,2)))
#'

dove3 <- function(formula,
                  data,
                  vaccine_infection_interaction = FALSE,
                  vaccine_knots = NULL,
                  vaccine_uninfected_knots = NULL,
                  vaccine_infected_knots = NULL,
                  prior_infection_knots = NULL,
                  related_vaccine_types = NULL,
                  last_piece_constant = FALSE,
                  reinfection_cutoff = 14,
                  plots = TRUE) {

  if (missing(x = formula)) {
    stop("a formula argument must be provided", call. = FALSE)
  }

  if (missing(x = data)) {
    stop("a data argument must be provided", call. = FALSE)
  }

  # reset options to allow for keeping na values
  opt <- options()
  options(na.action = 'na.pass')
  on.exit(options(opt))

  # add intercept from model if not provided to ensure that factors are handled properly
  if (attr(x = stats::terms(x = formula), which = "intercept") == 0L) {
    formula = update.formula(old = formula, new = .~. +1)
  }

  # try to obtain the model.frame
  mf <- tryCatch(expr = stats::model.frame(formula = formula, data = data),
                 error = function(e) {
                   message("unable to obtain model.frame")
                   stop(e$message, call. = FALSE)
                 })

  # extract covariates
  X <- suppressMessages(stats::model.matrix(object = mf, data = data))
  # remove intercept
  int <- attr(x = X, which = "assign") != 0L
  X <- X[,int, drop = FALSE]

  # identify the columns that correspond to the returns by exposure()

  lbl <- attr(x = stats::terms(x = formula), which = "term.labels")

  lbl1 <- paste0(lbl, "vaccine_time")
  lbl2 <- paste0(lbl, "vaccine_type")
  lbl3 <- paste0(lbl, "infection_time")
  lbl4 <- paste0(lbl, "infection_type")


  if (!any(lbl1 %in% colnames(x = X)) ||
      !any(lbl2 %in% colnames(x = X)) ||
      !any(lbl3 %in% colnames(x = X)) ||
      !any(lbl4 %in% colnames(x = X))) {
    stop("the RHS of formula did not contain an appropriate exposure() object",
         call. = FALSE)
  }

  i1 <- colnames(x = X) %in% lbl1
  i2 <- colnames(x = X) %in% lbl2
  i3 <- colnames(x = X) %in% lbl3
  i4 <- colnames(x = X) %in% lbl4

  # Vtime, Vtype, infection.time, infection.type
  Vtime <- X[,i1]
  Vtype <- X[,i2]
  infection.time <- X[,i3]
  infection.type <- X[,i4]
  # added drop=FALSE to properly handle case when only 1 covariate is in the model.
  X <- X[,-c(which(i1),which(i2),which(i3),which(i4)), drop = FALSE]

  if (ncol(x = X) == 0L) {
    stop("Model doesn't have any covariates",
         call. = FALSE)
  }

  # extract outcome-related variables that correspond to the returns by outcome()
  # subject_id, entry_time, event_time, censor_time
  dt <- suppressMessages(stats::model.response(data = mf))

  if(!all.equal(colnames(dt), c("subject_id", "entry_time", "event_time", "censor_time"))) {
    stop("the LHS of formula did not contain an appropriate outcome() object",
         call. = FALSE)
  }

  # subject.id, event.time, censor.time, entry.time
  subject.id <- dt[,1L]
  entry.time <- dt[,2L]
  event.time <- dt[,3L]
  censor.time <- dt[,4L]

  # prepare input data for CoxReg
  # subject.id, event.time, censor.time, entry.time, Vtime, Vtype,
  # infection.time, infection.type, X
  input_data <- data.frame("subject.id" = subject.id,
                           "event.time" = event.time,
                           "censor.time" = censor.time,
                           "entry.time" = entry.time,
                           "Vtime" = Vtime,
                           "Vtype" = Vtype,
                           "infection.time" = infection.time,
                           "infection.type" = infection.type)

  # sort input data by subject_id and infection time
  input_data <- cbind(input_data, X)
  input_data <- input_data[order(input_data$subject.id,
                                 input_data$infection.time),]

  # remove any cases that have NA in the input data
  use <- stats::complete.cases(input_data)

  if (all(!use, na.rm = TRUE)) {
    stop("input checks result in all NA -- verify inputs",
         call. = FALSE)
  }

  if (any(!use, na.rm = TRUE)) {
    input_data <- input_data[use]
    message(sum(!use, na.rm = TRUE),
            " cases removed from the analysis due to NA values")
  }

  # process knots
  tau = max(censor.time)
  Kvac = length(unique(Vtype))
  Kinf = length(unique(infection.type))

  if(is.null(prior_infection_knots)) {

    message("prior_infection_knots has zero length; default values will be used")
    prior_infection_knots = lapply(1:Kinf,
                                   function(i) c(reinfection_cutoff, seq(90,
                                                                         tau-min(infection.time[infection.type==i]),
                                                                         90)))
  }

  if(!vaccine_infection_interaction) {
    if(is.null(vaccine_knots)) {
      message("vaccine_knots has zero length; default values will be used")
      vaccine_knots = lapply(1:Kvac,
                             function(i) seq(30,
                                             tau-min(Vtime[Vtype==i]),
                                             30))
    }

    knots <- c(vaccine_knots,
               prior_infection_knots)

  } else{

    if(is.null(vaccine_uninfected_knots)) {
      message("vaccine_uninfected_knots has zero length; default values will be used")
      vaccine_uninfected_knots = lapply(1:Kvac,
                                        function(i) seq(30,
                                                        tau-min(Vtime[Vtype==i]),
                                                        30))
    }

    if(is.null(vaccine_infected_knots)) {
      message("vaccine_uninfected_knots has zero length; default values will be used")
      vaccine_infected_knots = lapply(1:Kvac,
                                      function(i) seq(60,
                                                      tau-min(Vtime[Vtype==i]),
                                                      60))
    }

    knots <- c(vaccine_uninfected_knots,
               vaccine_infected_knots,
               prior_infection_knots)

  }

  # number of parameters
  dbeta = ncol(X)
  if(!vaccine_infection_interaction) {
    dgamma1 = lengths(vaccine_knots)+1-last_piece_constant
    dgamma2 = lengths(prior_infection_knots)+2-last_piece_constant
    d = c(dbeta, dgamma1, dgamma2)
  } else {
    dgamma1 = lengths(vaccine_uninfected_knots)+1-last_piece_constant
    dgamma2 = lengths(vaccine_infected_knots)+1-last_piece_constant
    dgamma3 = lengths(prior_infection_knots)+2-last_piece_constant
    d = c(dbeta, dgamma1, dgamma2, dgamma3)
    temp = lapply(related_vaccine_types,
                  function(x) x+length(dgamma1))
    related_vaccine_types = c(related_vaccine_types, temp)
  }
  dcum = cumsum(d)

  # generate restriction matrix
  ResMat = diag(sum(d))

  if(length(related_vaccine_types)>0) {

    keeprow = rep(TRUE,sum(d))

    for(i in 1:length(related_vaccine_types)) {
      type = sort(unlist(related_vaccine_types[i]))
      index = c(dcum[type[1]]+1, dcum[type[2]]+1)
      ResMat[index[1], index[2]] = 1
      keeprow[index[2]] = FALSE
    }

    ResMat = ResMat[keeprow,]

  }


  ### main analysis

  res <- CoxReg(data = input_data,
                knots = knots,
                ResMat = ResMat,
                MaxIter = 500,
                constantVE = last_piece_constant,
                interact = vaccine_infection_interaction,
                cutoff = reinfection_cutoff,
                plots = plots)

  return( res )
}
