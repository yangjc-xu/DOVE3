#' @export
#' @import Rcpp
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

  if (ncol(x = X) == 0L) X <- NULL

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
  input_data <- cbind(input_data, X)

  # remove any cases that have NA in the
  # entry_time, event_time, event_status, covariates, or vacStatus
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
                                   function(i) c(reinfection_cutoff, seq(90,max((event.time-infection.time)[infection.type==i]), 90)))
  }

  if(!vaccine_infection_interaction) {
    if(is.null(vaccine_knots)) {
      message("vaccine_knots has zero length; default values will be used")
      vaccine_knots = lapply(1:Kvac,
                             function(i) seq(30, max((event.time-Vtime)[Vtype==i]), 30))
    }

    knots <- c(vaccine_knots,
               prior_infection_knots)

  } else{

    if(is.null(vaccine_uninfected_knots)) {
      message("vaccine_uninfected_knots has zero length; default values will be used")
      vaccine_uninfected_knots = lapply(1:Kvac,
                                        function(i) seq(30, max((event.time-Vtime)[Vtype==i]), 30))
    }

    if(is.null(vaccine_infected_knots)) {
      message("vaccine_uninfected_knots has zero length; default values will be used")
      vaccine_infected_knots = lapply(1:Kvac,
                                      function(i) seq(30, max((event.time-Vtime)[Vtype==i]), 60))
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
  keeprow = rep(TRUE,sum(d))
  for(i in 1:length(related_vaccine_types)) {
    type = sort(unlist(related_vaccine_types[i]))
    index = c(dcum[type[1]]+1, dcum[type[2]]+1)
    ResMat[index[1], index[2]] = 1
    keeprow[index[2]] = FALSE
  }
  ResMat = ResMat[keeprow,]


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
