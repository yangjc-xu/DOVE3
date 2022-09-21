#' @export
CoxReg = function(data, knots, Xg = NULL, ng = NULL,
                  ResMat, eps = 1e-4, MaxIter = 15,
                  constantVE = FALSE, conditional = TRUE,
                  interact = FALSE, cutoff = 14){
  ## data: subject.id, event.time, censor.time, vaccine.time, Vtime, Vtype, infection.time, infection.type, X
  tau = max(data$event.time)
  n = length(unique(data$subject.id))
  m = as.vector(table(data$subject.id))
  MaxReccur = max(m)
  index_dup = which(duplicated(data$subject.id))
  data_unique = data[-index_dup,]
  X = model.matrix(~., data_unique[,9:ncol(data)])[,-1]
  Time = data$event.time
  Delta = ifelse(data$event.time < data$censor.time, 1, 0)
  t = sort(unique(Time[which(Delta==1)]), decreasing = T)
  C = data_unique$censor.time
  VacTime = data_unique$vaccine.time
  FirstInfTime = data_unique$infection.time
  S = NULL
  V = NULL
  n_cohort = max(data$Vtype)
  n_inf = max(data$infection.type)
  interact_ind = 0
  infection_ind = 0
  S_vac = as.list(data_unique$Vtime)
  V_vac = as.list(data_unique$Vtype)

  if(conditional){
    Xg = matrix(0, nrow = 1, ncol = dim(X)[2])
    ng = 0
  }

  if(interact){
    S_vac = computeList(S_vac, S_vac, (1:n)-1)
    V_vac = computeList(V_vac, as.list(data_unique$Vtype + n_cohort), (1:n)-1)
    interact_ind = n_cohort + 1
    infection_ind = 2*n_cohort + 1
  } else{
    infection_ind = n_cohort + 1
  }

  for (i in 1:MaxReccur) {
    if(i == 1){
      S = computeList(S_vac, as.list(data_unique$infection.time), (1:n)-1)
      V = computeList(V_vac, as.list(data_unique$infection.type + infection_ind - 1), (1:n)-1)
    } else{
      subid_i = which(m == i)
      loc = which(data$subject.id%in%subid_i)
      loc = loc[seq(i,length(loc),i)]
      time_i = data$infection.time[loc]
      type_i = data$infection.type[loc]
      S_tmp = rep(tau, n)
      S_tmp[subid_i] = time_i
      V_tmp = rep(infection_ind, n)
      V_tmp[subid_i] = type_i + infection_ind - 1
      S = computeList(S, as.list(S_tmp), (1:n)-1)
      V = computeList(V, as.list(V_tmp), (1:n)-1)
    }
  }

  dimension = numeric(length(knots))
  for (i in 1:length(dimension)) {
    dimension[i] = length(knots[[i]]) + 1
  }
  if(constantVE){
    dimension = dimension - 1
  }

  result = Cox_general(Time, t, Delta, C, V,
                       S, X, m, Xg, ng, VacTime, FirstInfTime,
                       dimension, infection_ind, interact_ind,
                       knots, ResMat, eps, MaxIter,
                       constantVE, conditional, interact, cutoff)
  return(result)
}
