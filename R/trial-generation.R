
#' Simulation thresholds, parameters and structured containers
#' @export


gen_input <- function(pop, ln_cr, subtypes_cr, signatures_cr, return_to_env = FALSE){

  if (missing(pop)){
    stop("The 'pop' argument is missing. Provide the patients' population")
  }

  ln <- ln_cr
  subtypes <- subtypes_cr
  signatures <- signatures_cr

  sim_thld <- list(
    shape   = 1.05,       # shape parameter for the weibull distribution
    nmax    = 1000,        # max n patients per treatment x signature
    nval    = 100,         # minimum number for evaluating treatments x signature
    nupdate = 40,         # number of enrolled patients in the active arms before updating randomization prob
    nmin2   = 80,          # minimum number of enrolled patients in each treatment-subgroup combination for evaluation
    mxmonth = 36,         # n months follow-up
    k       = 60,         # n patients recruited per month
    pU      = 0.90,       # thresholds for early stopping: graduation
    pU2     = 0.85,        # thresholds subtypes graduation probability
    pU2n    = 0.50,       # same as pU2 but for n < nmin2
    pU2all  = 0.50,       # same as pU2 but for signature all
    pL      = 0.10,       # thresholds for early stopping: futility
    pL2     = 0.20,       # similar to PU2, but for futility
    pL2n    = 0.50,       # same as pL2 but for n < nmin2
    delta   = 0,          # clinical survival time difference (aft > delta)
    pow     = 2,           # power for adapting randomization probabilities
    dna_a   = 0.25,        # parameter of beta distribution for cdDNA
    dna_b   = 1.10        # parameter of beta
  )

  # percentage of treatments in the control group
  p_ctr <- switch(pop,
                  cr = c("Other" = 1, "ARSi" = 0, "Taxane_CT" = 0),
                  hs = c("Other" = .3, "ARSi" = .35, "Taxane_CT" = .35)
  )

  r_gx <- switch(pop,
                 cr = matrix(
                   c(
                     0.8 , 0.20, 0  , 0  ,
                     0.55, 0.45, 0  , 0  ,
                     0.4 , 0.2 , 0.1, 0.1,
                     0.3 ,	0.7 , 0 , 0  ,
                     0   , 0   , 0.5, 0.5,
                     0   , 0   , 0.5, 0.5,
                     0   , 0   , 0.5, 0.5,
                     0   , 0   , 0.5, 0.5,
                     0.3 , 0.7 , 0 , 0   ,
                     0.3 , 0.7 , 0 , 0   ,
                     0.2 , 0.6 , 0.1, 0.1,
                     0.2 , 0.6 , 0.1, 0.1,
                     0   , 0   , 0.5, 0.5,
                     0   , 0   , 0.5, 0.5,
                     0   , 0   , 0.5, 0.5,
                     0   , 0   , 0.5, 0.5
                   ), byrow = T, ncol = ln$X-1,
                   dimnames = list(ln$names_type, ln$X_names[-1])),
                 hs = matrix(
                   c(
                     0.6, 0.4, 0  ,
                     0.3, 0.6, 0.1,
                     0.3, 0.7, 0  ,
                     0.2, 0.7, 0.1,
                     0  , 0  , 1  ,
                     0  , 0  , 1  ,
                     0  , 0  , 1  ,
                     0  , 0  , 1
                   ), byrow = T, ncol = ln$X-1,
                   dimnames = list(ln$names_type, ln$X_names[-1]))
  )
  r_gx0 <- cbind(Control = apply(r_gx, 1, max), r_gx)
  r_gx0 <- r_gx0 + rowMeans(r_gx0)
  r_gx0 <- r_gx0/rowSums(r_gx0)
  r_gx0_sign <- t(apply(subtypes[, ln$names_sign], 2,
                        function(s) colSums(s*r_gx0)/sum(s)))
  # randomization probabilities
  r_prob <- list(
    type = r_gx0,
    sign = r_gx0_sign
  )
  lapply(r_prob, round, 2)

  # tuning parameters for the Bayesian (weibull) model
  tun_param <- list(
    prior = normal(0, .5),                   # prior on the beta parameters
    prior_aux = exponential(1),              # prior on the shape parameter
    prior_intercept = normal(switch(pop, "hs" = 4, "cr" = 3), 1),
    # prior on the intercept parameter
    basehaz = "weibull-aft",                 # specification of baseline hazard (note "weibull" is parametrized as proportional hazard model)
    algorithm = "sampling",                  # algorithm for fitting Bayesian model (sampling = full bayesian, meanfield/fullrank = approximate bayesian)
    chains = 1,                              # number of chains
    iter = 200,                             # number of iteration
    warmup = 100,                            # number of warmup
    postp_defna = rep(NA, ln$X - 1 + 2) %>%             # default names for coefficient of the Weibull model
      setNames(nm = c("(Intercept)", paste0("trt", ln$X_names[-1]), "weibull-shape"))
  )

  # results at stop
  at_stop <- list(
    stop = matrix(NA, ncol = ln$X - 1, nrow = ln$S,                # NA if no early stopping otherwise 1 = graduation, 2 = futility, 3 = maximum sample size
                  dimnames = list(ln$names_sign, ln$X_names[-1])),
    stop_type = matrix(NA, ncol = ln$X - 1, nrow = ln$G,           # similar to before but for subtypes (1 for early stops)
                       dimnames = list(ln$names_type, ln$X_names[-1])),
    t_trt = matrix(NA, ncol = ln$X - 1, nrow = ln$S,               # time by early stop
                   dimnames = list(ln$names_sign, ln$X_names[-1])),
    n_sign = matrix(NA, ncol = ln$X - 1, nrow = ln$S,              # n in the signature by early stop
                    dimnames = list(ln$names_sign, ln$X_names[-1])),
    n_sign_ctrl = matrix(NA, ncol = ln$X - 1, nrow = ln$S,         # n in the control group
                         dimnames = list(ln$names_sign, ln$X_names[-1])),
    post_p = matrix(NA, ncol = ln$X - 1, nrow = ln$S,              # prob of superiority at early stop
                    dimnames = list(ln$names_sign, ln$X_names[-1])),
    res_type = list()                                              # further results by subtypes
  )

  # simulated data and summary statistics
  sim_dat <- data.frame(
    id = numeric(),                                          # id var (useful for re-randomization)
    trt = factor(levels = ln$X_names),                       # treatment allocation
    trt_noctrl = factor(levels = c("Other", ln$X_names)),    # treatment allocation without control group
    line = numeric(),                                        # number of treatment line
    subtype = factor(levels = ln$names_type),                # biomarker subtype
    time = numeric(),                                        # vector of real times
    delta = numeric(),                                       # vector of censoring indicators
    month0 = numeric(),                                      # vector of recruitment times
    time_m = numeric(),                                       # vector of obs times
    ctDNA_fu = numeric()
  )
  sum_stat <- list(
    type = list(
      n = matrix(rep(0, ln$X*ln$G), ncol = ln$X, nrow = ln$G,      # patients enrolled
                 dimnames = list(ln$names_type, ln$X_names)),
      delta = matrix(rep(0, ln$X*ln$G), ncol = ln$X, nrow = ln$G,  # censoring indicator
                     dimnames = list(ln$names_type, ln$X_names)),
      PT = matrix(rep(0, ln$X*ln$G), ncol = ln$X, nrow = ln$G,     # observed person-months
                  dimnames = list(ln$names_type, ln$X_names)),
      t = matrix(rep(0, ln$X*ln$G), ncol = ln$X, nrow = ln$G,      # sum(PT^shape)
                 dimnames = list(ln$names_type, ln$X_names))
    ),
    sign = list(                                                   # same as before but for signatures
      n = matrix(rep(0, ln$X*ln$S), ncol = ln$X, nrow = ln$S,
                 dimnames = list(ln$names_sign, ln$X_names)),
      delta = matrix(rep(0, ln$X*ln$S), ncol = ln$X, nrow = ln$S,
                     dimnames = list(ln$names_sign, ln$X_names)),
      PT = matrix(rep(0, ln$X*ln$S), ncol = ln$X, nrow = ln$S,
                  dimnames = list(ln$names_sign, ln$X_names))
    )
  )

  # list of objects
  obj <- list(ln = ln, subtypes = subtypes, signatures = signatures,
              sim_thld = sim_thld, p_ctr = p_ctr, r_prob = r_prob, tun_param = tun_param,
              at_stop = at_stop, sim_dat = sim_dat, sum_stat = sum_stat)

  if (return_to_env == TRUE){
    list2env(obj, globalenv())
  } else {
    return(obj)
  }

}



#' Main sim function
#' @export

sim_trial <- function(ri_gx, ri_gx2, f_alter = 1, smart = FALSE,
                      intermediate = 0, progress = 0, input){

  res <- input$at_stop                                # results at early stop
  res_month <- list()                           # interim results to be filled if intermediate == TRUE
  month <- 0                                    # running calendar time in months

  sim_dat <- input$sim_dat
  sum_stat <- input$sum_stat
  repeat{

    # simulate new patients data ----
    if (any(input$r_prob$type[, -1] != 0)) sim_dat <-
        next_cohort(sim_dat, month, ri_gx, ri_gx2, input)
    # update data and summary statistics ----
    month    <- month + 1
    sim_dat  <- update_dat(sim_dat, month)
    sum_stat <- update_sumstat(sim_dat, input)


    # stop after mxmonth months
    if (month >= input$sim_thld$mxmonth){
      break
    }


    if (sum(sum_stat$sign$n[1, -1]) >= input$sim_thld$nupdate){

      # run t tests

      tres <- do.call(rbind, lapply(input$subtypes$subtypes, function(st) {

        sapply(input$ln$X_names[-1], function(trti) {

          idat <- subset(sim_dat, subtype == st & trt %in% c("Control", trti))
          if(!all(table(factor(idat$trt)) > 2) | length(unique(idat$trt)) != 2) {
            return(NA)
          } else {
            tres <- t.test(log(time) ~ trt, data = idat)
            tres$statistic
          }
        })
      }))
      rownames(tres) <- input$subtypes$subtypes
      colnames(tres) <- input$ln$X_names[-1]



      # evaluate treatments and update randomization probabilities ----
      res    <- eval_trt(res, tres, month, sum_stat, input)
      input$r_prob <- update_rand(input$r_prob, tres, res, f_alter = f_alter, input)
    }

  }

  return(list(res = res,
              r_prob = input$r_prob, sim_dat = sim_dat))

}



#' Simulate next cohort of patients
#' @export
next_cohort <- function(sim_dat, month, ri_gx, ri_gx2, input){

  # n <- structure(c(rmultinom(1, size = sim_thld$k, prob = subtypes$prev)),
  #                names = ln$names_type) # simulate k subtypes
  prev_subtype <- c(MCMCpack::rdirichlet(1, 100*input$subtypes$prev))
  n <- structure(c(rmultinom(1, size = input$sim_thld$k, prob = prev_subtype)),
                 names = input$ln$names_type) # simulate k subtypes
  # stop randomize patients when no treatment is available
  if (any(rowSums(input$r_prob$type[, -1]) == 0)){
    n[input$ln$names_type[rowSums(input$r_prob$type[, -1]) == 0]] <- 0
  }
  ngx <- Map(function(n, r) c(rmultinom(1, size = n, prob = r)),
             n, split(input$r_prob$type, 1:nrow(input$r_prob$type))) %>%
    do.call("rbind", .) %>% `colnames<-`(input$ln$X_names) # treatment randomization (subtypes)
  # which treatment has been assigned to the control
  n_ctr <- lapply(ngx[, "Control"], function(ni) c(rmultinom(1, size = ni, prob = input$p_ctr))) %>%
    do.call("rbind", .) %>% `colnames<-`(names(input$p_ctr))

  Tilist <- Map(function(m, n){ # simulate real PFS time
    #rweibull(n, shape = sim_thld$shape, scale = 1/m^(1/sim_thld$shape))
    rlnorm(n, sdlog = .3, meanlog = m)

  }, cbind(ri_gx[, colnames(n_ctr)], ri_gx[, colnames(ngx[, -1])]), cbind(n_ctr, ngx[, -1]))

  ctDNAlist <- Map(function(m, n){ # simulate real PFS time
    #rbeta(n, shape1 = sim_thld$dna_a, shape2 = sim_thld$dna_b / m)
    rnorm(n, mean = m, sd = .3)

  }, cbind(ri_gx2[, colnames(n_ctr)], ri_gx2[, colnames(ngx[, -1])]),
  cbind(n_ctr, ngx[, -1]))


  # record new simulated data
  if (sum(n) > 0){
    sim_dat[nrow(sim_dat) + seq(sum(n)), ] <- data.frame(
      id = nrow(sim_dat) + seq(sum(n)),
      trt = factor(rep(1:input$ln$X, colSums(ngx)), levels = 1:input$ln$X, labels = input$ln$X_names),
      trt2 = factor(rep(c(seq_along(input$p_ctr), 2:input$ln$X), colSums(cbind(n_ctr, ngx[, -1]))),
                    levels = 1:input$ln$X, labels = c("Other", input$ln$X_names[-1])),
      line = 1,
      subtype = factor(
        unlist(apply(cbind(n_ctr, ngx[, -1]), 2, function(ni) rep(1:input$ln$G, ni))),
        levels = 1:input$ln$G, labels = input$ln$names_type),
      time =  unlist(Tilist),
      delta = 0,
      month0 = month,
      time_m = 1,
      ctDNA_fu = unlist(ctDNAlist)
    )
  }


  return(sim_dat)
}



#' update simulate data based on current month
#' @export
update_dat <- function(sim_dat, month){

  sim_dat[month > sim_dat$month0 + sim_dat$time, "delta"] <- 1
  sim_dat[month > sim_dat$month0 + sim_dat$time, "time_m"] <-
    sim_dat[month > sim_dat$month0 + sim_dat$time , "time"]
  sim_dat[month <= sim_dat$month0 + sim_dat$time, "time_m"] <-
    month - sim_dat[month <= sim_dat$month0 + sim_dat$time , "month0"]

  return(sim_dat)
}


#' update summary statistics based on current month
#' @export
update_sumstat <- function(sim_dat, input){

  # by subtype
  input$sum_stat$type$n_unique <- with(sim_dat[sim_dat$line == 1, ], tapply(id, list(subtype, trt), function(x) length(x), default = 0))
  input$sum_stat$type$n <- with(sim_dat, tapply(id, list(subtype, trt), length, default = 0))
  input$sum_stat$type$delta <- with(sim_dat, tapply(delta, list(subtype, trt), sum, default = 0))
  input$sum_stat$type$PT <- with(sim_dat, tapply(time_m, list(subtype, trt), sum, default = 0))
  input$sum_stat$type$rate <- input$sum_stat$type$delta/input$sum_stat$type$PT
  # by signature
  input$sum_stat$sign$n_unique <- t(apply(input$subtypes[, input$ln$names_sign], 2,
                                          function(x) colSums(x*input$sum_stat$type$n_unique)))
  input$sum_stat$sign$n <- t(apply(input$subtypes[, input$ln$names_sign], 2, function(x) colSums(x*input$sum_stat$type$n)))
  input$sum_stat$sign$delta <- t(apply(input$subtypes[, input$ln$names_sign], 2, function(x) colSums(x*input$sum_stat$type$delta)))
  input$sum_stat$sign$PT <- t(apply(input$subtypes[, input$ln$names_sign], 2, function(x) colSums(x*input$sum_stat$type$PT)))
  input$sum_stat$sign$rate <- input$sum_stat$sign$delta/input$sum_stat$sign$PT

  return(input$sum_stat)
}


#' evaluate treatment for early stop (graduation, futility, max sample size)
#' @export
eval_trt <- function(res, tres, month, sum_stat, input){

  # superiority
  sub_sup <- which(!is.na(tres) & 1 - pnorm(tres) >= input$sim_thld$pU &
                     sum_stat$type$n[, -1] >= input$sim_thld$nval &
                     (is.na(res$stop_type) | res$stop_type == 0), arr.ind = T)

  # futility
  sub_inf <- which(!is.na(tres) & 1 - pnorm(tres) <= input$sim_thld$pL &
                     sum_stat$type$n[, -1] >= input$sim_thld$nval &
                     (is.na(res$stop_type) | res$stop_type == 0), arr.ind = T)

  # maximum sample size
  sub_max <- which(sum_stat$type$n[, -1] >= input$sim_thld$nmax & is.na(res$stop_type), arr.ind = T)
  if(nrow(sub_max) > 0) {
    for(i in 1:nrow(sub_max)) {
      s <- sub_max[i,]
      res$stop_type[s[1], s[2]] <- 0
    }}

  # saving results for graduation/futility
  sub_stop <- cbind(rbind(if (nrow(sub_sup) > 0) cbind(sub_sup, type = 1),
                          if (nrow(sub_inf) > 0) cbind(sub_inf, type = 2)))
  sub_stop <- sub_stop[sub_stop[, "type"] == 1, , drop = FALSE]
  if (!is.null(sub_stop) && nrow(sub_stop) > 0){
    for(i in 1:nrow(sub_stop)) {

      s <- sub_stop[i, ]
      res$stop_type[s[1], s[2]] <- 1
      res$t_trt[s[1], s[2]] <- month
      res$n_sign[s[1], s[2]] <- sum_stat$sign$n[, -1][s[1], s[2]]
      res$n_sign_ctrl[s[1], s[2]] <- sum_stat$sign$n[s[1], 1]

    }
  }

  return(res)

}



#' update randomization probabilities
#' @export
update_rand <- function(r_prob, tres, res, f_alter = 1, input){

  prob_sup <- list(type = 1 - pnorm(tres))
  prob_sup$type[is.na(prob_sup$type)] <- .5
  f_update <- if (f_alter == 1){
    (r_prob$type[, -1]/rowSums(r_prob$type[, -1]))^input$sim_thld$pow + prob_sup$type^input$sim_thld$pow
  } else if (f_alter == 2) {
    prob_sup$type
  } else if (f_alter == 3){
    r_prob$type[, -1]
  } else if (f_alter == 4){
    matrix(1, nrow = ln$G, ncol = ln$X - 1, dimnames = list(ln$names_type, ln$X_names[-1]))
  }
  f_update[!is.na(res$stop_type)] <- 0
  f_update[is.na(f_update)] <- 0 # safety check for errors
  r_prob$type <- cbind(Control = apply(f_update, 1, max), f_update)
  r_prob$type[rowSums(f_update) == 0, 1] <- 1 # in case all the f_update == 0
  # just preventive checks
  r_prob$type[is.na(r_prob$type)] <- 0
  r_prob$type <- r_prob$type/rowSums(r_prob$type) # re-normalize to 1

  # "equivalent" randomization probabilities (proportions) across signatures
  r_prob$sign <- t(sapply(input$subtypes[, input$ln$names_sign],
                          function(x) colSums(x*r_prob$type)/sum(x)))

  return(r_prob)

}



