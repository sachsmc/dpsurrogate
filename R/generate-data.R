#' Generate trial and individual-level data
#' @export
#' @param effect String indicating relationship between seff and yeff, can be "nonlinear", "linear", "null", etc
#' @param Zeffect Strength of association between trial level covariates and outcome
#' @param Ueffect Strength of association between random noise and outcome
#' @param Zgen Function to generate Zeffects (sample from a distribution)
#' @param Sgen Function to generate Seffects (sample from a distribution)

generate_data <- function(effect = "nonlinear", Zeffect = 0, Ueffect = 0,
                          Zgen = function(n) rnorm(n), Sgen = function(n) rnorm(n, sd = 1.5)) {
  ln_cr <-
    list(
      X = 5L,
      G = 16L,
      S = 5L,
      X_names = c("Control", "ARSi",
                  "Taxane_CT", "Platinum_CT", "PARPi"),
      names_type = c(
        "----",
        "---+",
        "--+-",
        "--++",
        "-+--",
        "-+-+",
        "-++-",
        "-+++",
        "+---",
        "+--+",
        "+-+-",
        "+-++",
        "++--",
        "++-+",
        "+++-",
        "++++"
      ),
      names_sign = c("all",
                     "TP53- & AR-", "TP53+", "HRD+", "TEfus+"),
      subtypes_sign = list(
        all = c(
          "----",
          "---+",
          "--+-",
          "--++",
          "-+--",
          "-+-+",
          "-++-",
          "-+++",
          "+---",
          "+--+",
          "+-+-",
          "+-++",
          "++--",
          "++-+",
          "+++-",
          "++++"
        ),
        `TP53- & AR-` = c("----", "---+", "-+--", "-+-+"),
        `TP53+` = c("--+-", "--++", "-++-", "-+++", "+-+-", "+-++",
                    "+++-", "++++"),
        `HRD+` = c("-+--", "-+-+", "-++-", "-+++",
                   "++--", "++-+", "+++-", "++++"),
        `TEfus+` = c("---+", "--++",
                     "-+-+", "-+++", "+--+", "+-++", "++-+", "++++")
      )
    )


  signatures_cr <-
    structure(
      list(
        signatures = c("all", "TP53- & AR-", "TP53+",
                       "HRD+", "TEfus+"),
        `----` = c("X", "X", "", "", ""),
        `---+` = c("X",
                   "X", "", "", "X"),
        `--+-` = c("X", "", "X", "", ""),
        `--++` = c("X",
                   "", "X", "", "X"),
        `-+--` = c("X", "X", "", "X", ""),
        `-+-+` = c("X",
                   "X", "", "X", "X"),
        `-++-` = c("X", "", "X", "X", ""),
        `-+++` = c("X",
                   "", "X", "X", "X"),
        `+---` = c("X", "", "", "", ""),
        `+--+` = c("X",
                   "", "", "", "X"),
        `+-+-` = c("X", "", "X", "", ""),
        `+-++` = c("X",
                   "", "X", "", "X"),
        `++--` = c("X", "", "", "X", ""),
        `++-+` = c("X",
                   "", "", "X", "X"),
        `+++-` = c("X", "", "X", "X", ""),
        `++++` = c("X",
                   "", "X", "X", "X"),
        prev = list(
          all = 1,
          `TP53- & AR-` = 0.4,
          `TP53+` = 0.45,
          `HRD+` = 0.28,
          `TEfus+` = 0.37
        )
      ),
      row.names = c(NA,
                    -5L),
      class = "data.frame"
    )

  subtypes_cr <-
    structure(
      list(
        ARA = c(
          "-",
          "-",
          "-",
          "-",
          "-",
          "-",
          "-",
          "-",
          "+",
          "+",
          "+",
          "+",
          "+",
          "+",
          "+",
          "+"
        ),
        HRD = c(
          "-",
          "-",
          "-",
          "-",
          "+",
          "+",
          "+",
          "+",
          "-",
          "-",
          "-",
          "-",
          "+",
          "+",
          "+",
          "+"
        ),
        TP53 = c(
          "-",
          "-",
          "+",
          "+",
          "-",
          "-",
          "+",
          "+",
          "-",
          "-",
          "+",
          "+",
          "-",
          "-",
          "+",
          "+"
        ),
        TEfus = c(
          "-",
          "+",
          "-",
          "+",
          "-",
          "+",
          "-",
          "+",
          "-",
          "+",
          "-",
          "+",
          "-",
          "+",
          "-",
          "+"
        ),
        prev = c(
          0.25,
          0.05,
          0.15,
          0.1,
          0.05,
          0.05,
          0.03,
          0.03,
          0.05,
          0.04,
          0.04,
          0.04,
          0.03,
          0.03,
          0.03,
          0.03
        ),
        all = c(1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        `TP53- & AR-` = c(1,
                          1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        `TP53+` = c(0,
                    0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
        `HRD+` = c(0,
                   0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
        `TEfus+` = c(0,
                     1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
        subtypes = c(
          "----",
          "---+",
          "--+-",
          "--++",
          "-+--",
          "-+-+",
          "-++-",
          "-+++",
          "+---",
          "+--+",
          "+-+-",
          "+-++",
          "++--",
          "++-+",
          "+++-",
          "++++"
        )
      ),
      out.attrs = list(
        dim = c(
          ARA = 2L,
          HRD = 2L,
          TP53 = 2L,
          TEfus = 2L
        ),
        dimnames = list(
          ARA = c("ARA=-", "ARA=+"),
          HRD = c("HRD=-", "HRD=+"),
          TP53 = c("TP53=-", "TP53=+"),
          TEfus = c("TEfus=-", "TEfus=+")
        )
      ),
      row.names = c(NA, -16L),
      class = "data.frame"
    )

  input <- gen_input("cr", ln_cr, subtypes_cr, signatures_cr, return_to_env = FALSE)


  trtZ <- Zgen(64) #rnorm(64)
  seff <- Sgen(64) #rnorm(64, sd = 1.45)

  UU <- rnorm(64)


  if(effect == "nonlinear") {

    sbeff <- bs(seff, knots = c(-1, 0, 1), degree = 1, Boundary.knots = c(-8, 8))
    yeff <- -1 + sbeff %*% c(-.05, .6, 3, 3.5) +
      + Zeffect * abs(trtZ) + Ueffect * UU
    #plot(yeff ~ seff)

  } else if(effect == "linear") {

    yeff <- -1 + seff * 1 + Zeffect * abs(trtZ) + Ueffect * UU

  } else if(effect == "simple") {

    yeff <- -1 + seff * 1 + Zeffect * trtZ + Ueffect * UU

  } else if(effect == "null") {

    yeff <- 0 + Zeffect * abs(trtZ) + Ueffect * UU

  } else if(effect == "inter") {

    yeff <- ifelse(trtZ < 0, 0, 1) * seff + Ueffect * UU

  } else if(effect == "interhide") {

    yeff <- ifelse(UU < 0, 0, 1) * seff + Zeffect * trtZ

  } else if(effect == "onetrt") {

    trttmp <- rep(input$ln$X_names[-1],  each = 16)
    yeff <- ifelse(trttmp == "ARSi", 1, 0) * seff + Zeffect * trtZ + Ueffect * UU

  } else if(effect == "twotrt") {

    trttmp <- rep(input$ln$X_names[-1],  each = 16)
    trtZ <- rnorm(64, ifelse(trttmp == "ARSi", -2, ifelse(trttmp == "Taxane_CT", 2, 0)))
    yeff <- ifelse(trttmp %in% c("ARSi", "Taxane_CT"), seff - min(seff), 0) +
      Zeffect * trtZ + Ueffect * UU

  } else if(effect == "manybiom") {

    biomtmp <- rep(input$ln$names_type, 4)
    biomnum <- as.numeric(as.factor(biomtmp))

    biocut <- as.numeric(cut(biomnum, c(0, 6, 11, 16), right = TRUE))

    seff <- rnorm(64, mean = c(-1, 2, 0)[biocut], sd = .75)

    trtZ <- rnorm(64, mean = biocut - 2, sd = .25)

    yeff <- ifelse(biocut == 1,  .75 * seff^2 + min(seff) + 4,
                   ifelse(biocut == 2, 2 * seff + min(seff) - 2,
                          0)) +
      Zeffect * trtZ + Ueffect * UU


  }

  # treatment effects
  ri_gx <- matrix(
    c(rep(3, 16), 3 - yeff),
    ncol = input$ln$X,
    nrow = input$ln$G,
    byrow = FALSE,
    dimnames = list(input$ln$names_type, c("Other", input$ln$X_names[-1]))
  )

  ri_gx2 <-
    matrix(
      c(rep(3, 16), 3 - seff),
      ncol = input$ln$X,
      nrow = input$ln$G,
      byrow = FALSE,
      dimnames = list(input$ln$names_type, c("Other", input$ln$X_names[-1]))
    )

  J <- paste0(rownames(ri_gx), ":", rep(colnames(ri_gx)[-1], each = 16))

  test <- sim_trial(
    ri_gx,
    ri_gx2,
    f_alter = 1,
    smart = FALSE,
    intermediate = 0,
    progress = 0,
    input
  )


  test$sim_dat$sig.tp53.arneg <-
    ifelse(sapply(as.character(test$sim_dat$subtype), function(x) {
      signatures_cr[2, x]
    }) == "X", 1, 0)
  test$sim_dat$sig.tp53pos <-
    ifelse(sapply(as.character(test$sim_dat$subtype), function(x) {
      signatures_cr[3, x]
    }) == "X", 1, 0)
  test$sim_dat$sig.hrdpos <-
    ifelse(sapply(as.character(test$sim_dat$subtype), function(x) {
      signatures_cr[4, x]
    }) == "X", 1, 0)
  test$sim_dat$sig.tefus <-
    ifelse(sapply(as.character(test$sim_dat$subtype), function(x) {
      signatures_cr[5, x]
    }) == "X", 1, 0)

  tdat <- test$sim_dat

  txs <- c("ARSi", "Taxane_CT", "Platinum_CT", "PARPi")
  ldat <- NULL
  for(j in txs) {
    tmp <- subset(tdat, trt %in% c(j, "Control"))
    ldat <- rbind(ldat, data.frame(subtype = tmp$subtype, time = tmp$time, ctDNA_fu = tmp$ctDNA_fu,
                                   trt = ifelse(tmp$trt == "Control", 0, 1), trttype = j))

  }
  ldat$J <- factor(paste0(ldat$subtype, ":", ldat$trttype))

  rdat <- data.frame(
    seff = -seff,
    yeff = -yeff,
    trtZ = trtZ,
    J = J
  )
  tenn <- table(ldat$J)
  rdat <- merge(rdat, data.frame(ngrp = as.vector(tenn), J = names(tenn)), by = "J")

  list(ldat = ldat, rdat = rdat)

}
