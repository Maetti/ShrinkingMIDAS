midas_test <- function() {


   #   ____________________________________________________________________________
   #   Set Inputs                                                              ####

   library(dplyr)
   library(ggplot2)

   ## set logging file
   logger::log_appender(logger::appender_tee(here::here("inst", "simulation", "logging", "model_check",
                                                        glue::glue("{gsub('-', '_', Sys.Date())}_model_check_log"))))
   logger::log_info("Setting input variables")

   ## simulations
   nSimulation <- 15
   vSeed <- 1001:(1001 + nSimulation)

   bForceNewData <- FALSE


   ## specifiy general amount
   nY <- 250
   nTrain <- 200

   ## predictor variable setup
   vBeta <- c(1.2, -0.5, 3, 1.1, 0.5, -2.5, -1.5, 0.8, 3.5, rep(0, 41))
   nK = length(vBeta) ## amount of predictors
   nFreq = 3
   nLag = 6
   nP = 3

   ## data generating process variables
   nMu = 0.1
   nRho = 0.5
   nVar = 1
   nWithin = 0.5
   nBetween = 0

   ## data generation for response
   .sPoly = "beta"
   nT1 = 1
   nT2 = 3


   bolTrue <- vBeta != 0
   vTrueValue <- vBeta[bolTrue]

   nG <- length(vBeta)
   nGroupSize <- rep(nLag, nG) # group Index

   # dgp_exp_almon_lag(nLag = vLag, 0.005, -0.05)

   ## set seed
   nSeed <- 1001:(1001 + nSimulation)

   ## direction
   sDir <- "inst/simulation"
   # sDir <- here::here("inst/data/simulation")

   nSim <- here::here(sDir, "input")
   nSim <- nSim[1:10]


   #   ____________________________________________________________________________
   #   Create Models (if necessary)                                            ####

   logger::log_info("Creating Stan models")

   lCMDmodels <- vector("list", length(stanmodels))
   dirStan <- here::here("inst/stan")
   sStanFiles <- dir(dirStan)[grep("stan", dir(dirStan))]

   for (i in seq_along(sStanFiles)) {
      sStanModel <- glue::glue("{dirStan}/{sStanFiles[i]}")
      lCMDmodels[[i]] <- cmdstanr::cmdstan_model(stan_file = sStanModel)
      #stanc_options = list("01"))
   }

   sModelName <- stringr::str_to_title(gsub("_", " ", stringr::str_extract(sStanFiles, "[^\\.]+")))

   names(lCMDmodels) <- sModelName

   logger::log_info("Stan models creation done: {length(sModelName)} models")

   lCMDmodels <- lCMDmodels[-which(names(lCMDmodels) == "Horseshoe Plus Wo Group")]
   sModelName <- sModelName[!(sModelName == "Horseshoe Plus Wo Group")]

   #   ____________________________________________________________________________
   #   Prepare Posterior Beta                                                  ####


   #   ____________________________________________________________________________
   #   Create Models (if necessary)                                            ####

   logger::log_info("Creating Stan models")

   lCMDmodels <- vector("list", length(stanmodels))
   dirStan <- here::here("inst/stan")
   sStanFiles <- dir(dirStan)[grep("stan", dir(dirStan))]

   for (i in seq_along(sStanFiles)) {
      sStanModel <- glue::glue("{dirStan}/{sStanFiles[i]}")
      lCMDmodels[[i]] <- cmdstanr::cmdstan_model(stan_file = sStanModel)
      #stanc_options = list("01"))
   }

   sModelName <- stringr::str_to_title(gsub("_", " ", stringr::str_extract(sStanFiles, "[^\\.]+")))

   names(lCMDmodels) <- sModelName

   logger::log_info("Stan models creation done: {length(sModelName)} models")

   lCMDmodels <- lCMDmodels[-which(names(lCMDmodels) == "Horseshoe Plus Wo Group")]
   sModelName <- sModelName[!(sModelName == "Horseshoe Plus Wo Group")]

   #   ____________________________________________________________________________
   #   Prepare Posterior Beta                                                  ####

   arrRawOut <- arrow::open_dataset(here::here("inst", "simulation",
                                               "output", "02_raw_extracted"),
                                    partitioning = c("model", "simulation"))

   ## models to compare
   vMDcompare <- list.files(here::here("inst", "simulation", "output", "02_raw_extracted"))
   lSeriesModel <- vector("list", length = length(vMDcompare))


   #   ____________________________________________________________________________
   #   MIDAS TESTING                                                           ####


   i <- 1
   ## get simulated data
   nSimRunSave <- helper_create_number_name(i)

   sDirInputLoad <- here::here(sDir, "input", glue::glue("{nSimRunSave}_rawdata.rds"))

   lData <- readRDS(sDirInputLoad)

   y_train <- lData$model_data$y_train
   y_test <- lData$model_data$y_test

   x_train <- lData$x_raw[1:(6 + 3*nTrain), ]
   x_test <- lData$x_raw[(6 + 3*nTrain + 1):nrow(lData$x_raw), ]


   midasr::mls(x_train[, 1], 0:5, 3, midasr::nealmon)

   midasr::fmls(x_train[, 1], k = 5, m = 3, midasr::nealmon)

   head(midasr::mls(x_train[, 1], 0:5, 3, midasr::nealmon), n = 10)

   mX1 <- midasr::mls(x_train[, 1], 0:5, 3, midasr::nealmon)[2:201, ]
   mX2 <- lData$x_align[[1]][1:200, ]

   all(mX1 == mX2)



   ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
   ### Estimate MIDAS and Prediction                                           ####

   # x <- x_train[, 1]
   # eq_r <- midasr::midas_r(formula = formula(y_train2 ~ mls(x, 0:5, 3, nealmon) - 1),
   #                         start = list("x" = c(1, -0.5)))

   library(midasr)

   ## prepare data
   y_train2 <- as.matrix(c(NA, y_train, NA))
   sFormula <- as.formula(paste0("y_train2 ~ ", paste0("mls(x_", 1:50, ", 0:5, 3, nealmon)", collapse = " + "), " - 1"))
   lStart <- lapply(1:50, function(x) {c(1, -0.5)})
   names(lStart) <- paste0("x_", 1:50)

   for (i in 1:50) {
      assign(x = paste0("x_", i), x_train[, i])
      assign(x = paste0("x_test_", i), x_test[, i])
   }

   # md_midas <- midasr::midas_r(formula = sFormula, start = lStart,
   #                             Ofunction = "optim", method = "Nelder-Mead",
   #                             weight_gradients = list(nealmon = nealmon_gradient))


   ## estimate model
   md_midas <- midasr::midas_r(formula = sFormula, start = lStart)

   ## check y response
   checkY <- md_midas$model[, 1]
   checkYTrue <- lData$model_data$y_train

   all(checkY == checkYTrue)

   ## check design matrix
   checkDesign <- unname(md_midas$model[, -1])
   checkXTrue <- lData$model_data$x_train
   dim(checkDesign)
   dim(checkXTrue)
   checkDesign[1:5, 1:6]
   checkXTrue[1:5, 1:6]

   all(checkDesign == checkXTrue, na.rm = TRUE)

   ## check if coefficient are the same
   # - yes they are -> can use midas_coefficients on the whole test matrix
   nAlmon1 <- md_midas$coefficients[1:2]

   lCoefTest <- vector("list", length = length(md_midas$coefficients) / 2)
   for (i in 1:(length(md_midas$coefficients) / 2)) {
      lCoefTest[[i]] <- nealmon(p = md_midas$coefficients[(2*i - 1):(2 * i)], d = 6)
   }
   nCoefTest <- unlist(lCoefTest)

   all(nCoefTest == md_midas$midas_coefficients)



   ## prediction data
   xTest <- lData$x_align[[1]][201:250, ]

   all(lData$model_data$x_test[, 1:6] == xTest)

   ## prediction method #1
   mPred1 <- lData$model_data$x_test %*% md_midas$midas_coefficients

   ## prediction method #2
   lTestData <- lapply(1:50, function(x) {lData$x_raw[601:753, x]})
   names(lTestData) <- paste0("x_", 1:50)
   mPred2 <- forecast(md_midas, newdata = lTestData)

   all(round(mPred1, 5) == round(mPred2$mean[-1], 5))
   round(cbind(mPred1, mPred2$mean[-1], lData$model_data$y_test), 4)

   # x_test <- lData$x_raw[(6 + 3*nTrain - 6):nrow(lData$x_raw), ]
   # dim(x_test)
   #
   # midasr::mls(x_test[, 1], 0:5, 3, nealmon)
   #
   # length(lData$x_raw[, 1])
   # length(lData$x_raw[, 1]) / 6
   # xTest2 <- midasr::mls(lData$x_raw[1:756, 1], 0:5, 3, nealmon)
   # xTest2[201:252, ]
   #
   # midasr::mls(lData$x_raw[609:756, 1], 0:5, 3, nealmon)
   #
   #
   #
   #
   # midasr::mls(lTestData[[10]], 0:5, 3, nealmon)
   #
   # midasr::mls(lData$x_raw[598:753, 1], 0:5, 3, nealmon)
   #
   # midasr::mls(lData$x_raw[604:753, 1], 0:5, 3, nealmon)
   #
   #
   # ## calculate own prediction to check with the forecasting method from midasr
   # lAlignTest <- lapply(lData$x_align, function(x) x[201:250, ])
   # vCoef <- coef(midas_model)
   # midas_pred <- function(midas_model, lAlignTest) {
   #
   #    lPredReg <- vector("list", length = length(lAlignTest))
   #    for (i in seq_along(lAlignTest)) {
   #       maTest <- lAlignTest[[i]]
   #       nCoef <- vCoef[(2*i - 1):(2*i)]
   #       nPoly <- nealmon(p = nCoef, d = 6)
   #       lPredReg[[i]] <- maTest %*% nPoly
   #    }
   #
   #    rowSums(do.call("cbind", lPredReg))
   # }

   #   ____________________________________________________________________________
   #   Example from Midas Package                                              ####


   set.seed(1001)
   n <- 250
   trend <- 1:n
   x <- rnorm(4 * n)
   z <- rnorm(12 * n)
   fn_x <- midasr::nealmon(p = c(1, -0.5), d = 8)
   fn_z <- midasr::nealmon(p = c(2, 0.5, -0.1), d = 17)
   y <- 2 + 0.1 * trend + midasr::mls(x, 0:7, 4) %*% fn_x + midasr::mls(z, 0:16, 12) %*% fn_z + rnorm(n)

   eq_u <- midasr::midas_r(y ~ trend + midasr::mls(x, 0:7, 4) + midasr::mls(z, 0:16, 12))

   library(midasr)
   eq_r <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))

   eq_r <- zz_midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon),
                      start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)), Ofunction = "nls")




   eval_f <- function(x, y, z) {
      k <- y * 100
      u <- z + k
      print(u)
      print(k)
      return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
   }


   eval_grad_f <- function(x) {
      return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
                 200 * (x[2] - x[1] * x[1]) ) )
   }

   eval_grad_f(x0)

   eval_grad_f2 <- function(x) {
      numDeriv::grad(func = eval_f, x = x)
   }


   eval_grad_f2(x0)

   x0 <- c( -1.2, 1 )

   opts <- list("algorithm" = "NLOPT_LD_LBFGS",
                "xtol_rel" = 1.0e-8)


   res <- nloptr::nloptr( x0=x0,
                  eval_f=eval_f,
                  eval_grad_f=eval_grad_f,
                  opts=opts, y = 100, z = 20)


   res2 <- nloptr::nloptr( x0=x0,
                          eval_f=eval_f,
                          eval_grad_f=eval_grad_f2,
                          opts=opts)




   #   ____________________________________________________________________________
   #   custom MIDAS                                                            ####

   lInput <-
      lapply(1:50, function(x) {
         list(
            lag = 6,
            theta = c(1, -0.5)
         )
      })


   vTheta <- rep(c(1, -0.5), 50)

   midas_coef <- function(vTheta) {

      sInd <- which(1:length(vTheta) %% 2 == 1)
      nLag <- 6

      for (i in seq_along(sInd)) {

         nT1 <- vTheta[[sInd[[i]]]]
         nT2 <- vTheta[[sInd[[i + 1]]]]

         lCoef[[i]] <- dgp_exp_almon_lag(nLag, nT1, nT2)
      }

      unlist(lCoef)

   }


   midas_optim <- function(vTheta) {

      coeffs <- midas_coef(vTheta)

      y_hat <- X %*% coeffs
      r <- y - y_hat
      sum(r^2)

   }


}



zz_midas_r <- function (formula, data, start, Ofunction = "optim", weight_gradients = NULL, ...) {

      Zenv <- new.env(parent = environment(formula))

      if (missing(data)) {
            ee <- NULL
      }
      else {
            ee <- data_to_env(data)
            parent.env(ee) <- parent.frame()
      }
      if (missing(start)) {
            stop("Please supply starting values.")
      }
      assign("ee", ee, Zenv)
      formula <- as.formula(formula)
      cl <- match.call()
      mf <- match.call(expand.dots = FALSE)
      mf$formula <- formula
      m <- match(c("formula", "data"), names(mf), 0L)
      mf <- mf[c(1L, m)]
      mf[[1L]] <- as.name("model.frame")
      mf[[3L]] <- as.name("ee")
      mf[[4L]] <- as.name("na.omit")
      names(mf)[c(2, 3, 4)] <- c("formula", "data", "na.action")
      itr <- checkARstar(terms(eval(mf[[2]], Zenv)))
      if (!is.null(itr$lagsTable))
            mf[[2]] <- itr$x
      mf <- eval(mf, Zenv)
      mt <- attr(mf, "terms")
      args <- list(...)
      y <- as.numeric(model.response(mf, "numeric"))
      X <- model.matrix(mt, mf)
      if (is.null(ee)) {
            yy <- eval(formula[[2]], Zenv)
      }
      else {
            yy <- eval(formula[[2]], ee)
      }
      y_index <- 1:length(yy)
      if (!is.null(attr(mf, "na.action"))) {
            y_index <- y_index[-attr(mf, "na.action")]
      }
      if (length(y_index) > 1) {
            if (sum(abs(diff(y_index) - 1)) > 0)
                  warning("There are NAs in the middle of the time series")
      }
      ysave <- yy[y_index]
      if (inherits(yy, "ts")) {
            class(ysave) <- class(yy)
            attr(ysave, "tsp") <- c(time(yy)[range(y_index)], frequency(yy))
      }
      if (inherits(yy, c("zoo", "ts"))) {
            y_start <- index2char(index(ysave)[1], frequency(ysave))
            y_end <- index2char(index(ysave)[length(ysave)], frequency(ysave))
      }
      else {
            y_start <- y_index[1]
            y_end <- y_index[length(y_index)]
      }
      prepmd <- prepmidas_r(y, X, mt, Zenv, cl, args, start, Ofunction,
                            weight_gradients, itr$lagsTable)
      prepmd <- c(prepmd, list(lhs = ysave, lhs_start = y_start,
                               lhs_end = y_end))

      class(prepmd) <- "midas_r"
      midas_r.fit(prepmd)
}


prepmidas_r <- function (y, X, mt, Zenv, cl, args, start, Ofunction, weight_gradients,
                         lagsTable, unrestricted = NULL, guess_start = TRUE, tau = NULL) {


      start <- start[!sapply(start, is.null)]
      if (is.null(weight_gradients))
            use_gradient <- FALSE
      else use_gradient = TRUE
      if (!is.null(args$guess_start)) {
            guess_start <- args$guess_start
            args$guess_start <- NULL
      }
      terms.lhs <- as.list(attr(mt, "variables"))[-2:-1]
      dterm <- function(fr, ltb = NULL) {
            term_name <- as.character(fr)[1]
            weight_name <- ""
            rf <- function(p) p
            grf <- function(p) diag(1)
            start <- 0
            freq <- 1
            lagstruct <- 0
            if (term_name %in% c("mls", "fmls", "dmls", "mlsd")) {
                  type <- term_name
                  term_name <- as.character(fr[[2]])
                  wpos <- 5
                  if (type == "mlsd") {
                        freq <- NA
                  }
                  else {
                        freq <- eval(fr[[4]], Zenv)
                  }
                  lags <- eval(fr[[3]], Zenv)
                  nol <- switch(type, fmls = lags + 1, dmls = lags +
                                      1, mls = length(lags), mlsd = length(lags))
                  lagstruct <- switch(type, fmls = 0:lags, dmls = 0:lags,
                                      mls = lags, mlsd = lags)
                  start <- rep(0, nol)
                  grf <- function(p) diag(nol)
                  if (length(fr) > wpos - 1 && fr[[wpos]] != "*") {
                        mf <- fr[-wpos]
                        mf[[1]] <- fr[[wpos]]
                        weight_name <- as.character(fr[[wpos]])
                        noarg <- length(formals(eval(fr[[wpos]], Zenv)))
                        if (noarg < 2)
                              stop("The weight function must have at least two arguments")
                        mf <- mf[1:min(length(mf), noarg + 1)]
                        if (length(mf) > 3) {
                              start_eval <- 4
                              if (type == "mlsd")
                                    start_eval <- 5
                              if (length(mf) >= start_eval) {
                                    for (j in start_eval:length(mf)) {
                                          mf[[j]] <- eval(mf[[j]], Zenv)
                                    }
                              }
                        }
                        mf[[3]] <- ifelse(is.null(ltb), nol, sum(ltb[,
                                                                     1]))
                        rf <- function(p) {
                              mf[[2]] <- p
                              eval(mf, Zenv)
                        }
                        if (use_gradient) {
                              gmf <- mf
                              if (weight_name %in% names(weight_gradients)) {
                                    weight_gradient_name <- paste0(as.character(fr[[2]]),
                                                                   "_tmp_gradient_fun")
                                    gmf[[1]] <- as.name(weight_gradient_name)
                                    assign(weight_gradient_name, weight_gradients[[weight_name]],
                                           Zenv)
                              }
                              else {
                                    gmf[[1]] <- as.name(paste0(weight_name, "_gradient"))
                              }
                              grf <- function(p) {
                                    gmf[[2]] <- p
                                    eval(gmf, Zenv)
                              }
                        }
                        else grf <- NULL
                  }
            }
            list(weight = rf, term_name = term_name, gradient = grf,
                 start = start, weight_name = weight_name, frequency = freq,
                 lag_structure = lagstruct)
      }
      if (is.null(lagsTable)) {
            ltb <- rep(list(NULL), length(terms.lhs))
      }
      else {
            ltb <- lagsTable
            if (attr(mt, "intercept") == 1) {
                  ltb <- ltb[-1]
            }
      }
      rfd <- mapply(dterm, terms.lhs, ltb, SIMPLIFY = FALSE)
      if (attr(mt, "intercept") == 1) {
            intc <- dterm(expression(1))
            intc$term_name <- "(Intercept)"
            rfd <- c(list(intc), rfd)
      }
      rf <- lapply(rfd, "[[", "weight")
      names(rf) <- sapply(rfd, "[[", "term_name")
      weight_names <- sapply(rfd, "[[", "weight_name")
      weight_inds <- which(weight_names != "")
      weight_names <- names(rf)[weight_names != ""]
      start_default <- lapply(rfd, "[[", "start")
      names(start_default) <- names(rf)
      if (length(weight_names) == 0)
            Ofunction <- "lm"
      else {
            if (is.null(start)) {
                  cl$formula <- update_weights(cl$formula, setNames(lapply(1:length(weight_names),
                                                                           function(x) ""), weight_names))
                  warning("Since the start = NULL, it is assumed that U-MIDAS model is fitted")
                  return(eval(cl, Zenv))
            }
            else {
                  if (any(!weight_names %in% names(start)))
                        stop("Starting values for weight parameters must be supplied")
            }
      }
      start_default[names(start)] <- start
      np <- cumsum(sapply(start_default, length))
      pinds <- midasr:::build_indices(np, names(start_default))
      for (i in 1:length(start_default)) names(start_default[[i]]) <- NULL
      if (!is.null(lagsTable)) {
            inones <- function(ones, intro) {
                  ones[ones == 1] <- intro
                  ones
            }
            yname <- all.vars(mt[[2]])
            nms <- names(pinds)
            all_coef2 <- function(p) {
                  pp <- lapply(pinds, function(x) p[x])
                  cr <- c(1, -p[pinds[[yname]]])
                  res <- mapply(function(fun, cf, tb) {
                        restr <- fun(cf)
                        if (is.null(tb)) {
                              restr
                        }
                        else {
                              mltp <- rowSums(apply(tb, 2, inones, restr) %*%
                                                    diag(cr))
                              mltp[rowSums(tb) != 0]
                        }
                  }, rf, pp, lagsTable, SIMPLIFY = FALSE)
                  return(res)
            }
      }
      else {
            all_coef2 <- function(p) {
                  pp <- lapply(pinds, function(x) p[x])
                  res <- mapply(function(fun, param) fun(param), rf,
                                pp, SIMPLIFY = FALSE)
                  return(res)
            }
      }
      initial_midas_coef <- all_coef2(unlist(start_default))
      if (sum(is.na(unlist(initial_midas_coef))) > 0)
            stop("Check your starting values, NA in midas coefficients")
      npx <- cumsum(sapply(initial_midas_coef, length))
      xinds <- midasr:::build_indices(npx, names(start_default))
      if (length(weight_names) > 0 && guess_start) {
            wi <- rep(FALSE, length(rf))
            wi[weight_inds] <- TRUE
            Xstart <- mapply(function(cf, inds, iswhgt) {
                  if (iswhgt) {
                        X[, inds, drop = FALSE] %*% cf
                  }
                  else X[, inds, drop = FALSE]
            }, initial_midas_coef, xinds, wi, SIMPLIFY = FALSE)
            npxx <- cumsum(sapply(Xstart, function(x) {
                  ifelse(is.null(dim(x)), 1, ncol(x))
            }))
            xxinds <- midasr:::build_indices(npxx, names(start_default))
            XX <- do.call("cbind", Xstart)
            prec <- suppressWarnings(lsfit(XX, y, intercept = FALSE))
            lmstart <- lapply(xxinds, function(x) coef(prec)[x])
            names(lmstart) <- names(xxinds)
            for (i in 1:length(lmstart)) names(lmstart[[i]]) <- NULL
            nms <- !(names(start_default) %in% names(start))
            start_default[nms] <- lmstart[nms]
            for (ww in which(wi)) {
                  normalized <- FALSE
                  if (rfd[[ww]]$weight_name %in% c("nealmon", "nbeta",
                                                   "nbetaMT", "gompertzp", "nakagamip", "lcauchyp"))
                        normalized <- TRUE
                  else {
                        normalized <- is_weight_normalized(rf[[ww]],
                                                           start_default[[ww]])
                  }
                  if (normalized) {
                        start_default[[ww]][1] <- lmstart[[ww]]
                  }
            }
      }
      starto <- unlist(start_default)
      all_coef <- function(p) unlist(all_coef2(p))
      mdsrhs <- function(p) {
            coefs <- all_coef(p)
            X %*% coefs
      }
      fn0 <- function(p, ...) {
            r <- y - mdsrhs(p)
            sum(r^2)
      }
      if (!is.null(tau)) {
            fn0 <- function(p, ...) {
                  r <- y - mdsrhs(p)
                  sum(tau * pmax(r, 0) + (tau - 1) * pmin(r, 0))
            }
      }
      if (!use_gradient) {
            gradD <- function(p) jacobian(all_coef, p)
            gr <- function(p) grad(fn0, p)
      }
      else {
            grf <- sapply(rfd, "[[", "gradient")
            pp0 <- lapply(pinds, function(xx) starto[xx])
            grmat0 <- mapply(function(fun, param) fun(param), grf,
                             pp0, SIMPLIFY = FALSE)
            colnos <- sapply(grmat0, ncol)
            rownos <- sapply(grmat0, nrow)
            np <- length(colnos)
            ccol <- cumsum(colnos)
            rrow <- cumsum(rownos)
            pindm <- cbind(c(1, rrow[-np] + 1), rrow, c(1, ccol[-np] +
                                                              1), ccol)
            pindm <- apply(pindm, 1, function(x) list(row = x[1]:x[2],
                                                      col = x[3]:x[4]))
            if (is.null(lagsTable)) {
                  gradD <- function(p) {
                        pp <- lapply(pinds, function(x) p[x])
                        grmat <- mapply(function(fun, param) fun(param),
                                        grf, pp, SIMPLIFY = FALSE)
                        if (length(grmat) == 1) {
                              res <- grmat[[1]]
                        }
                        else {
                              res <- matrix(0, nrow = sum(rownos), ncol = sum(colnos))
                              for (j in 1:length(grmat)) {
                                    ind <- pindm[[j]]
                                    res[ind$row, ind$col] <- grmat[[j]]
                              }
                        }
                        res
                  }
            }
            else {
                  expandD <- function(grm, ltb, cr) {
                        if (is.null(ltb)) {
                              return(grm)
                        }
                        else {
                              el <- lapply(data.frame(ltb), inones2, grm)
                              mltp <- Reduce("+", mapply(`*`, el, cr, SIMPLIFY = FALSE))
                              return(mltp[rowSums(ltb) != 0, ])
                        }
                  }
                  inones2 <- function(ones, intro) {
                        m <- matrix(0, nrow = length(ones), ncol = ncol(intro))
                        if (sum(ones) != nrow(intro))
                              stop("Wrong gradient for AR* term")
                        m[ones == 1, ] <- intro
                        m
                  }
                  expandD2 <- function(fun, param, ltb, nparam = 1) {
                        cf <- fun(param)
                        if (is.null(ltb))
                              return(matrix(0, nrow = length(cf), ncol = nparam))
                        else {
                              mltp <- -apply(ltb, 2, inones, cf)
                              return(mltp[rowSums(ltb) != 0, -1, drop = FALSE])
                        }
                  }
                  dind <- which(names(pinds) == yname)
                  cr <- c(1, -starto[pinds[[dind]]])
                  pp <- lapply(pinds, function(x) starto[x])
                  grmat1 <- mapply(function(fun, param) fun(param),
                                   grf, pp, SIMPLIFY = FALSE)
                  egrmat1 <- mapply(expandD, grmat1, lagsTable, SIMPLIFY = FALSE,
                                    MoreArgs = list(cr))
                  colnos <- sapply(egrmat1, ncol)
                  rownos <- sapply(egrmat1, nrow)
                  np <- length(colnos)
                  ccol <- cumsum(colnos)
                  rrow <- cumsum(rownos)
                  pindm <- cbind(c(1, rrow[-np] + 1), rrow, c(1, ccol[-np] +
                                                                    1), ccol)
                  pindm <- apply(pindm, 1, function(x) list(row = x[1]:x[2],
                                                            col = x[3]:x[4]))
                  gradD <- function(p) {
                        cr <- c(1, -p[pinds[[dind]]])
                        pp <- lapply(pinds, function(x) p[x])
                        grmat <- mapply(function(fun, param) fun(param),
                                        grf, pp, SIMPLIFY = FALSE)
                        egrmat <- mapply(expandD, grmat, lagsTable, SIMPLIFY = FALSE,
                                         MoreArgs = list(cr))
                        res <- matrix(0, nrow = sum(rownos), ncol = sum(colnos))
                        gr_star <- do.call("rbind", mapply(expandD2,
                                                           rf, pp, lagsTable, SIMPLIFY = FALSE, MoreArgs = list(length(pinds[[dind]]))))
                        res[, pinds[[dind]]] <- gr_star
                        for (j in 1:length(egrmat)) {
                              ind <- pindm[[j]]
                              res[ind$row, ind$col] <- egrmat[[j]]
                        }
                        res
                  }
            }
            gr <- function(p) {
                  XD <- X %*% gradD(p)
                  resid <- y - X %*% all_coef(p)
                  as.vector(-2 * apply(as.vector(resid) * XD, 2, sum))
            }
      }
      hess <- function(x) numDeriv::hessian(fn0, x)
      if (is.null(unrestricted)) {
            if (ncol(X) < nrow(X)) {
                  if (attr(mt, "intercept") == 1) {
                        unrestricted <- lm(y ~ ., data = data.frame(cbind(y,
                                                                          X[, -1]), check.names = FALSE))
                  }
                  else {
                        unrestricted <- lm(y ~ . - 1, data = data.frame(cbind(y,
                                                                              X), check.names = FALSE))
                  }
            }
      }
      control <- c(list(Ofunction = Ofunction), args)
      if (!("method" %in% names(control)) & Ofunction == "optim") {
            control$method <- "BFGS"
      }
      term_info <- rfd
      names(term_info) <- sapply(term_info, "[[", "term_name")
      term_info <- mapply(function(term, pind, xind) {
            term$start <- NULL
            term$coef_index <- pind
            term$midas_coef_index <- xind
            term
      }, term_info, pinds[names(term_info)], xinds[names(term_info)],
      SIMPLIFY = FALSE)
      if (!is.null(tau)) {
            gr <- NULL
            hess <- NULL
      }
      list(coefficients = starto, midas_coefficients = all_coef(starto),
           model = cbind(y, X), unrestricted = unrestricted, term_info = term_info,
           fn0 = fn0, rhs = mdsrhs, gen_midas_coef = all_coef, opt = NULL,
           argmap_opt = control, start_opt = starto, start_list = start,
           call = cl, terms = mt, gradient = gr, hessian = hess,
           gradD = gradD, Zenv = Zenv, use_gradient = use_gradient,
           nobs = nrow(X), tau = tau)
}

checkARstar <- function (trms) {
      vars <- as.list(attr(trms, "variables"))[-2:-1]
      env <- environment(trms)
      idx <- which(sapply(vars, function(y) if (length(y) >= 2)
            y[[2]]) == trms[[2]])
      lagsTable <- NULL
      if (length(idx) > 0 && length(vars[[idx]]) >= 5 && vars[[idx]][[5]] ==
          "*") {
            fs <- lapply(sapply(vars, function(y) if (length(y) >=
                                                      4)
                  y[[4]]), eval, env)
            if (length(unique(unlist(fs))) > 1) {
                  lags <- eval(vars[[idx]][[3]], env)
                  push <- lapply(fs, "*", lags)
                  lagsTable <- lapply(1:length(vars), function(w) {
                        z <- vars[[w]]
                        if (length(z) >= 4 && eval(z[[4]], env) != 1) {
                              l <- eval(z[[3]], env)
                              if (length(l) == 1 & as.character(z[1]) %in%
                                  c("fmls", "dmls"))
                                    l <- 0:l
                              tp <- matrix(0, ncol = length(lags) + 1, nrow = max(l) +
                                                 max(push[[w]]) + 1)
                              tp[l + 1, 1] <- 1
                              for (r in 2:ncol(tp)) tp[l + push[[w]][r -
                                                                           1] + 1, r] <- 1
                              tp
                        }
                  })
                  shortSeq <- function(s) {
                        wt <- which(!diff(s) == 1)
                        idx <- c(1, 1 + c(wt, wt - 1), length(s))
                        ams <- s[intersect(1:length(s), idx)]
                        fc <- cumsum(c(TRUE, !round(diff(ams)/2 + head(ams,
                                                                       -1)) %in% s))
                        out <- lapply(split(ams, fc), function(x) if (length(x) ==
                                                                      2)
                              do.call("call", c(":", as.list(x)))
                              else x)
                        names(out) <- NULL
                        out
                  }
                  vars <- lapply(1:length(vars), function(w) {
                        z <- vars[[w]]
                        if (length(z) >= 4 && eval(z[[4]], env) != 1) {
                              fun <- as.character(z[1])
                              l <- eval(z[[3]], env)
                              if (fun %in% c("fmls", "dmls")) {
                                    if (length(l) == 1)
                                          l <- 0:l
                                    else stop("fmls and dmls are not used with a vector of lag orders")
                              }
                              nl <- sort(unique(l + rep(c(0, push[[w]]),
                                                        each = length(l))))
                              if (fun == "mls")
                                    z[[3]] <- do.call("call", c("c", shortSeq(nl)))
                              else if (all(diff(nl) == 1))
                                    z[3] <- max(nl)
                              else if (fun == "fmls") {
                                    z[1] <- call("mls")
                                    z[[3]] <- do.call("call", c("c", shortSeq(nl)))
                              }
                              else stop("Use fmls or mls instead of dmls")
                        }
                        z
                  })
                  icp <- attr(trms, "intercept") == 1
                  trms <- formula(paste(trms[[2]], "~", paste(vars,
                                                              collapse = " + ")), env)
                  if (!icp)
                        trms <- update.formula(trms, . ~ . - 1)
                  else lagsTable <- c(list(NULL), lagsTable)
            }
      }
      list(x = trms, lagsTable = lagsTable)
}


midas_r.fit <- function (x) {

      args <- x$argmap_opt
      function.opt <- args$Ofunction
      args$Ofunction <- NULL
      if (!(function.opt %in% c("optim", "spg", "optimx", "lm",
                                "nls", "dry_run")))
            stop("The optimisation function is not in the supported functions list. Please see the midasr:::midas_r.fit code for the supported function list")
      if (function.opt == "optim" | function.opt == "spg") {
            args$par <- x$start_opt
            args$fn <- x$fn0
            if (x$use_gradient) {
                  args$gr <- x$gradient
            }
            opt <- try(do.call(function.opt, args), silent = TRUE)
            if (inherits(opt, "try-error")) {
                  stop("The optimisation algorithm of MIDAS regression failed with the following message:\n",
                       opt, "\nPlease try other starting values or a different optimisation function")
            }
            par <- opt$par
            names(par) <- names(coef(x))
            x$convergence <- opt$convergence
      }
      if (function.opt == "optimx") {
            args$par <- x$start_opt
            args$fn <- x$fn0
            if (x$use_gradient) {
                  args$gr <- x$gradient
            }
            opt <- try(do.call(function.opt, args), silent = TRUE)
            if (inherits(opt, "try-error")) {
                  stop("The optimisation algorithm of MIDAS regression failed with the following message:\n",
                       opt, "\nPlease try other starting values or a different optimisation function")
            }
            bmet <- which.min(opt$value)
            par <- as.numeric(opt[bmet, 1:length(args$par)])
            names(par) <- names(coef(x))
            x$convergence <- opt$convcode[bmet]
      }
      if (function.opt == "lm") {
            if (is.null(x$unrestricted))
                  stop("Not possible to estimate MIDAS model, more parameters than observations")
            par <- coef(x$unrestricted)
            names(par) <- names(coef(x))
            opt <- NULL
            x$convergence <- 0
      }
      if (function.opt == "nls") {
            rhs <- x$rhs
            if (x$use_gradient) {
                  orhs <- rhs
                  rhs <- function(p) {
                        res <- orhs(p)
                        attr(res, "gradient") <- x$model[, -1] %*% x$gradD(p)
                        res
                  }
            }
            y <- x$model[, 1]
            args$formula <- formula(y ~ rhs(p))
            args$start <- list(p = x$start_opt)
            browser()
            opt <- try(do.call("nls", args), silent = TRUE)
            if (inherits(opt, "try-error")) {
                  stop("The optimisation algorithm of MIDAS regression failed with the following message:\n",
                       opt, "\nPlease try other starting values or a different optimisation function")
            }
            par <- coef(opt)
            names(par) <- names(coef(x))
            x$convergence <- opt$convInfo$stopCode
      }
      if (function.opt == "dry_run") {
            opt <- NULL
            par <- x$start_opt
      }
      x$opt <- opt
      x$coefficients <- par
      names(par) <- NULL
      x$midas_coefficients <- x$gen_midas_coef(par)
      x$fitted.values <- as.vector(x$model[, -1] %*% x$midas_coefficients)
      x$residuals <- as.vector(x$model[, 1] - x$fitted.values)
      x
}
