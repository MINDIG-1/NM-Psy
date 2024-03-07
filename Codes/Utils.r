library(gamlss)
library(ggplot2)
library(dplyr)
library(splitTools)
library(caret)
library(doParallel)
library(foreach)


#-------------------------------Feature Specific
feature <- "FC"  #or "FC"  "Spec"

if (feature == "Spec"){

  df_path <- "/home/ubuntu/MEEG-normative-modeling/NModel_tr_te/Featuresdf/df_rel_pwr_avgch_allsites_EO.csv"
  all_ch_path <- "/home/ubuntu/MEEG-normative-modeling/NModel_tr_te/AllchAllsites"
  column_names <- paste0("X", 1:124)
  col_names_to_drop <- c("X125","X126", "X127", "X128", "subject_id", "handedness", "scores")


} else if (feature == "FC"){
  df_path <- "/home/ubuntu/MEEG-normative-modeling/NModel_FC/Featuresdf/df_FC_avg.csv"
  all_ch_path <- "/home/ubuntu/MEEG-normative-modeling/NModel_FC/AllchAllsites"
  col_names_to_drop <- c("subject_id", "handedness", "scores")
  column_names <- paste0("X", 0:2277)

}


#------------------------------Shared
result_root <- "/home/ubuntu/NM-Psy/Results"
min_age <- 5
max_age <- 18

data_all <- read.csv(df_path)
data_all <- subset(data_all, age < max_age)
data_all <- subset(data_all, age > min_age)
data_HC_all_fb <- subset(data_all, group == "HC")
if (c("handedness", "scores") %in% colnames(data_HC_all_fb))  data_HC_all_fb <- subset(data_HC_all_fb, select = -c(handedness, scores))

f_bands <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")

dist_list <- c("PE","GT","GG","GB2","BCPE","BCT","exGAUS","JSU","SEP1","SEP2","SEP3","SEP4","SHASH","SHASHo","ST3","ST4")
# dist_list <- c("SHASH", "GG")

site_list <- c("ABCCT", "MIPDB",  "lausanneASD" ,"femaleASD"  , "HBN"  )
sex_list <- c("M","F")



cov_test <- function(m_formula, si_formula, nu_formula, DistFam, cov, data_tr, fb){
  formulas <- c("mu_formula", "si_formula", "nu_formula")
  i <- 1
  m0 <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula, family = DistFam, data = data_tr, n.cyc = 50)
  
  score_df <- data.frame(m_i = i-1, AIC_score = m0$aic, BIC_score = m0$sbc, mu_formula = deparse(mu_formula), si_formula=deparse(si_formula), nu_formula=deparse(nu_formula), DistFam=DistFam)
                         
  for (form in formulas){
    tryCatch({
        eval(parse(text = paste(form, "<- update(", form, ", ~ . + ",cov,")")))
        m <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula, family = DistFam, data = data_tr, n.cyc = 50)
        
        
        temp_df <- data.frame(m_i= i, AIC_score = m$aic, BIC_score = m$sbc, mu_formula = deparse(mu_formula), si_formula=deparse(si_formula), nu_formula=deparse(nu_formula), DistFam=DistFam )

        
        score_df <- rbind(score_df,temp_df)

    }, error = function(e) {
        warning(paste("Error in iteration", i," in ",cov, ":", conditionMessage(e)))
        temp_df <- data.frame(m_i= 1e10, AIC_score = 1e10, BIC_score = 1e10, mu_formula = deparse(mu_formula), si_formula=deparse(si_formula), nu_formula=deparse(nu_formula), DistFam=DistFam,
                            expv_test= 0, msll_test=0, smse_test=0 )

    })
    i <- i + 1
    
  }
  return(score_df)

}

cent_cov <- function(m, mu_formula, cent ){
    if (grepl('sex', mu_formula) && grepl('dataset', mu_formula)){
        mu <- 0
        for (q in seq_along(site_list)){
            test_data['dataset'] <- site_list[q]
            for ( k in seq_along(sex_list)){
                test_data['sex'] <- sex_list[k]
                mu_temp <- calcent(m,cent,test_data)
                mu <- mu + mu_temp
            }}
            mu <- mu/(length(site_list) * length(sex_list))

    } else if (grepl('sex', mu_formula)){
            mu <- 0
            for ( k in seq_along(sex_list)){
                test_data['sex'] <- sex_list[k]
                mu_temp <- calcent(m,cent,test_data)
                mu <- mu + mu_temp
            }
            mu <- mu/length(sex_list)

    } else if (grepl('dataset', mu_formula)) {
        mu <- 0
        for ( k in seq_along(site_list)){
            test_data['dataset'] <- site_list[k]
            mu_temp <- calcent(m,cent,test_data)
            mu <- mu + mu_temp
          }
          mu <- mu/length(site_list)
    } else {
        mu <- calcent(m,cent,test_data)
    }
    return (mu)
}

df_split <- function(df, strata_columns, p){
    ir <- df[strata_columns]
    y <- multi_strata(ir)
    inds <- partition(y, p = c(keep = p, notkeep = 1-p), split_into_list = FALSE)
    df_keep <- df[inds=='keep', ]
    df_notkeep <- df[inds=='notkeep', ]
    return(list(df1 = df_keep, df2 = df_notkeep))

}

data_sub <- function (df, col) {
    df <- dplyr::select(df, global_id,age, sex, dataset,group, all_of(col))
    colnames(df)[colnames(df) == col] <- "y"
    return (df)
}

gZscore <- function (obj, newdata)
  {
  predData <- newdata[, -which(names(newdata) == "y")]

  if ("mu"%in%obj$parameters )
      {if ( is.null(obj$mu.fix))
        mu <- predict(obj,what = "mu",   newdata = predData, type = "response")
    else if (obj$mu.fix==TRUE) mu <- rep(fitted(obj, "mu")[1], length(xvar))
    }
  if ("sigma"%in%obj$parameters)
    {if (is.null(obj$sigma.fix))
        sigma <- predict(obj,what = "sigma",   newdata = predData, type = "response")
    else if (obj$sigma.fix==TRUE) sigma <- rep(fitted(obj, "sigma")[1], length(xvar))
    }
  if ("nu"%in%obj$parameters )
    { if  (is.null(obj$nu.fix))
        nu <- predict(obj,what = "nu",   newdata = predData, type = "response")
    else if (obj$nu.fix==TRUE) nu <- rep(fitted(obj, "nu")[1], length(xvar))
    }
  if ("tau"%in%obj$parameters )
    { if (is.null(obj$tau.fix))
        tau <- predict(obj,what = "tau",   newdata = predData, type = "response")
    else if (obj$tau.fix==TRUE) tau <- rep(fitted(obj, "tau")[1], length(xvar))
    }


    yval <- newdata$y
    lpar <- length(obj$parameters)
    fname <- obj$family[1]
    pfun <- paste("p",fname,sep="")

    if(lpar==1)
            {newcall <-call(pfun,yval,mu=mu) }
    else if(lpar==2)
            {newcall <-call(pfun,yval,mu=mu,sigma=sigma) }
    else if(lpar==3)
            {newcall <-call(pfun,yval,mu=mu,sigma=sigma,nu=nu) }
    else
            {newcall <-call(pfun,yval,mu=mu,sigma=sigma,nu=nu,tau=tau) }

    cdf <- eval(newcall)
    rqres <- qnorm(cdf)
    return (rqres)
}


calcent <- function(obj, cent, predData) {

  if ("mu"%in%obj$parameters )
      {if ( is.null(obj$mu.fix))
        mu <- predict(obj,what = "mu",   newdata = predData, type = "response")
    else if (obj$mu.fix==TRUE) mu <- rep(fitted(obj, "mu")[1], length(xvar))
    }
  if ("sigma"%in%obj$parameters)
    {if (is.null(obj$sigma.fix))
        sigma <- predict(obj,what = "sigma",   newdata = predData, type = "response")
    else if (obj$sigma.fix==TRUE) sigma <- rep(fitted(obj, "sigma")[1], length(xvar))
    }
  if ("nu"%in%obj$parameters )
    { if  (is.null(obj$nu.fix))
        nu <- predict(obj,what = "nu",   newdata = predData, type = "response")
    else if (obj$nu.fix==TRUE) nu <- rep(fitted(obj, "nu")[1], length(xvar))
    }
  if ("tau"%in%obj$parameters )
    { if (is.null(obj$tau.fix))
        tau <- predict(obj,what = "tau",   newdata = predData, type = "response")
    else if (obj$tau.fix==TRUE) tau <- rep(fitted(obj, "tau")[1], length(xvar))
    }


    lpar <- length(obj$parameters)
      fname <- obj$family[1]
      qfun <- paste("q",fname,sep="")
  xvar <- predData$age
      o <- order(xvar)
        mat <- xvar[o]
        cent <- cent
          for (var in cent) {
              if (lpar == 1) {        
                  newcall <- call(qfun, var/100, mu = mu[o])
              }
              else if (lpar == 2) {
                  newcall <- call(qfun, var/100, mu = mu[o], 
                    sigma = sigma[o])
              }
              else if (lpar == 3) {
                  newcall <- call(qfun, var/100, mu = mu[o], 
                    sigma = sigma[o], nu = nu[o])
              }
              else {
                  newcall <- call(qfun, var/100, mu = mu[o], 
                    sigma = sigma[o], nu = nu[o], 
                    tau = tau[o])
              }
              ll <- eval(newcall)
              mat <- ll
          }

          return(mat)
}









