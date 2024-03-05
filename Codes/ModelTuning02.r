vars <- ls()
rm(list = vars)

source("/home/ubuntu/NM-Psy/Codes/Utils.r")

stage <- "ModelTuning"

score_path <- paste(result_root, stage, "Distnpoly", sep = "/")
result_path <- paste(result_root, stage, "CovEffect", sep = "/") 
if (!file.exists(result_path)) dir.create(result_path)

param_files <- list.files(score_path)
bic_files <- param_files[grep("^BIC", param_files)]
bic_files <-  bic_files[bic_files != "BIC_scores_all.csv"]



phase2 <- TRUE

for (fb in f_bands){

    cov_results_df <- data.frame()
    cov_sex_df <- data.frame()

    data_HC_all <- subset(data_HC_all_fb, f_band == fb)

    data_tr <- data_HC_all
  
    params <- read.csv(paste(score_path, bic_files[grep(fb,bic_files)], sep="/"))
    params <- params[1,]
    mu_formula <- as.formula(params$mu_formula)
    si_formula <- as.formula(params$si_formula)
    nu_formula <- as.formula(params$nu_formula)
    DistFam <- params$DistFamily
    # formulas <- c("mu_formula", "si_formula", "nu_formula")



  #------------------------------------------------Covariate :: Sex :: Fixed
    cov <- "factor(sex)"
    aic_sex <- cov_test(mu_formula,si_formula, nu_formula,DistFam, cov, data_tr, fb)
    sorted_aic_sex <- aic_sex[order(aic_sex$AIC_score),]
    sorted_bic_sex <- aic_sex[order(aic_sex$BIC_score),]

    sorted_bic_sex <- sorted_bic_sex [1,]

    cov_sex_df <- rbind(cov_sex_df,sorted_bic_sex)
    cov_file_name <- paste0("cov_sex_",fb,".csv")
    write.csv(cov_sex_df, paste(result_path, cov_file_name, sep="/"))


  #------------------------------------------------Covariate :: Site :: Fixed

    mu_formula <- as.formula(sorted_bic_sex$mu_formula)
    si_formula <- as.formula(sorted_bic_sex$si_formula)
    nu_formula <- as.formula(sorted_bic_sex$nu_formula)

    cov <- "factor(dataset)"
    aic_site <- cov_test(mu_formula,si_formula, nu_formula,DistFam, cov, data_tr, fb)
    sorted_aic_site <- aic_site[order(aic_site$AIC_score),]
    sorted_bic_site <- aic_site[order(aic_site$BIC_score),]

    sorted_bic_site <- sorted_bic_site [1,]


  # ------------------------------------------------Covariate :: Site :: Random


    mu_formula <- as.formula(sorted_bic_site$mu_formula)
    si_formula <- as.formula(sorted_bic_site$si_formula)
    nu_formula <- as.formula(sorted_bic_site$nu_formula)
    

    cov <- "random(factor(site))"
    aic_r_site <- cov_test(mu_formula,si_formula, nu_formula, DistFam, cov, data_tr, fb)
    sorted_aic_r_site <- aic_r_site[order(aic_r_site$AIC_score),]
    sorted_bic_r_site <- aic_r_site[order(aic_r_site$BIC_score),]

    sorted_bic_r_site <- sorted_bic_r_site [1,]

  #-------------------------------------------------Final model
    cov_results_df <- rbind(cov_results_df,sorted_bic_r_site)


    if (phase2)
    {

        mu_formula <- as.formula(cov_results_df$mu_formula)
        si_formula <- as.formula(cov_results_df$si_formula)
        nu_formula <- as.formula(cov_results_df$nu_formula)

        tryCatch ({
            for (d in dist_list)
            {
                model <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula,family = d, data = data_tr, n.cyc = 50)
                if (model$sbc < cov_results_df$BIC_score)
                {
                    cov_results_df$DistFam <- model$family[1]
                }
                }
            }, error = function(e) {

            })
    }

    cov_file_name <- paste0("cov_result_",fb,".csv")
    write.csv(cov_results_df, paste(result_path, cov_file_name, sep="/"))

    

  

}

