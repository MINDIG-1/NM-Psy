vars <- ls()
rm(list = vars)

source("/home/ubuntu/NM-Psy/Codes/Utils.r")

stage <- "ModelTuning"

if (!file.exists(result_root)) dir.create(result_root)
result_path <- paste(result_root,stage, sep = "/")
if (!file.exists(result_path)) dir.create(result_path)
result_path <- paste(result_root,stage, "Distnpoly", sep = "/")
if (!file.exists(result_path)) dir.create(result_path)



dist_list <- c("ZAGA","ZAIG", "BEINF","PE","GT","GG","GB2","GB1","BCPE","BCT","exGAUS","JSU","NET","SEP1","SEP2","SEP3","SEP4","SHASH","SHASHo","ST1","ST2","ST3","ST4","ST5","TF")
dist_list <- c("SHASH", "GG")

npoly_list <- matrix(c(1,1,0,
                    2,1,0,
                    2,1,1,
                    2,2,1,
                    2,2,2,
                    3,1,0,
                    3,1,1,
                    3,2,1,
                    3,2,2,
                    3,3,0,
                    3,3,1,
                    3,3,2),
                    byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma","nu")))

npoly_list <- matrix(c(1,0,0,
                   1,1,0),
                   byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma","nu")))

n_c <- 50

df_all <- data.frame(f_band = character(), iFP = character(), DistFamily = character(),
                     AIC_Score = numeric(), BIC_Score = numeric(), mdlDev = numeric(), mu_formula=character(), si_formula= character(), nu_formula=character(), stringsAsFactors = FALSE)


for (fb in f_bands){

  data_HC_all <- subset(data_HC_all_fb, f_band == fb)

  df <- data.frame(f_band = character(), npoly = character(), DistFamily = character(),
                   AIC_Score = numeric(), BIC_Score = numeric(), mdlDev = numeric(), mu_formula=character(), si_formula= character(), nu_formula=character(),stringsAsFactors = FALSE)

  for (d in dist_list){

      for(np in 1:NROW(npoly_list) ) {

            if(npoly_list[np,"mu"]>0){
            mu_formula <- as.formula(paste("y ~ fp(age, ", npoly_list[np,"mu"], ")"))
            }else{
            mu_formula <- as.formula(paste("y ~ 1 "))
            }
            if(npoly_list[np,"sigma"]>0){
            sigma_formula <- as.formula(paste("~ fp(age, ", npoly_list[np,"sigma"], ")"))
            }else{
            sigma_formula <- as.formula(paste("~ 1 "))
            }
            if(npoly_list[np,"nu"]>0){
            nu_formula <- as.formula(paste("~ fp(age, ", npoly_list[np,"nu"], ")"))
            }else{
            nu_formula <- as.formula(paste("~ 1"))
            }

            tryCatch({
                print(paste(fb, d))
                model <- gamlss(mu_formula, sigma.formula = sigma_formula, nu.formula = nu_formula, family = d, data = data_HC_all, n.cyc = n_c)
                temp_df <- data.frame(f_band = fb, npoly=paste(as.character(npoly_list[np,]), collapse = ""),
                                    DistFamily = d, AIC_Score = model$aic, BIC_Score = model$sbc,mdlDev = model$G.deviance,
                                    mu_formula=deparse(mu_formula), si_formula= deparse(sigma_formula), nu_formula=deparse(nu_formula),
                                    stringsAsFactors = FALSE)
                df <- rbind(df, temp_df)
                df_all <- rbind(df_all, temp_df)

            }, error = function(e) {
                error_list <- conditionMessage(e)
            })
        }
    }
  
  sorted_BIC <- df[order(df$BIC_Score), ]
  df_all <- rbind(df_all, sorted_BIC)
  write.csv(sorted_BIC, paste(result_path,paste0("BIC_scores_", fb, ".csv"), sep = "/"))

}

