vars <- ls()
rm(list = vars)

source("/home/ubuntu/NM-Psy/Codes/Utils.r")

stage <- "ModelTuning"

score_path <- paste(result_root,stage, "CovEffect", sep = "/")
result_path <- paste(result_root,stage, "Bootstrap", sep = "/") 
csv_path <- paste(result_root,stage, "Bootstrap", "0_CSVs",sep = "/")  

if (!dir.exists(csv_path)) dir.create(csv_path)
if (!file.exists(result_path)) dir.create(result_path)

param_files <- list.files(score_path)
csv_files <- param_files[grep("^cov_result", param_files)]

for (fb in f_bands) {
  res_path <- paste(result_path,fb, sep ="/")
  if (!dir.exists(res_path))  dir.create(res_path)
  csv_path2 <- paste(csv_path,fb, sep ="/")
  if (!dir.exists(csv_path2))  dir.create(csv_path2)
}


num_resamples <- 1000

test_data <- data.frame(age = seq(min_age, max_age, by = 0.05) )

for (fb in f_bands){

    params <- read.csv(paste(score_path, csv_files[grep(fb, csv_files)], sep = "/"))
    mu_formula <- as.formula(params$mu_formula)
    si_formula <- as.formula(params$si_formula)
    nu_formula <- as.formula(params$nu_formula)
    DistFam <- params$DistFam


    data_HC_all <- subset(data_HC_all_fb, f_band == fb)

    num_observations <- nrow(data_HC_all)
    data_tr <- data_HC_all

    mu_values <- list()
    mu_cplx_b <- list()
    for (j in 1:num_resamples){
        print(j)
        indices <- sample(1:num_observations, replace = TRUE)
        resampled_df <- data_tr[indices, ]
        m <- try(gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula, family = DistFam, data = resampled_df, n.cyc = 50,  trace = TRUE), silent = TRUE)
        if (!inherits(m, "try-error")) {
            mu <- cent_cov(m,params$mu_formula, 50)
        }

        mu_values[[j]] <- mu
      
      }

    mu_matrix <- do.call(cbind, mu_values)
    mu_mat <- t(mu_matrix)
    confidence_interval <- apply(mu_mat, 2, function(col) quantile(col, c(0.025, 0.975)))

    m <- gamlss(mu_formula, sigma.formula = si_formula,nu.formula = nu_formula, family = DistFam, data = data_tr, n.cyc =50)

    mu <- cent_cov(m, params$mu_formula, 50)

    plot_data <- data.frame(x = test_data$age, mu = mu, lower = confidence_interval[1, ], upper = confidence_interval[2, ])

    p <- ggplot(plot_data) +
        geom_point(data = data_tr, aes(x = age, y = y), color = "gray", alpha = 0.5) +  
        geom_ribbon(aes(x = x, ymin = lower, ymax = upper), fill = "gray", alpha = 0.8)+
        geom_line(aes(x = x, y = mu), color = "#8B285A", size = 2) + 
        geom_line(aes(x = x, y= lower), color = "#F6C0DC", linetype = "dotted", size = 1) +  
        geom_line(aes(x = x, y= upper), color = "#F6C0DC", linetype = "dotted", size = 1) +
        labs(x = "Age", y = paste0("Rel Power ",fb)) +  # Labels for axes 
        ylim(0, 1)+
        theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), panel.background = element_rect(fill = "white"))

    print(p)



    cpname <- paste(result_path,fb, paste0("bpMedian",fb ,".png"), sep = "/")
    ggsave(filename = cpname, plot = p, width = 10, height = 8)
    write.csv(plot_data, paste(csv_path,fb,"bpMedian.csv", sep='/'), row.names = FALSE)


    }
    