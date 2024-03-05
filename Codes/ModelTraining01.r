vars <- ls()
rm(list = vars)

source("/home/ubuntu/NM-Psy/Codes/Utils.r")


stage <- "ModelTraining"

score_path <- paste(result_root, "ModelTuning", "CovEffect", sep = "/")
result_path <- paste(result_root,stage, sep = "/")
if (!file.exists(result_path)) dir.create(result_path)
result_path <- paste(result_root, stage, "Zs-Parallel", sep = "/") 
model_path <-paste(result_root, stage, "Models", sep = "/")
data_split_path <- paste(result_root, stage, "DataSplit", sep = "/")
if (!file.exists(result_path)) dir.create(result_path)

if (!file.exists(data_split_path)) dir.create(data_split_path)
if (!file.exists(model_path)) dir.create(model_path)
if (!file.exists(result_path)) dir.create(result_path)


registerDoParallel(cores = detectCores()-1)

set.seed(42) 

folder_l <- c("DataSplit", "Models")
for (folder_n in folder_l){
for (fb in f_bands){
  dir_n <- paste(result_root,stage, folder_n, fb, sep = "/")
  if (!file.exists(dir_n)) dir.create(dir_n)
}}

param_files <- list.files(score_path)
csv_files <- param_files[grep("^cov_result", param_files)]


dataset_list <- list.files(all_ch_path)

model_not_trained <- list()

n_cyc <- 300
 

for (fb in f_bands)
{
    gDev_df_list <- list()

    df = read.csv(paste(all_ch_path, dataset_list[grep(match(fb, f_bands), dataset_list)], sep = '/'))  
    df <- df[, -which(names(df) %in% col_names_to_drop)]
    df <- subset(df, age < max_age)
    df <- subset(df, age > min_age)

    data_HC <- subset(df, group == "HC")


    strata_columns <- c("age", "sex", "dataset")
    data_HC<- df_split(data_HC,strata_columns, 0.8)
    data_HC_tr <- data_HC$df1
    data_HC_te <- data_HC$df2


    write.csv(data_HC_tr, paste(data_split_path, fb, "data_HC_tr.csv",sep = "/"), row.names = FALSE)
    write.csv(data_HC_te, paste(data_split_path, fb, "data_HC_te.csv",sep = "/"), row.names = FALSE)


    params <- read.csv(paste(score_path, csv_files[grep(fb, csv_files)], sep = "/"))
    mu_formula <- as.formula(params$mu_formula)
    si_formula <- as.formula(params$si_formula)
    nu_formula <- as.formula(params$nu_formula)

    DistFam <- params$DistFam

    models_list <- list()

    results <- foreach(j = seq_along(column_names), .combine = c, .export = c("column_names", "data_HC_tr", "data_HC_te")) %dopar% {


    col <- column_names[j]

    data_tr <- data_sub(data_HC_tr, col)
    data_te <- data_sub(data_HC_te,col)

    m <- try(gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula, family = DistFam, data = data_tr, n.cyc = n_cyc), silent = TRUE)


    if (inherits(m, "try-error")) {
        m <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula, family = "SHASH", data = data_tr, n.cyc = n_cyc)
        model_not_trained <- append(model_not_trained, list(list(fb, col, "SHASH")))
    }
    saveRDS(m, paste(model_path,fb, paste0("m",j,".rds"), sep = "/"))


    }

  
  stopImplicitCluster()

}

  