vars <- ls()
rm(list = vars)

source("/home/ubuntu/NM-Psy/Codes/Utils.r")


stage <- "ModelTraining"

result_path <- paste(result_root, stage, "Zs-Parallel", sep = "/") 
score_path <- paste(result_root, "ModelTuning", "CovEffect", sep = "/")
model_path <-paste(result_root, stage, "Models", sep = "/")
data_split_path <- paste(result_root, stage, "DataSplit", sep = "/")
if (!file.exists(result_path)) dir.create(result_path)


param_files <- list.files(score_path)
csv_files <- param_files[grep("^cov_result", param_files)]

dataset_list <- list.files(all_ch_path)

feature <- "Spec"  #or "FC"

if (feature == "FC"){

  df_fb <- read.csv(all_ch_path)
  df_fb <- df_fb[, -which(names(df_fb) %in% col_names_to_drop)]
  df_fb <- subset(df_fb, age < max_age)
  df_fb <- subset(df_fb, age > min_age)
}



for (fb in f_bands)
    {
    if (feature == "Spec"){
    df <- read.csv(paste(all_ch_path, dataset_list[grep(match(fb, f_bands), dataset_list)], sep = '/')) 
    df <- df[, -which(names(df) %in% col_names_to_drop)]
    df <- subset(df, age < max_age)
    df <- subset(df, age > min_age)

    }else if (feature == "FC"){

    df <- subset(df_fb, f_band == fb)
    }


    data_HC <- subset(df, group == "HC")
    data_asd <- subset(df, group == "ASD")
    data_adhd <- subset(df, group == "ADHD")
    data_anx <- subset(df, group == "Anxiety")
    data_ld <- subset(df, group == "Learning")

    #Downsampling ADHD, ANX, LD
    strata_columns <- c("age", "sex")
    data_adhd <- df_split(data_adhd,strata_columns, 0.5 )$df1
    data_ld <- df_split(data_ld,strata_columns, 0.99)$df1
    data_anx <- df_split(data_anx,strata_columns, 0.7)$df1


    #HC train test split
    data_HC_tr <- read.csv(paste(data_split_path, fb, "data_HC_tr.csv",sep = "/"))
    data_HC_te <- read.csv(paste(data_split_path, fb, "data_HC_te.csv",sep = "/"))


    temp_HC_tr <- dplyr::select(data_HC_tr, global_id,age, sex, dataset,group)
    temp_HC_te <-dplyr::select(data_HC_te, global_id,age, sex, dataset,group)
    temp_asd <-dplyr::select(data_asd, global_id,age, sex, dataset,group)
    temp_adhd <-dplyr::select(data_adhd, global_id,age, sex, dataset,group)
    temp_anx <-dplyr::select(data_anx, global_id,age, sex, dataset,group)
    temp_ld <-dplyr::select(data_ld, global_id,age, sex, dataset,group)

    # params <- read.csv(paste(score_path, csv_files[grep(fb, csv_files)], sep = "/"))

    

    for (j in seq_along(column_names)){

        col <- column_names[j]
        m <- readRDS(paste(model_path,fb, paste0("m",j,".rds"), sep = "/"))
        DistFam <- m$family[0]

        data_tr <- data_sub(data_HC_tr, col)
        data_te <- data_sub(data_HC_te, col)
        data_te_asd <- data_sub(data_asd, col)
        data_te_adhd <- data_sub(data_adhd, col)
        data_te_anx <- data_sub(data_anx, col)
        data_te_ld <- data_sub(data_ld, col)

        col_name <- paste0("zs_", col)
        temp_HC_tr[col_name] <- gZscore(m,data_tr)
        temp_HC_tr$group <- "HC_tr"
        temp_HC_te[col_name] <- gZscore(m,data_te)
        temp_HC_te$group <- "HC_te"
        temp_asd[col_name] <- gZscore(m,data_te_asd)

        temp_adhd[col_name] <- gZscore(m,data_te_adhd)
        temp_anx[col_name] <- gZscore(m,data_te_anx)
        temp_ld[col_name] <- gZscore(m,data_te_ld)

    }
    combined_data <- rbind(temp_HC_tr, temp_HC_te, temp_asd, temp_adhd, temp_anx, temp_ld)

    zs_save <- paste(result_path,paste0("Allch-", fb), sep = "/")
    if (!dir.exists(zs_save))
    {
      dir.create(zs_save)
    }
    write.csv(combined_data, paste(zs_save,"zsAllCh.csv", sep='/'))



}

