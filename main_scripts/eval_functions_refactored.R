# +
## Helper Functions
# -

# Function to perform adjusted PC calculation
perform_adj_PC <- function(dat, score_name, age_sex) {
  formula <- if (age_sex) {
    paste(score_name, "~ age + is_male + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  } else {
    paste(score_name, "~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  }
  mod <- lm(as.formula(formula), data = dat)
  return(residuals(mod))
}

R2ObservedLinear <- function(data, phecode, score_num, n, ncase, ncont, covs){
  if(is.na(score_num)){
      lmv <- lm(as.formula(paste(phecode, paste(c(covs), collapse=" + "), sep=" ~ ")), data = data)
  }else{
      score_char <- paste0("scale(SCORE", score_num, "_SUM)")
      lmv <- lm(as.formula(paste(phecode, paste(c(score_char, covs), collapse=" + "), sep=" ~ ")), data = data)
  }
  
  R2_lin <- var(lmv$fitted.values)/(ncase/n*ncont/n)  
  return(R2_lin)}

R2Liability<- function(data, phecode,score_num, covs){
  n <- dim(data)[1]
  ncase <- sum(data[[phecode]] == T, na.rm = T)
  ncont <- sum(data[[phecode]] == F, na.rm = T)
  P = ncase/n
  thresh <- -qnorm(P,0,1) #The value along normal distribution that captures the proportion in the population with condition 
  zval <- dnorm(thresh) #the y coordinate of this threshold, point on density curve
  R2_lin <- R2ObservedLinear(data, phecode,score_num,n,ncase,ncont, covs)
  R2_obs <- R2_lin*(P*(1-P))/(zval^2)
  return(c(R2_lin, R2_obs))
}

inormal <- function(x){qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))}

getVals <- function(cormat, outcome, measures){
  return(cormat[outcome, measures])
}

# Define a helper function to perform linear regression and write beta lines
perform_lm_and_write_beta <- function(dat_anc, score_col, phecode, pop, beta_file) {
  lm_model <- lm(as.formula(paste("mean1", paste(c(paste0("scale(", score_col, ")")), collapse = " + "), sep = " ~ ")), data = dat_anc)
  write_beta_line(lm_model, score_col, phecode, pop, beta_file)
  return(summary(lm_model)$r.squared)
}

write_beta_line <- function(model, score_col, phecode, pop, file){
  write(paste(c(phecode, pop, score_col, round(summary(model)$coefficients[2, c(1, 2,3)],4), format.pval(summary(model)$coefficients[2, 4], digits = 2)), collapse = " ") , file = file, append = T)
}


combine_res <- function(dir, pattern){
    allfiles <- list.files(path = dir, pattern = pattern)
    readdata <- function(fn){ dt_temp <- fread(paste0(dir, fn)); return(dt_temp)}
    mylist <- lapply(allfiles, readdata)
    mydata <- do.call('rbind',mylist) 
    return(mydata)
}

combine_res_prscsx <- function(dir, pattern) 
{
    allfiles <- list.files(path = dir, pattern = pattern)
    readdata <- function(fn) {
        dt_temp <- fread(paste0(dir, fn))
        return(dt_temp)
    }
    mylist <- lapply(allfiles, readdata)
    mylist <- lapply(mylist, add_colnames_prscsx)
    mydata <- rbindlist(mylist, fill = TRUE)
    return(mydata)
}
add_colnames_prscsx <- function(dat){
    if(ncol(dat) == 10){
        colnames(dat) <-  c("phecode", "anc", "cases", "controls", "baseline", "AFRx", "EURx", "AMRx", "EASx", "prscsx")
    }else{
        colnames(dat) <-  c("phecode", "anc", "cases", "controls", "baseline", "AFRx", "EURx", "AMRx", "prscsx")
    }
   return(dat)
}

#confirmed is equivalent to before, will use this!
# Function to compare R2 from result file with adjusted PC
compare_scores_PCadj_return_list <- function(score_data, phen_data, covars, age_sex = FALSE) {
  score_num <- length(colnames(score_data)) - 4
  print(paste("Score Number:", score_num))
  dat <- phen_data %>% 
    left_join(covars) %>% 
    left_join(score_data, by = c(person_id = "IID")) %>% 
    filter(unrel == 1)
  print(paste("Data Dimensions:", dim(dat)))
  
  res <- list()
  for (pop in c("afr", "eur", "amr")) {
    print(paste("Population:", pop))
    dat_anc <- dat %>% filter(ancestry_pred_other == pop) %>% na.omit() 
    dat_anc$mean <- inormal(dat_anc$mean)
    dat_anc$mean1 <- perform_adj_PC(dat_anc, "mean", age_sex)
    scores <- lapply(1:score_num, function(i) perform_adj_PC(dat_anc, paste0("SCORE", i, "_SUM"), age_sex))
    res[[pop]] <- cbind(dat_anc$mean1, bind_cols(scores))
  }
  return(res)
}

#confirmed is equivalent to before, will use this!
make_outcome_list <- function (score_name, phen_name, age_sex = FALSE){
  ldpred <- fread(paste0("computed_scores/", score_name, "_LDpred.txt"))
  prscs <- fread(paste0("computed_scores/", score_name, "_PRScs.txt"))
  prscsx <- fread(paste0("computed_scores/", score_name, "_PRScsx.txt"))
  covars <- as.data.frame(fread("AoU_98K_covariates.tsv"))
  phen <- fread(paste0("outcomes/", phen_name, ".txt"))
    
  ldp_list <- compare_scores_PCadj_return_list(ldpred, phen, covars, age_sex)
  cs_list <- compare_scores_PCadj_return_list(prscs, phen, covars, age_sex)
  x_list <- compare_scores_PCadj_return_list(prscsx, phen,  covars, age_sex)
  alist <- Map(cbind, cs_list, ldp_list, x_list)
    
  new_col_name <- c("mean", paste0("prscs_score", 1:6), "mean2", 
                    paste0(c("ldp_inf", "ldp_auto"), rep(1:6, each = 2)), 
                    "mean3", paste0("x_score", 1:5))
  alist <- lapply(alist, setNames, nm = new_col_name)
  cormat_list <- lapply(alist, cor)
  cormat_df <- as.data.frame(do.call(rbind, lapply(cormat_list, 
                                                   getVals, outcome = "mean", measures = c(paste0("prscs_score", 1:6), 
                                                                                           paste0("ldp_auto", 1:6), 
                                                                                           "x_score5"))))
  cormat_df$anc <- rownames(cormat_df)
  colnames(cormat_df) <- c(prscs_names, ldpred_names, x_names, 
                           "anc")
  cormat_df$outcome <- phen_name
  return(list(cormat_df, cormat_list, alist))
}

##to evaluate scores for continous traits
compare_method_PCadj <- function(dir, covars, score_name, phen_name, method, age_sex = FALSE, sire = NULL) {
  phecode <- score_name
  
  score_data <- fread(paste0("computed_scores/", score_name, "_", method, ".txt"))
  phen_data <- fread(paste0("outcomes/", phen_name, ".txt"))
  
  score_num <- length(colnames(score_data)) - 4
  print(score_num)
  
  dat <- phen_data %>% 
    left_join(covars) %>% 
    left_join(score_data, by = c(person_id = "IID"))
  
  if (!is.null(sire)) {
    dat <- dat %>%
      left_join(sire) %>%
      filter(unrel == 1)
  } else {
    dat <- dat %>%
      filter(unrel == 1)
  }
  
  grouping_name <- ifelse(is.null(sire), "ancestry", "sire")
  adj_level <- ifelse(age_sex == FALSE, "PC-adj", "PCagesex-adj")
  R2_file <- paste0(dir, phecode,  "_", method, "_", grouping_name, "_R2_", adj_level, "_unrel.txt")
  beta_file <- paste0(dir, phecode,  "_", method, "_", grouping_name, "_beta_", adj_level, "_unrel.txt")
  
  grouping_var <- ifelse(is.null(sire), "ancestry_pred_other", "sire")
  for (pop in unique(dat[[grouping_var]])[!is.na(unique(dat[[grouping_var]]))] ) {
    dat_anc <- dat[dat[[grouping_var]] == pop, ]
    dat_anc <- na.omit(dat_anc)
    dat_anc$mean <- inormal(dat_anc$mean)
    for (i in 1:score_num) {
      score_col <- paste0("SCORE", i, "_SUM")
      dat_anc[[paste0("score", i)]] <- perform_adj_PC(dat_anc, score_col, age_sex)
    }
    dat_anc$mean1 <- perform_adj_PC(dat_anc, "mean", age_sex)
    
    r_squared_values <- numeric(score_num)
    for (i in 1:score_num) {
      score_col <- paste0("score", i)
      r_squared_values[i] <- perform_lm_and_write_beta(dat_anc, score_col, phecode, pop, beta_file)
    }
    write(paste(phecode, pop, dim(dat_anc)[1], paste(r_squared_values, collapse = ";"), sep = ";"), file = R2_file, append = TRUE)
  }
}

##to evaluate scores for continous traits
compare_method_PCadj_boot <- function(dir, covars, score_name, phen_name, method, bootstrap = FALSE, age_sex = FALSE) {
  phecode <- score_name
  
  score_data <- fread(paste0("computed_scores/", score_name, "_", method, ".txt"))
  phen_data <- fread(paste0("outcomes/", phen_name, ".txt"))
  
  score_num <- length(colnames(score_data)) - 4
  print(score_num)
  
  dat <- phen_data %>% 
    left_join(covars) %>% 
    left_join(score_data, by = c(person_id = "IID"))
  
  dat <- dat %>%
      filter(unrel == 1)
  
  grouping_name <- "ancestry"
  adj_level <- ifelse(age_sex == FALSE, "PC-adj", "PCagesex-adj")
  R2_file <- paste0(dir, phecode,  "_", method, "_", grouping_name, "_R2_", adj_level, "_unrel.txt")
  beta_file <- paste0(dir, phecode,  "_", method, "_", grouping_name, "_beta_", adj_level, "_unrel.txt")
  boot_file <- paste0(dir, phecode,  "_", method, "_", grouping_name, "_boot_", adj_level, "_unrel.txt")
  
  grouping_var <- "ancestry_pred_other"
  for (pop in c("eur", "afr", "amr")) {
    dat_anc <- dat[dat[[grouping_var]] == pop, ]
    dat_anc <- na.omit(dat_anc)
    dat_anc$mean <- inormal(dat_anc$mean)
    for (i in 1:score_num) {
      score_col <- paste0("SCORE", i, "_SUM")
      dat_anc[[paste0("score", i)]] <- perform_adj_PC(dat_anc, score_col, age_sex)
    }
    dat_anc$mean1 <- perform_adj_PC(dat_anc, "mean", age_sex)
    
    r_squared_values <- numeric(score_num)
    for (i in 1:score_num) {
      score_col <- paste0("score", i)
      r_squared_values[i] <- perform_lm_and_write_beta(dat_anc, score_col, phecode, pop, beta_file)
      if(bootstrap == TRUE){
          boot_res <- boot(dat_anc, extract_r_squared, 
                           R = 1000, 
                           score_col = score_col, 
                           phecode = phecode)
          boot_dat <- paste(round(boot_res$t0, 5), 
                            paste0(round(quantile(boot_res$t, probs=c(.025, .975)), 5), 
                                   collapse = ";")  , sep = ";")
          write(paste(phecode, pop, dim(dat_anc)[1], score_col, boot_dat, sep = ";"), file = boot_file, append = TRUE)
          
      }
    }
    write(paste(phecode, pop, dim(dat_anc)[1], paste(r_squared_values, collapse = ";"), sep = ";"), file = R2_file, append = TRUE)
  }
}

# +
perform_lm <- function (dat_anc, score_col, phecode) 
{
    lm_model <- lm(as.formula(paste("mean1", paste(c(paste0("scale(", 
        score_col, ")")), collapse = " + "), sep = " ~ ")), data = dat_anc)
    return(summary(lm_model)$r.squared)
}

extract_r_squared <- function(data, indices, score_col, phecode) {
    d <- data[indices, ]  # Resample data
    perform_lm(d, score_col, phecode)  # Compute R-squared
}
# -

compare_method_popprev <- function(phecode, dir, method, age_sex = FALSE, sire = NULL) {
  print(phecode)
  score_data <- fread(paste0("computed_scores/", phecode, "_", method, ".txt"))
  score_num <- length(colnames(score_data)) - 4
  
  dat <- icd %>%
    select("person_id", all_of(phecode)) %>%
    left_join(covars) %>%
    left_join(score_data, by = c(person_id = "IID")) 
  
  if (!is.null(sire)) {
    dat <- dat %>%
      left_join(sire) %>%
      filter(unrel == 1)
  } else {
    dat <- dat %>%
      filter(unrel == 1)
  }
  
  grouping_var <- ifelse(is.null(sire), "ancestry_pred_other", "sire")
  grouping_name <- ifelse(is.null(sire), "ancestry", "sire")
  if (age_sex) {
    covs <- c("age", "is_male", paste0("PC", 1:10))
    adj_level <- "PCagesex-adj"
  } else {
    covs <- paste0("PC", 1:10)
    adj_level <- "PC-adj"
  }
  
  cstat_file <- paste0(dir, phecode,  "_", method, "_", grouping_name, "_cstat_", adj_level, "_unrel.txt")
  beta_file <- paste0(dir, phecode,  "_", method, "_", grouping_name, "_beta_", adj_level, "_unrel.txt")
  R2liab_file <- paste0(dir, phecode,  "_", method, "_", grouping_name, "_R2liab_", adj_level, "_unrel.txt") 
  print(cstat_file)
  
  for (pop in unique(dat[[grouping_var]])[!is.na(unique(dat[[grouping_var]]))]) {
    dat_anc <- dat[dat[[grouping_var]] == pop, ]
    dat_anc <- na.omit(dat_anc)
    ncases <- sum(dat_anc[[phecode]] == T, na.rm = T)
    ncontrols <- sum(dat_anc[[phecode]] == F, na.rm = T)
    if(ncases > 100){
        baseline <- glm(as.formula(paste(phecode, paste0(covs, collapse = "+"), sep = " ~ ")), family = "binomial", data = dat_anc)
        baseline_R2liab <- R2Liability(dat_anc, phecode, NA, covs)
    
        R2obs_values <- numeric(score_num)
        R2liab_values <- numeric(score_num)
        Cstat_values <- numeric(score_num)
        for (i in 1:score_num) {
          score_col <- paste0("SCORE", i, "_SUM")
          model <- glm(as.formula(paste(phecode, paste(c(paste0("scale(", score_col, ")"), covs), collapse = " + "), sep = " ~ ")), family = "binomial", data = dat_anc)
      
          write_beta_line(model, i, phecode, pop, beta_file)
          R2liab <- R2Liability(dat_anc, phecode, i, covs)
          R2obs_values[i] <- R2liab[1]
          R2liab_values[i] <- R2liab[2]
          Cstat_values[i] <- Cstat(model)
        }

        write(paste(phecode, pop, ncases, ncontrols, baseline_R2liab[1], paste(R2obs_values, collapse = ";"), sep = ";" ), file = R2liab_file, append = TRUE)
        write(paste(phecode, pop, ncases, ncontrols, baseline_R2liab[2], paste(R2liab_values, collapse = ";"), sep = ";"  ), file = R2liab_file, append = TRUE)
        write(paste(phecode, pop, ncases, ncontrols, Cstat(baseline), paste(Cstat_values, collapse = ";"),sep = ";"  ), file = cstat_file, append = TRUE)
    }
  }
}

compare_method_popprev_boot <- function (phecode, dir, method, age_sex = FALSE, sire = NULL, bootstrap = FALSE) 
{
    print(phecode)
    score_data <- fread(paste0("computed_scores/", phecode, "_", 
        method, ".txt"))
    score_num <- length(colnames(score_data)) - 4
    dat <- icd %>% select("person_id", all_of(phecode)) %>% left_join(covars) %>% 
        left_join(score_data, by = c(person_id = "IID"))
    if (!is.null(sire)) {
        dat <- dat %>% left_join(sire) %>% filter(unrel == 1)
    }
    else {
        dat <- dat %>% filter(unrel == 1)
    }
    grouping_var <- "ancestry_pred_other"
    grouping_name <- "ancestry"
    if (age_sex) {
        covs <- c("age", "is_male", paste0("PC", 1:10))
        adj_level <- "PCagesex-adj"
    }
    else {
        covs <- paste0("PC", 1:10)
        adj_level <- "PC-adj"
    }
    cstat_file <- paste0(dir, phecode, "_", method, "_", grouping_name, 
        "_cstat_", adj_level, "_unrel.txt")
    beta_file <- paste0(dir, phecode, "_", method, "_", grouping_name, 
        "_beta_", adj_level, "_unrel.txt")
    R2liab_file <- paste0(dir, phecode, "_", method, "_", grouping_name, 
        "_R2liab_", adj_level, "_unrel.txt")
    R2liab_boot_file <- paste0(dir, phecode, "_", method, "_", grouping_name, 
        "_R2liab_boot_", adj_level, "_unrel.txt")

    for (pop in c("eur", "afr", "amr")) {
        dat_anc <- dat[dat[[grouping_var]] == pop, ]
        dat_anc <- na.omit(dat_anc)
        ncases <- sum(dat_anc[[phecode]] == T, na.rm = T)
        ncontrols <- sum(dat_anc[[phecode]] == F, na.rm = T)
        if (ncases > 100) {
            baseline <- glm(as.formula(paste(phecode, paste0(covs, 
                collapse = "+"), sep = " ~ ")), family = "binomial", 
                data = dat_anc)
            baseline_R2liab <- R2Liability(dat_anc, phecode, 
                NA, covs)
            R2obs_values <- numeric(score_num)
            R2liab_values <- numeric(score_num)
            Cstat_values <- numeric(score_num)
            for (i in 1:score_num) {
                score_col <- paste0("SCORE", i, "_SUM")
                model <- glm(as.formula(paste(phecode, paste(c(paste0("scale(", 
                  score_col, ")"), covs), collapse = " + "), 
                  sep = " ~ ")), family = "binomial", data = dat_anc)
                write_beta_line(model, i, phecode, pop, beta_file)
                R2liab <- R2Liability(dat_anc, phecode, i, covs)
                R2obs_values[i] <- R2liab[1]
                R2liab_values[i] <- R2liab[2]
                Cstat_values[i] <- Cstat(model)
                if(bootstrap == TRUE){
                    boot_res <- boot(dat_anc, extract_R2diff, 
                                     R = 1000,
                                     phecode = phecode,
                                     score_num = i,
                                     covs = covs)
                    boot_dat <- paste(round(boot_res$t0, 5), 
                                      paste0(round(quantile(boot_res$t, probs=c(.025, .975)), 5), 
                                             collapse = ";")  , sep = ";")
                    write(paste(phecode, pop, ncases, ncontrols, score_col, boot_dat, sep = ";"), file = R2liab_boot_file, append = TRUE)
          
      }
            }
            write(paste(phecode, pop, ncases, ncontrols, baseline_R2liab[1], 
                paste(R2obs_values, collapse = ";"), sep = ";"), 
                file = R2liab_file, append = TRUE)
            write(paste(phecode, pop, ncases, ncontrols, baseline_R2liab[2], 
                paste(R2liab_values, collapse = ";"), sep = ";"), 
                file = R2liab_file, append = TRUE)
            write(paste(phecode, pop, ncases, ncontrols, Cstat(baseline), 
                paste(Cstat_values, collapse = ";"), sep = ";"), 
                file = cstat_file, append = TRUE)
            
            
        }
    }
}

get_R2diff <- function(dat_anc, phecode, score_num, covs){
    R2liab <- R2Liability(dat_anc, phecode, score_num, covs)
    baseline_R2liab <- R2Liability(dat_anc, phecode, NA, covs)
    return(R2liab[2] - baseline_R2liab[2])
}
extract_R2diff <- function(data, indices,  phecode, score_num, covs) {
    d <- data[indices, ]  # Resample data
    get_R2diff(d, phecode, score_num, covs)
}

perform_adj_PC_specify_covs <- function (dat, score_name, covs) 
{
    formula <- paste(score_name, "~", covs)
    mod <- lm(as.formula(formula), data = dat)
    return(residuals(mod))
}

extract_meta_R2diff <- function (data, indices, method, phecode, covs) 
{
    d <- data[indices, ]
    if(method == "PRScs"){
        avg <- mean(c(get_R2diff(d, phecode, 1, covs), 
                      get_R2diff(d, phecode, 2, covs), 
                      get_R2diff(d, phecode, 3, covs)))
        return(avg)
    }else if(method == "LDpred"){
        avg <- mean(c(get_R2diff(d, phecode, 2, covs), 
                      get_R2diff(d, phecode, 4, covs), 
                      get_R2diff(d, phecode, 6, covs)))
        return(avg)
    }
}

extract_meta_r_squared <- function (data, indices, method, phecode) 
{
    d <- data[indices, ]
    if(method == "PRScs"){
        avg <- mean(c(perform_lm(d, paste0("score", 1), phecode), 
                    perform_lm(d, paste0("score", 2), phecode), 
                    perform_lm(d, paste0("score", 3), phecode)))
        return(avg)
    }else if(method == "LDpred"){
        avg <- mean(c(perform_lm(d, paste0("score", 2), phecode), 
                    perform_lm(d, paste0("score", 4), phecode), 
                    perform_lm(d, paste0("score", 6), phecode)))
        return(avg)
    }
}
