# --- Fonction pour retrouver les valeurs d'origines
back_norm <- function(x_norm, xmin = min(data[,vi], na.rm = T), xmax = max(data[,vi], na.rm = T)){
  
  x <- x_norm * (xmax - xmin) + xmin
  return(x)
}

# CHECK ERROR FUNCTION FOR PEDIGREE ---------------------------------------
#' The function [checking_error_pedigree] checks for :
#' (i) All IDs are represetend in the dataset i.e. if missing IDs from the whole population
#' (ii) The presence of non unique id in 'individual_id'
#' (iii) Presence of an emergence date OR a capture date
#' (iv) Mismatch between emergence date and first capture date
#' (v) Social parent in individual_id modalities
#' (vi) Presence of sexe
#' (vii) Mismatch parents sexe e.g. if a mother is a female
#' 
#' @param all_id : complete list of individual ID e.g. reference list (default = NULL)
#' @param individual_id : dataset individual unique IDs vector
#' @param mother_id : mother IDs vector 
#' @param father_id : father IDs vector 
#' @param emergence_date : emergence dates vector (default = NULL)
#' @param capture_date : capture dates vector (default = NULL)
#' @param sexe : individual's sex vector (default = NULL)

# 
# all_id = capture_individual_id$individual_id
# individual_id = pedigree_actu$individual_id
# mother_id = pedigree_actu$mere_soc
# father_id = pedigree_actu$pere_soc
# emergence_date = pedigree_actu$date_emergence
# capture_date = pedigree_actu$date_capture
# sexe = pedigree_actu$sexe


checking_error_pedigree <- function(all_id = NULL,
                                   individual_id,
                                   mother_id,
                                   father_id,
                                   emergence_date = NULL,
                                   capture_date = NULL,
                                   sexe = NULL){
  
  # DATAFRAME
  data <- data.frame(
    individual_id = individual_id,
    mother_id = mother_id,
    father_id = father_id)
  
  # --- conditional columns
  if(length(emergence_date) > 0){data["emergence_date"] <- emergence_date}
  if(length(capture_date) > 0){data["capture_date"] <- capture_date}
  if(length(sexe) > 0){data["sexe"] <- sexe}
  
  
  # ERROR IDENTIFICATION
  
  #' (i) All IDs present in the dataset ?
  if(length(all_id) > 0){
    missing_id <- !(all_id %in% individual_id)  # check if all ID are presents in the pedigree
    id_missing_id <- unique(all_id[which(missing_id)])  # get missing id
  } else {id_missing_id = character(0)}
  
  
  #' (ii) Uniqueness of ID
  id_distribution <- table(individual_id) > 1
  nb_non_unique_id <- sum(id_distribution)
  non_unique_id <- names(which(id_distribution > 0))
  
  
  #' (iii) Presence of emergence date or capture date
  if(!length(emergence_date) == 0 & !length(capture_date) == 0){
    presence_emergence_or_capture <- is.na(emergence_date) & is.na(capture_date)
    pb_rows_presence_emergence_or_capture <- which(presence_emergence_or_capture)
    pb_id_presence_emergence_or_capture <- data[pb_rows_presence_emergence_or_capture, "individual_id"]  # ID of problematic data
  } else {pb_id_presence_emergence_or_capture = character(0)}
  
  
  #' (iv) emergence_date < date_1st_capture
  if(!length(emergence_date) == 0 & !length(capture_date) == 0){
    mismatch_emergence_first_capture <- emergence_date > capture_date
    pb_rows_emergence_capture <- which(mismatch_emergence_first_capture)  #  rows with mismatch
    nb_mismatch_emergence_first_capture <- sum(mismatch_emergence_first_capture, na.rm = T)  # total number of mismatch
    pb_id_row_capture <- data[pb_rows_emergence_capture, "individual_id"]  # ID of mismatch rows 
  } else {pb_id_row_capture = character(0)}

  
  #' (v) Social parents not in ID (no corresponding row)
  # mother
  non_in_id_mom <- sum(!unique(mother_id) %in% individual_id) - 1  # -1 to remove NA
  pb_rows_id_mom <- which((!mother_id %in% individual_id) & !is.na(mother_id))
  pb_id_mom <- unique(data[pb_rows_id_mom, "mother_id"])
  
  # father
  non_in_id_dad <- sum(!unique(father_id) %in% individual_id) - 1  # -1 to remove NA
  pb_rows_id_dad <- which((!father_id %in% individual_id) & !is.na(father_id))
  pb_id_dad <- unique(data[pb_rows_id_dad, "father_id"])
  
  
  if(length(sexe) > 0){
    #' (vi) Missing sexe
    na_sexe <- is.na(sexe)
    pb_id_sexe_na <- pedigree_actu[which(na_sexe),"individual_id"]
    
    
    #' (vii) mother == female and father == male
    # mother
    unique_mother_id <-  unique(na.omit(mother_id))
    mother_sexe <- data[data$individual_id %in% unique_mother_id, c("individual_id", "sexe")]
    pb_id_mismtach_mother_sexe <- mother_sexe[which(mother_sexe$sexe == "m"), "individual_id"]
    
    # father
    unique_father_id <-  unique(na.omit(father_id))
    father_sexe <- data[data$individual_id %in% unique_father_id, c("individual_id", "sexe")]
    pb_id_mismtach_father_sexe <- mother_sexe[which(father_sexe$sexe == "f"), "individual_id"]
  } else {pb_id_sexe_na = character(0) ; pb_id_mismtach_mother_sexe = character(0) ; pb_id_mismtach_father_sexe = character(0)}
  
  
  
  
  
  # OUTPUT MESSAGES 
  # if any error...
  if (nb_non_unique_id > 0 |
      length(pb_id_presence_emergence_or_capture) > 0 |
      length(pb_id_row_capture) |
      non_in_id_mom > 0 |
      non_in_id_dad > 0 |
      length(pb_id_sexe_na) > 0 |
      length(pb_id_mismtach_mother_sexe) > 0 |
      length(pb_id_mismtach_father_sexe) > 0 |
      length(id_missing_id) > 0) {
        
        
    cat("\n====================================================================== ")
    cat("\n========================== PROBLEM DETECTED ==========================")
    cat("\n======================================================================\n\n ")
    
    #(i)
    if(length(id_missing_id) > 0){
      cat("MISSING IDs IN THE DATASET : \n", id_missing_id, "\n\n")
    }
    
    # (ii) Uniqueness of ID
    if(nb_non_unique_id > 0){
      cat("NON-UNIQUE ID : \n", non_unique_id, "\n\n")
    }
    
    
    # (iii) Presence of emergence date or capture date
    if(!(length(emergence_date) == 0 & length(capture_date) == 0)){
      
      if(length(pb_id_presence_emergence_or_capture > 0)){
        cat("MISSING EMERGENCE DATE AND CAPTURE DATE, ID : \n", pb_id_presence_emergence_or_capture, "\n\n")
      }
    }
    
    
    # (iv) emergence_date < date_1st_capture
    if(!(length(emergence_date) == 0 & length(capture_date) == 0)){
      
      if(any(length(emergence_date) == 0, length(capture_date) == 0)){
        cat("MISMATCH EMERGENCE DATE and FIRST CAPTURE DATE, ID : \n")
        cat("! Both emergence date and capture date have to be added !\n\n")
      }
      
      if(length(pb_id_row_capture > 0)){
        cat("MISMATCH EMERGENCE DATE and FIRST CAPTURE DATE, ID : \n", pb_id_row_capture, "\n\n")
      }
    }
    
    
    # (v) Social parents not in ID (no corresponding row)
    if(non_in_id_mom > 0){
      cat("MOTHER MISSING IN INDIVIDUAL_ID : \n", pb_id_mom, "\n\n")
    } 
    
    if(non_in_id_dad > 0){
      cat("FATHER MISSING IN INDIVIDUAL_ID : \n", pb_id_dad, "\n\n")
    }
    
    
    # (vi) Missing sexe
    if(length(pb_id_sexe_na) > 0){
      cat("MISSING SEXE, ID : \n", pb_id_sexe_na, "\n\n")
    }
    
    
    # (vii) mother == female and father == male
    # mother
    if(length(pb_id_mismtach_mother_sexe) > 0){
      cat("MOTHER REPORTED AS MALE : \n", pb_id_mismtach_mother_sexe, "\n\n")
    }
    
    # father
    if(length(pb_id_mismtach_father_sexe) > 0){
      cat("FATHER REPORTED AS MALE : \n", pb_id_mismtach_father_sexe, "\n")
    }
  

    
    cat("\n\n[You can access problematic IDs by creating an object with the function's results using object$...]\n\n")
    
  }
  
  # No problem 
  if(!nb_non_unique_id > 0 &&
     !length(pb_id_presence_emergence_or_capture) > 0 &&
     !length(pb_id_row_capture) &&
     !non_in_id_mom > 0 &&
     !non_in_id_dad > 0 &&
     !length(pb_id_sexe_na) > 0 &&
     !length(pb_id_mismtach_mother_sexe) > 0 &&
     !length(pb_id_mismtach_father_sexe) > 0 &&
     !length(id_missing_id) > 0){
    cat("No problem detected in the pedigree.")
  }
  
  
  return(
    error = list(
      MISSING_ID = id_missing_id,
      NON_UNIQUE_ID = non_unique_id,
      MISSING_EMERGENCE_CAPTURE = pb_id_presence_emergence_or_capture,
      MISMATCH_DATE_ID = pb_id_row_capture,
      MOTHER_MISSING_ID = pb_id_mom,
      FATHER_MISSING_ID = pb_id_dad,
      MISSING_SEXE = pb_id_sexe_na,
      MISMATCH_SEX_MOTHER = pb_id_mismtach_mother_sexe,
      MISMATCH_SEX_FATHER = pb_id_mismtach_father_sexe
    )
  )
}  # ... end if error







# CHECK ERROR FUNCTION FOR PHENOTYPES & SOCIAL STRUCTURE DATA SET ---------------------------------------
#' The function [checking_error_data] checks for :
#' (i) The uniqueness of capture IDs
#' (ii) The presence of all individual IDs from the dataset in the pedigree
#' (iii) The presence of all individual IDs from the pedigree in the dataset
#' (iv) Mismatch between the year of capture and the year of territory observation
#' (v) Presence of measures at emergence while litter size == 0
#' (vi) Group size is not equal to the sum of adults + young (1 year) + yearling
#' (vii) Group size == NA while group composition is available
#' (viii) Presence of some NAs (not all) in group composition
#' 
#' @param capture_id : a vector containing all the capture IDs 
#' @param pedigree_individual_id : a vector containing all the individual IDs of the pedigree
#' @param data_individual_id : a vector containing all the individual IDs of the dataset 
#' @param capture_year : a vector containing the capture's years 
#' @param territoire_year : a vector containing the observation's territory years
#' @param litter_size_at_birth : a vector containing the litter size at birth
#' @param df_group_composition : a data_frame containing the group size in first column and one year adults and yearlings for the other

checking_error_data <- function(capture_id,
                                pedigree_individual_id,
                                data_individual_id,
                                capture_year,
                                territoire_year,
                                litter_size_at_birth,
                                df_group_composition
                                ){
  
  # ERROR IDENTIFICATION
  
  #' (i) Uniqueness of capture ID ?
  capture_id_distribution <- table(capture_id) > 1
  non_unique_capture_id <- names(which(capture_id_distribution))
  
  
  #' (ii) All individual IDs of the data present in the pedigree
  missing_id_ped <- !(data_individual_id %in% pedigree_individual_id)  # check if all ID are presents in the pedigree
  id_missing_id_ped <- unique(data_individual_id[which(missing_id_ped)])  # get missing id
  
  
  #' (iii) All individual IDs of the pedigree present the data set IDs 
  missing_id_data <- !(pedigree_individual_id %in% data_individual_id)  # check if all ID are presents in the pedigree
  id_missing_id_data <- unique(pedigree_individual_id[which(missing_id_data)])  # get missing id
  
  
  #' (iv) Year of capture == year of territory IDs
  mismatch_capture_terr_year <- capture_year != territoire_year
  mismatch_capture_terr_year_capture_id <- capture_id[which(mismatch_capture_terr_year)]
  
  
  #' (v) Presence of measures at emergence while litter size == 0
  litter_size_at_birth_mismtach <- litter_size_at_birth == 0  # which litter_size == 0
  litter_size_at_birth_capture_id <- capture_id[which(litter_size_at_birth_mismtach)]
  
  
  #' (vi) Group size is not equal to the sum of adults + young (1 year) + yearling
  group_size_mismatch <- df_group_composition[1] != apply(df_group_composition[2:length(df_group_composition)], 1, function(x) sum(x, na.rm = T))
  group_size_mismatch_capture_id <- capture_id[group_size_mismatch & !is.na(group_size_mismatch)]
  

    #' (vii) Group size == NA while group composition is available
  group_size_na_mismatch <- is.na(df_group_composition[1]) & apply(df_group_composition[2:length(df_group_composition)], 1, function(x) all(!is.na(x)))
  group_size_na_mismatch_capture_id <- capture_id[group_size_na_mismatch]
  
  
  #' (viii) Presence of some NAs (not all) in group composition
  composition_na_mismatch <- apply(df_group_composition[2:length(df_group_composition)], 1, function(x) any(is.na(x)) & !all(is.na(x)))
  composition_na_mismatch_capture_id <- capture_id[composition_na_mismatch]
  
  
  
  # OUTPUT MESSAGES 
  # if any error...
  if (length(non_unique_capture_id) > 0 |
      length(id_missing_id_ped) > 0 |
      length(id_missing_id_data) > 0 |
      length(mismatch_capture_terr_year_capture_id) > 0 |
      length(litter_size_at_birth_capture_id) > 0 |
      length(group_size_mismatch_capture_id) > 0 |
      length(group_size_na_mismatch_capture_id) > 0 |
      length(composition_na_mismatch_capture_id) > 0
      ) {
    
    
    
    cat("\n====================================================================== ")
    cat("\n========================== PROBLEM DETECTED ==========================")
    cat("\n======================================================================\n\n ")
    
    #' (i) Uniqueness of capture ID ?
    if(length(non_unique_capture_id) > 0){
      cat("NON-UNIQUE CAPTURE IDs IN THE DATASET : \n", non_unique_capture_id, "\n\n")
    }
    
    
    #' (ii) All individual IDs of the data present in the pedigree
    if(length(id_missing_id_ped) > 0){
      cat("INDIVIDUAL IDs IN DATASET BUT NOT IN PEDIGREE : \n", id_missing_id_ped, "\n\n")
    }
    
    
    #' (iii) All individual IDs of the pedigree present the data set IDs 
    if(length(id_missing_id_data) > 0){
      cat("INDIVIDUAL IDs IN PEDIGREE BUT NOT IN DATASET : \n", id_missing_id_data, "\n\n")
    }
    

    #' (iv) Year of capture == year of territory IDs
    if(length(mismatch_capture_terr_year_capture_id) > 0){
      cat("MISMTACH YEAR OF CAPTURE AND YEAR OF TERRITORY, CAPTURE ID : \n", mismatch_capture_terr_year_capture_id, "\n\n")
    }
    
    
    #' (v) Presence of measures at emergence while litter size == 0
    if(length(litter_size_at_birth_capture_id) > 0){
      cat("PRESENCE OF MEASURE AT EMERGENCE WHILE LITTER SIZE == 0, CAPTURE ID : \n", litter_size_at_birth_capture_id, "\n\n")
    }
    
    
    #' (vi) Group size is not equal to the sum of adults + young (1 year) + yearling
    if(length(group_size_mismatch_capture_id) > 0){
      cat("GROUP SIZE IS NOT EQUAL TO THE SUM OF INDIVIDUALS, CAPTURE ID : \n", group_size_mismatch_capture_id, "\n\n")
    }
    
    
    #' (vii) Group size == NA while group composition is available
    if(length(group_size_na_mismatch_capture_id) > 0){
      cat("GROUP SIZE IS NA BUT COMPOSITION AVAILABLE, CAPTURE ID : \n", group_size_na_mismatch_capture_id, "\n\n")
    }
    
    
    #' (viii) Presence of some NAs (not all) in group composition
    if(length(composition_na_mismatch_capture_id) > 0){
      cat("PRESENCE OF NAs IN GROUP COMPOSITION (!= ALL NAs), CAPTURE ID : \n", composition_na_mismatch_capture_id, "\n\n")
    }
    
    
    
    cat("\n\n[You can access problematic rows by creating an object with the function's results using object$...]\n\n")
    
  }
  
  
  
  
  
  # No problem 
  if(!length(non_unique_capture_id) > 0 &&
     !length(id_missing_id_ped) > 0 &&
     !length(id_missing_id_data) > 0 &&
     !length(mismatch_capture_terr_year_capture_id) > 0 &&
     !length(litter_size_at_birth_capture_id > 0) &&
     !length(group_size_mismatch_capture_id) > 0 &&
     !length(group_size_na_mismatch_capture_id) > 0 &&
     !length(composition_na_mismatch_capture_id) > 0
    ){
    cat("No problem detected in the dataset")
  }
  
  
  return(
    error = list(
      NON_UNIQUE_CAPTURE_ID = which(capture_id %in% non_unique_capture_id),
      MISSING_INDIVIDUAL_ID_IN_PEDIGREE = which(data_individual_id %in% id_missing_id_ped),
      MISSING_INDIVIDUAL_ID_IN_DATASET = which(pedigree_individual_id %in% id_missing_id_data),
      MISMATCH_CAPTURE_TERRITORY_YEAR = which(capture_id %in% mismatch_capture_terr_year_capture_id),
      MISMATCH_MEASURES_LITTER_SIZE_0 = which(capture_id %in% litter_size_at_birth_capture_id),
      MISMATCH_GROUP_SIZE_SUM = which(capture_id %in% group_size_mismatch_capture_id),
      GROUP_SIZE_NA_BUT_COMPOSITION_AVAILABLE = which(capture_id %in% group_size_na_mismatch_capture_id),
      COMPOSITION_NA = which(capture_id %in% composition_na_mismatch_capture_id)
      
    )
  )
}  # ... end if error






# FONCTIONS SPLINES 1 [from Christophe] -------------------------------------------------------
#' @param cov - variable explicative sur laquelle la spline cubique sera construite
#' @param nknot - nombre de noeuds de la spline cubique, déterminé automatiquement si non spécifié
#' 
spline.prep <- function(cov, nknot = NA) {
  
  if (is.na(nknot)) {
    n.knots <- max(5, min(round(length(unique(cov)) / 4), 35))
    
  } else {
    n.knots <- nknot
  }
  
  # Calcul des probabilités pour les quantiles des nœuds
  prob.tmp <- seq(0, 1, length = n.knots + 2)
  prob <- prob.tmp[-c(1, length(prob.tmp))]  # Exclusion des extrêmes (0 et 1)
  
  # Détermination des positions des nœuds basée sur les quantiles de "cov"
  knots <- quantile(unique(cov), probs = prob)
  
  # Création de la matrice de conception pour les effets fixes (X)
  X <- cbind(rep(1, length(cov)), cov, cov ^ 2) # colonne d'intercept, cov et cov^2
  
  # Création de la matrice de conception pour les effets aléatoires (Z)
  # Calcul des différences en cube entre chaque valeur de "cov" et chaque nœud
  Z.tmp <- (abs(outer(X = cov, Y = knots, "-"))) ^ 3
  
  # Calcul de la matrice de covariance entre les nœuds
  omega.all <- (abs(outer(knots, knots, "-"))) ^ 3
  
  # Décomposition en valeurs singulières (SVD) de la matrice de covariance des nœuds
  svd.omega.all <- svd(omega.all)
  
  # Transformation de la matrice des nœuds pour les effets aléatoires
  sqrt.omega.all <- t(svd.omega.all$v %*% (t(svd.omega.all$u) * sqrt(svd.omega.all$d)))
  Z <- t(solve(sqrt.omega.all, t(Z.tmp)))  # Application de la transformation
  
  # Retour des résultats sous forme de liste
  return(list(
    cov = cov,
    knots = knots,
    X = X,
    Z = Z
  ))
}



# Fonction de production du modèle matriciel spline
spl.X <- function(x, xk)
{
  q <- length(xk) + 2  # nombre de paramètres
  n <- length(x)  # nombre de données
  X <- matrix(data = 1, nrow = n, ncol = q) # Modèle matriciel inititalisé avec la première colonne contennt les intercepts
  X[,2] <- x  # la colonne deux contient les données
  X[,3:q] <- outer(X = x, Y = xk, FUN = rk)  # les autres colonnes R(x, xk)
  return(X)
}



# R(x, z) pour la fonction cubic splines
rk <- function(x, z)  
{ ((z - 0.5) ^2 - 1/12) * ((x - 0.5) ^2 -1/12)/4 - 
    ((abs(x-z) - 0.5)^4 - (abs(x-z) - 0.5)^2/2+7/240) /24
}







# FONCTIONS SPLINES 2 [From Woods, 2006] -------------------------------------------------------
#' @param covariable - A vector containing the covariable from which we want to apply spline function and calculate knots position
#' @param nknot - number of knots (~inflection points)
#' @param prediction - boolean, if output is the prediction data 
#' @param prediction_data - a vector containing the data from which we want prediction
#' 
#' @return X - model matrix with column 1 the intercept, column 2 the x values and 
#'             the other nknot columns the r(x, x_j') values with x_j' the position of knot j 
#' @return knots_position - the position of the knots calculated with the quantile of x

cubic_spline <- function(covariable = NULL, n.knots = 6, prediction = F, prediction_data = NULL, knot_position = NULL, covariable_origin = data$annee_emergence)
{
  
  # Calculate the cubic spline basis model (formula from Wood, 2006 p.126)
  rk <- function(x, z)  
  { ((z - 0.5) ^2 - 1/12) * ((x - 0.5) ^2 -1/12)/4 - 
      ((abs(x-z) - 0.5)^4 - (abs(x-z) - 0.5)^2/2+7/240) /24
  }
  
  # --- Fonction pour retrouver les valeurs d'origines
  back_norm <- function(x_norm, xmin = min(data[,vi], na.rm = T), xmax = max(data[,vi], na.rm = T)){
    
    x <- x_norm * (xmax - xmin) + xmin
    return(x)
  }
  
  
  
  if(is.null(knot_position)){
    # Knots position definition by quantile
    prob.tmp <- seq(0, 1, length = n.knots + 2)  # q (number of parameters) = n.knots + 2
    prob <- prob.tmp[-c(1, length(prob.tmp))]  # removing 0 and 1 to calculate knot probability
    knots <- quantile(unique(covariable), probs = prob)  # Knots position 
  } else {
    # Knots position defining by knot_position argument
    # --- Extract modalities to associate each normalized values with each original value
    ref <- covariable_origin[!duplicated(covariable_origin)]
    row.ref <- which(!duplicated(covariable_origin))

    # --- Extract the modality row where we want to place the knot 
    index.temp <- which(ref %in% knot_position)
    row.knots <- row.ref[index.temp]
    
    # --- Define knots position
    knots <- covariable[row.knots]
  }

  
  
  # Generation of the matricial model for spline
  if(prediction == F){
    data <- covariable
  } else if (prediction == T & !is.null(prediction_data)){
    data <- prediction_data
  } else { stop(print("Missing prediction_data while prediction == T."))}
  
  q <- n.knots + 2  # number of parameters
  n <- length(data)  # sample size
  X <- matrix(data = 1, nrow = n, ncol = q, dimnames = list(c(1:n), paste0("V", 1:q)))  # initialization of the model matrix, with first column being 1 the intercepts
  X[,2] <- data  # ... second column being the covariable OR prediction values
  X[,3:q] <- outer(X = data, Y = knots, FUN = rk)   # ... the others the R(X, Z) 
  
  
  return(list(X = X, knots_position = knots))
  
}


#' # FONCTIONS SPLINES 2 [From Woods, 2006] -------------------------------------------------------
#' #' @param covariable - A vector containing the covariable from which we want to apply spline function and calculate knots position
#' #' @param nknot - number of knots (~inflection points)
#' #' @param prediction - boolean, if output is the prediction data 
#' #' @param prediction_data - a vector containing the data from which we want prediction
#' #' 
#' #' @return X - model matrix with column 1 the intercept, column 2 the x values and 
#' #'             the other nknot columns the r(x, x_j') values with x_j' the position of knot j 
#' #' @return knots_position - the position of the knots calculated with the quantile of x
#' 
#' cubic_spline <- function(covariable = NULL, n.knots = 6, prediction = F, prediction_data = NULL)
#' {
#'   
#'   # Calculate the cubic spline basis model (formula from Wood, 2006 p.126)
#'   rk <- function(x, z)  
#'   { ((z - 0.5) ^2 - 1/12) * ((x - 0.5) ^2 -1/12)/4 - 
#'       ((abs(x-z) - 0.5)^4 - (abs(x-z) - 0.5)^2/2+7/240) /24
#'   }
#'   
#'   
#'   # Knots position definition by quantile
#'   prob.tmp <- seq(0, 1, length = n.knots + 2)  # q (number of parameters) = n.knots + 2
#'   prob <- prob.tmp[-c(1, length(prob.tmp))]  # removing 0 and 1 to calculate knot probability
#'   knots <- quantile(unique(covariable), probs = prob)  # Knots position 
#'   
#'   
#'   # Generation of the matricial model for spline
#'   if(prediction == F){
#'     data <- covariable
#'   } else if (prediction == T & !is.null(prediction_data)){
#'     data <- prediction_data
#'   } else { stop(print("Missing prediction_data while prediction == T."))}
#'   
#'   q <- n.knots + 2  # number of parameters
#'   n <- length(data)  # sample size
#'   X <- matrix(data = 1, nrow = n, ncol = q, dimnames = list(c(1:n), paste0("V", 1:q)))  # initialization of the model matrix, with first column being 1 the intercepts
#'   X[,2] <- data  # ... second column being the covariable OR prediction values
#'   X[,3:q] <- outer(X = data, Y = knots, FUN = rk)   # ... the others the R(X, Z) 
#'   
#'   
#'   return(list(X = X, knots_position = knots))
#'   
#' }


# MONOMOLECULAR MODEL -----------------------------------------------------
plotmono <- function(y0, r, maxt) {
  curve(
    1 - (1 - y0) * exp(-r * x),
    from = 0,
    to = maxt,
    xlab = 'Time',
    ylab = 'Disease Incidence',
    col = 'mediumblue'
  )
}

# plotmono(0.0017, 0.00242, 6000)




# PARTIAL RESIDUALS COMPUTATON --------------------------------------------
#' #' @param data - Data used for modelling
#' #' @param mod - lmer model result
#' #' 
#' #' @return the partial residuals of the variable "annee_emergence"
compute_partial_residual <- function(data, mod, marginal = F){
  
  # Get residual from model
  res <- residuals(mod)
  # extract estimates
  fixed_coef <- fixef(mod)  # coefficient of fixed effect
  coef_annee_emergence <- fixed_coef[grep("annee_emergence", names(fixed_coef))]  # get coefficient associated with cohort only
  
  # add random effects
  if(marginal == F){
    random_coef <- ranef(mod)
    df_random <- data.frame(mother = rownames(random_coef$mother), intercept = random_coef$mother[,1])
  }
  
  # add reference level of factor
  ref.name <- levels(drop.levels(data$annee_emergence))[1]  # get the name of the reference
  coef_annee_emergence <- c(coef_annee_emergence,  1)  # 1 is to add the reference level
  names(coef_annee_emergence)[length(coef_annee_emergence)] <- ref.name  # add name to the reference level in vector
  
  
  # Computation of partial residuals for each individual
  # --- i.e. we remove the effect of cohort on residuals (as it was not included in model)
  partial_res <- c()
  for(i in 1:nrow(data)){
    index.coef.temp <- grep(pattern = data[i, "annee_emergence"], names(coef_annee_emergence))  # get associated estimate
    
    if(marginal == T){
      partial_res[i] <- 
        res[i] + coef_annee_emergence[index.coef.temp] +  # effect of the cohort (+ because model as factor)
        fixed_coef[1] # add intercept
      
    } else if(marginal == F){
      index.coef.random.temp <- grep(pattern = data[i, "mother"], df_random$mother)  # get associated estimate with random effect level
      partial_res[i] <- 
        res[i] + coef_annee_emergence[index.coef.temp] +  # effect of the cohort
        fixed_coef[1] + # intercept
        df_random[index.coef.random.temp, "intercept"] # variance from random effect
    }
    
 
  }
  
  return(partial_res)
  
}





# FITNESS  ---------------------------------------------------------------------
# individual_id = pedigree_data$individual_id
# mother_id = pedigree_data$mother
# father_id = pedigree_data$father
# sexe = pedigree_data$sexe
# cohort = pedigree_data$cohort

# fitness <-
#   function(individual_id = pedigree_data$individual_id,
#            mother_id = pedigree_data$mother,
#            father_id = pedigree_data$father,
#            sexe = pedigree_data$sexe,
#            cohort = pedigree_data$cohort) {
# 
#   # GENERATE DF TO FACILITATE DATA MANIPULATION
#   df_temp = data.frame(individual_id = individual_id,
#                        mother_id = mother_id,
#                        father_id = father_id,
#                        sexe = sexe,
#                        cohort = cohort)
# 
# 
#   # COMPUTE LRS
#   # Function to count total number of offspring for an individual
#   LRS <- function(id) {
#     offspring_count <- sum(mother_id == id | father_id == id, na.rm = TRUE)
#     return(offspring_count)
#   }
# 
#   # Apply the function to each individual ID
#   LRS <-
#     (sapply(X = individual_id,
#             FUN =  LRS)  # for the argument of the function count_offspring
#     )
# 
#   # Add to df
#   df_temp["LRS"] <- LRS
# 
# 
# 
#   # MEAN LRS PER SEXE PER COHORT
#   df_temp %>%
#     group_by(sexe, cohort) %>%
#     summarise(n = sum(LRS > 0), w = mean(LRS, na.rm = T)) %>% 
#     print(n=100)
# 
#   }
# 
# 


# COMPUTE GENETIC PARAMETERS FOR UNIVARIATE MODEL -------------------------------------------------------
#' @param mod_result - MCMCglmm model result
#' @param re - names of random effect to display
#' 
#' @return Median and credibility intervalle of the specified random effects re
#' 
compute_stats_re <- function(mod_result = mod_result, re) {
  for (i in re) {
    params <- mod_result$VCV[, i] / rowSums(mod_result$VCV)
    cat(paste0(
      i, " : ", round(median(params), 3), " [",
      paste0(round(c(HPDinterval(params)[1], HPDinterval(params)[2]), 3), collapse = "; "),
      "]", "\n"
    ))
  }
}




# PLOT POSTERIOR -------------------------------------------------------
#' @param mod_result - MCMCglmm model result
#' @param re - names of random effect to display
#' 
#' @return Posterir distribution of specified random effects re
#' 
plot_posterior_re <- function(mod_result = mod_result, re) {
  for (i in re) {
    par(mfrow = c(1, 3))
    pos.temp <- mod_result$VCV[, i]
    plot(mcmc(pos.temp), main = i)
  }
}



# PLOT LME4 RESIDUALS -------------------------------------------------------
#' @param mod_result_lme4 - lme4 model result
#' @param y - vector of observed values for response variable
#' 
#' @return Residuals evaluation of lme4 model
#' 
plot_lme4_residuals <- function(mod_result_lme4, y) {
    par(mfrow = c(1, 3))
    
    # PREDICTION VERSUS OBSERVED
    plot(fitted(mod_result_lme4) ~ y, xlab = "Observed", ylab = "Fitted")
    abline(c(0, 1), col = "red", lwd = 2)
    
    # RESIDUALS
    plot(residuals(mod_result_lme4), ylab = "Residuals")
    abline(h = 0, lwd = 2, col = "red")
    
    # HITOGRAM RESIDUALS
    hist(residuals(mod_result_lme4), breaks = 100, xlab = "Residuals", main = NULL)
    par(mfrow = c(1, 1))
}

# PLOT MCMCglmm RESIDUALS -------------------------------------------------------
#' @param mod_result - lme4 model result
#' @param data - data frame used for model (for prediction)
#' @param y - vector of observed values for response variable  

#' @return Residuals evaluation of MCMCglmm model
#' 
plot_mcmc_residuals <- function(mod_result, data, y) {
  par(mfrow = c(1, 3))
  
  # RESIDUALS COMPUTATIO
  pred <- MCMCglmm::predict.MCMCglmm(object = mod_result)
  res <- y - pred
  
  # PLOT PREDICTION ~ DATA
  plot(pred ~ y, xlab = "Observed", ylab = "Prediction")
  abline(c(0, 1), lwd = 2, col = "red")
  
  # PLOT RESIDUALS
  plot(res, xlab = "Rank", ylab = "Residuals")
  abline(h = 0, lwd = 2, col = "red")
  
  # PLOT RESIDUALS
  hist(res, breaks = 50, xlab = "Residuals", main = NULL)
  
  par(mfrow = c(1, 1))
}





# GAM ---------------------------------------------------------------------
#' @param mod_result - lme4 model result
#' @param data - data frame used for model 
#' @param periode - periode on which to compute phenotype variation 
#' @param vi - predictor we wish to 'GAM' (should be "annee_emergence" or "annee_capture") 
#' @param trait - trait to be predicted on vi
#' @param xlab - label of the x axis
#' @param ylab - label of the y axis
#' @param title - title of the plot
#' 
#' @return Return a plot of the temporal phenotype variation of the trait

plot_gam <- function(data = data_mcmc_masse_marmotton,
                     n.knots = 1, 
                     periode = NULL,
                     knot_position = 2005,
                     covariable_origin = data_mcmc_masse_marmotton$annee_emergence,
                     vi = "annee_emergence",
                     trait = "masse",
                     xlab = "Cohort",
                     ylab = "Masse",
                     title = "Pups",
                     ylim = NULL){
  
  # SELECT PERIOD OF YEAR
  if(!is.null(periode)){
    data <- data[  data[,vi] %in% periode,]
  }
  
  # data[,vi] <- drop.levels(  data[,vi])  # drop unused levels
  
  # CONVERT VI TO NUMERIC
  data[,vi] <- as.numeric(as.character(data[,vi]))
  
  
  # NORMALIZATION
  data[,"vi_gam"]  = data[,vi] - min(data[,vi])
  data[,"vi_gam"] <- data[,"vi_gam"] / max(data[,"vi_gam"])
  
  
  
  # GENERATION OF THE CUBIC MODEL
  X <- cubic_spline(covariable = data[,"vi_gam"],
                    n.knots = n.knots,
                    knot_position = knot_position,
                    covariable_origin =  data[,vi])$X
  
  # DATASET FOR MODEL
  data_mod <- data.frame(trait = data[,trait], X)  
  
  # Model creation to get the estimates associated with the b_j(x) basis functions
  mod.1 <- lm(trait ~ . - 1, data = data_mod)  # - 1 to remove intercept
  coef(mod.1)
  
  x_prediction <- 0:100/100  # values to predict
  
  # prediction
  new <- cubic_spline(covariable = data[,"vi_gam"],  # needed to calculate knots position 
                      n.knots = n.knots,
                      prediction = T,
                      prediction_data = x_prediction,
                      knot_position = knot_position,
                      covariable_origin =  data[,vi])$X
  
  new <- as.data.frame(new)  # format for predict() function
  names(new) <-  colnames(model.matrix(mod.1))  # Matching new_data names with model's predictors names
  
  
  # CONFIDENCE INTERVAL
  pred_ic <- predict(object = mod.1,
                     se.fit = T,  # get standard error associated with fit
                     newdata = new,
                     interval = c("confidence"),
                     level = 0.975  # levle of confidence
  )
  
  # PREDICTION INTERVAL
  pred_prediction <- predict(object = mod.1,
                             se.fit = T,  # get standard error associated with fit
                             newdata = new,
                             interval = c("prediction"),
                             level = 0.975  # levle of confidence
  )
  
  # DATAFRAME WITH SPLINED PREDICTION AND INTERVALS
  data_predict <-
    data.frame(
      x_prediction = x_prediction,  # data used for prediction
      ic_neg = pred_ic$fit[, "lwr"],  # lower IC
      ic_pos = pred_ic$fit[, "upr"],  # upper IC
      fitted = pred_ic$fit[, "fit"],  # fitted values
      se = pred_ic$se.fit,  # stadard error associated with fitted values
      pred_neg = pred_prediction$fit[, "lwr"],
      pred_pos = pred_prediction$fit[, "upr"]
    )
  
  
  # POSITION OF KNOTS 
  # X POSITION 
  X_knot_location <- cubic_spline(covariable = data[,"vi_gam"],  # needed to calculate knots position 
                                  n.knots = n.knots,
                                  prediction = T,
                                  prediction_data = x_prediction,
                                  knot_position = knot_position,
                                  covariable_origin =  data[,vi])$knots
  
  
  X_knots <- cubic_spline(  # matrix of prediction
    covariable = data[,"vi_gam"],
    n.knots = n.knots,
    prediction = T,
    prediction_data = c(min(data[,"vi_gam"]), X_knot_location, max(data[,"vi_gam"])),
    knot_position = knot_position,
    covariable_origin =  data[,vi]
  )
  
  Y_knots <- as.vector(t(X_knots$X %*% coef(mod.1)))  # y position of knots
  
  
  # CLASSIC REGRESSION
  mod_lm <- lm(trait ~ V2, data = data_mod)
  
  # --- prediction
  pred_ic_lm <- predict(object = mod_lm,
                             se.fit = F,  # get standard error associated with fit
                             newdata = data.frame(V2 = x_prediction),  # get linear only
                             interval = c("prediction"),
                             level = 0.975  # levle of confidence
  )
  
  # --- store IC interval and prediction
  data_predict_lm <-
    data.frame(
      x_prediction = x_prediction,  # data used for prediction
      ic_neg = pred_ic_lm[, "lwr"],  # lower IC
      ic_pos = pred_ic_lm[, "upr"],  # upper IC
      fitted = pred_ic_lm[, "fit"]  # fitted values
    )
  
  # PLOT
  # NAMES OF LABELS
  lab_lm <- TeX("$\\hat{y} = \\beta_0 + \\beta_1 x$")
  lab_gam <- TeX("$\\hat{y} = \\sum^q_{i=1} f(x_i)$")
  
  
  
  # NAMES FOR THE X AXIS (FROM 0 / 1 TO REAL VALUES OF AGE)
  # --- Fonction pour retrouver les valeurs d'origines
  back_norm <- function(x_norm, xmin = min(data[,vi], na.rm = T), xmax = max(data[,vi], na.rm = T)){
    
    x <- x_norm * (xmax - xmin) + xmin
    return(x)
  }
  
  
  # COMPUTE MEAN and IC FOR EACH YEAR
  compute_ci <- function(sd){
    n = length(sd)
    return(qnorm(p = .95) * sd / sqrt(n))
  }
  
  data_vi <- 
    data %>% 
    group_by(!!sym(vi)) %>% 
    summarise(m = mean(!!sym(trait)), x = mean(!!sym(vi)),
              ymin = mean(!!sym(trait)) - compute_ci(sd = sd(!!sym(trait))),
              ymax = mean(!!sym(trait)) + compute_ci(sd = sd(!!sym(trait))))
  
  
  data_vi["period"] <- ifelse(test = data_vi$annee_emergence <=2005, "periode1", "period2")
  
  
  
  
  
  
  
  # PLOT VISUALISATION
plot <- ggplot() +
    # Data weight ~ age
    geom_point(
      data = data,
      mapping = aes(x = !!sym(vi), y = !!sym(trait)),
      alpha = 0.1,
      size = .8
    ) +
    
    # MEAN PER AGE YEAR --------------
  geom_point(
    data = data_vi,
    mapping = aes(x = !!sym(vi), y = m, col = period),
    alpha = 1,
    # color = "#171717",
    size = 2
  ) +
  
  scale_color_manual(values = c("#4e9785", "#cc3300"))  + #  c("#4e9785", "#cc3300", "#36465d")
  
  geom_errorbar(data = data_vi, mapping = aes(x = !!sym(vi), ymin = ymin, ymax = ymax, col = period),
                width = 0.5,
                # col = "#171717",
                size = 0.5
                ) +
    
    xlim(c(min(periode), max(periode))) +
    
    # geom_errorbar(data = mean_trait_year,
    #               mapping = aes(x = age_year, ymin = m - se, ymax = m + se),
    #               width = .2) +
    

   # # # IC interval for regression line
   # geom_ribbon(
   #   data = data_predict_lm,
   #   mapping = aes(x = back_norm(x_prediction), ymin = ic_neg, ymax = ic_pos),
   #   alpha = 0.1,
   #   col = "darkred",
   #   # linetype = "longdash",
   #   linewidth = NA,
   #   fill = "darkred"
   # ) +
    
    # geom_label(aes(x = 8, y = 6.5), label = lab_lm, fill = "darkred", parse = T,
    #            hjust = 0, vjust = 0, size = 6, color = "white", family = "mono") +
    

  # # X LABELS OF DATE
  # # Set custom x-axis tick labels
  scale_x_continuous(
    breaks = seq(1997, 2023, 4),
    labels = seq(1997, 2023, 4),
    expand = c(0.05, 0.05),
    limits = c(1997, 2023)
  ) +
  
  
  # geom_label(aes(x = 12.5, y = 3.8), label = lab_gam, fill = "black", parse = T,
  #            hjust = 0, vjust = 0, size = 6, color = "white", family = "mono") +
  

    
  tidyquant::theme_tq() +
  
  
    theme(
      # panel.grid = element_line(color = "#EEEEE0", linetype = "dashed"),  # grid of plot
      # background color
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      # axis.line = element_line(color = "black"),
      # box around plot
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      
      
      # Increase space between axis and labels
      axis.title.y = element_text(size = 18, margin = margin(0, 9, 0, 0), face = "bold"),  # space between ylab and plot
      axis.title.x = element_text(face = "bold"),
      axis.text.x = element_text(margin = margin(t = 6)),
      axis.text.y = element_text(margin = margin(r = 6)),
      plot.title = element_text(size = 18, margin = margin(b = 6)),
      # Text size
      text = element_text(size = 18)
    ) +
  
    # Label xlab, ylab and title
    labs(x = xlab,
         y = ylab,
         title = title) 
  # theme(aspect.ratio=1)  # shape plot as squared
  
  



# ADD GAM modeling to plot if "tibia"
if(mean(data_vi$m, na.rm = T) < 100){  # bidouillage
  size_multiple <- 0.8  # to adjust size of points
  
  plot <- plot +     # GAM -------------------------
  geom_line(
    data = data_predict,
    mapping = aes(x = back_norm(x_prediction), y = fitted), #104E8B
    linewidth = 1, col = "#545454" 
  ) +
    
    # # # IC interval for GAM
    # geom_ribbon(
    #   data = data_predict,
    #   mapping = aes(x = back_norm(x_prediction), ymin = ic_neg, ymax = ic_pos),
    #   alpha = 0.3,
    #   col = "#292929",
    #   # linetype = "longdash",
    #   linewidth = NA,
    #   fill = "#292929"
    # ) +
    
    # KNOTS --------------------------
  geom_point(
    aes(x = c(min(data[,vi], na.rm = T), back_norm(X_knots$knots_position), max(data[,vi], na.rm = T)), y = Y_knots),
    size = c(4*size_multiple, rep(4.5*size_multiple, length(X_knots$knots_position)), 2*size_multiple),
    col = "#545454",
    shape = 15
  ) +
    
    geom_point(
      aes(x = c(min(data[,vi], na.rm = T), back_norm(X_knots$knots_position), max(data[,vi], na.rm = T)), y = Y_knots),
      size = c(2.5*size_multiple, rep(2.5*size_multiple, length(X_knots$knots_position)), 4*size_multiple),
      col = c("#545454", rep("white", length(X_knots$knots_position)), "#545454"),
      shape = 15
    ) +
  
  # ADD 95% PREDICTION INTERVAL
  geom_line(data = data_predict,
            mapping = aes(x = back_norm(x_prediction), y = pred_neg),
            linetype = 3,
            col = "#545454",
            linewidth = 1) +
    
    geom_line(data = data_predict,
              mapping = aes(x = back_norm(x_prediction), y = pred_pos),
              linetype = 3,
              col = "#545454",
              linewidth = 1)
  
  
    
} else {
  
  # AJOUT LINEAIRE SI MASSE
  plot <- plot +
    # LM CLASSIC -----------------
  geom_smooth(
    data = data,
    mapping = aes(x =as.numeric(!!sym(vi)), y = !!sym(trait)),
    method = "lm", 
    linewidth = 1,
    col = "#545454",
    # fill = "#292929",
    alpha = 0.5,
    se = F
  ) +
    
    
    # ADD 95% PREDICTION INTERVAL
    geom_line(data = data_predict_lm,
              mapping = aes(x = back_norm(x_prediction), y = ic_neg),
              linetype = 3,
              col = "#545454",
              linewidth = 1) +
    
    geom_line(data = data_predict_lm,
              mapping = aes(x = back_norm(x_prediction), y = ic_pos),
              linetype = 3,
              col = "#545454",
              linewidth = 1)
    
}

  if(!is.null(ylim)){
    plot <- plot +
      ylim(c(ylim[1], ylim[2]))
  }
  
  return(plot)
}




# CHOICE KNOT NUMBER FOR SIZE AND MASSE -----------------------------------
#' @param data - data frame used to test number of knot needed
#' @param knot_tested - maximum number of knot tested
#' 
#' @return Return 2x2 plot for tibia / masse with the r squared of model as a function of knot number

compute_r_squared_knots <- 
  function(data, knot_tested = 6){
  # Storage of R²
  r_mass <- c()
  r_tibia <- c()
  r_mass_full <- c()
  r_tibia_full <- c()
  
  # ADD R FOR CLASSIC REGRESSION WITH NO GAM
  fixed_effects <- c(  # fixed effect
    "sexe",
    "capture_time",
    "litter_size_naissance",
    "age_d",
    "as.numeric(annee_emergence)",
    "taille_groupe_naissance")
  
  # MASSE
  vd <- "masse"
  formula_fixed <- as.formula(paste0(vd,  # response variable
                                     "~",
                                     paste0(fixed_effects, collapse = " + "),
                                     "+ (1|mother) + (1|annee_emergence)"
  ))
  
  mod_mass <- lme4::lmer(formula_fixed,
                         data = data[!is.na(data[, vd]),])
  
  r.temp <- r.squaredGLMM(mod_mass)
  r_mass <- c(r_mass, r.temp[1])
  r_mass_full <- c(r_mass_full, r.temp[2])
  
  
  # TIBIA
  fixed_effects <- c(  # fixed effect
    "sexe",
    "capture_time",
    "litter_size_naissance",
    "age_d",
    "as.numeric(annee_emergence)",
    "taille_groupe_naissance")
  
  
  # FORMULA FOR FIXED EFFECT IN THE MODEL
  vd <- "tibia"
  formula_fixed <- as.formula(paste0(vd,  # response variable
                                     "~",
                                     paste0(fixed_effects, collapse = " + "),
                                     "+ (1|mother) + (1|annee_emergence)"
  ))
  
  mod_tibia <- lme4::lmer(formula_fixed,
                          data = data[!is.na(data[, vd]),])
  
  r.temp <- r.squaredGLMM(mod_tibia)
  r_tibia <- c(r_tibia, r.temp[1])
  r_tibia_full <- c(r_tibia_full, r.temp[2])
  
  
  # NUMBER OF KNOTS TESTED
  knot_tested <- knot_tested
  
  
  
  for(i in 1:knot_tested){
    # CUBIC SPLINES
    n.knots <- i
    X <- cubic_spline(covariable = as.numeric(data[,"annee_emergence_gam"]), n.knots = n.knots)$X
    
    data <- cbind(data, X)
    
    # GET MODEL WITH CUBIC SPLINE
    bx <- paste0("V", 2:(n.knots+2))
    
    fixed_effects <- c(  # fixed effect
      "sexe",
      "capture_time",
      "litter_size_naissance",
      "age_d",
      "taille_groupe_naissance")
    
    
    # MASS --------------------------------------------------
    # FORMULA FOR FIXED EFFECT IN THE MODEL
    vd <- "masse"
    formula_fixed <- as.formula(paste0(vd,  # response variable
                                       "~",
                                       paste0(bx, collapse = " + "),  # cubic age
                                       " + ",
                                       paste0(fixed_effects, collapse = " + "),
                                       "+ (1|mother) + (1|annee_emergence)"
    ))
    
    mod_mass <- lme4::lmer(formula_fixed,
                           data = data[!is.na(data[, vd]),]
    )
    
    
    r.temp <- r.squaredGLMM(mod_mass)
    r_mass <- c(r_mass, r.temp[1])
    r_mass_full <- c(r_mass_full, r.temp[2])
    
    
    
    # TIBIA --------------------------------------------------
    # GET MODEL WITH CUBIC SPLINE
    bx <- paste0("V", 2:(n.knots+2))
    
    fixed_effects <- c(  # fixed effect
      "sexe",
      "capture_time",
      "litter_size_naissance",
      "age_d",
      "taille_groupe_naissance")
    
    
    # FORMULA FOR FIXED EFFECT IN THE MODEL
    vd <- "tibia"
    formula_fixed <- as.formula(paste0(vd,  # response variable
                                       "~",
                                       paste0(bx, collapse = " + "),  # cubic age
                                       " + ",
                                       paste0(fixed_effects, collapse = " + "),
                                       "+ (1|mother) + (1|annee_emergence)"
    ))
    
    mod_tibia <- lme4::lmer(formula_fixed,
                            data = data[!is.na(data[, vd]),]
    )
    
    
    r.temp <- r.squaredGLMM(mod_tibia)
    r_tibia <- c(r_tibia, r.temp[1])
    r_tibia_full <- c(r_tibia_full, r.temp[2])
    
  }
  
  return(list(r_mass = r_mass, r_tibia = r_tibia, r_mass_full = r_mass_full, r_tibia_full = r_tibia_full))
}












#' # GAM ---------------------------------------------------------------------
#' #' @param mod_result - lme4 model result
#' #' @param data - data frame used for model 
#' #' @param periode - periode on which to compute phenotype variation 
#' #' @param vi - predictor we wish to 'GAM' (should be "annee_emergence" or "annee_capture") 
#' #' @param trait - trait to be predicted on vi
#' #' @param xlab - label of the x axis
#' #' @param ylab - label of the y axis
#' #' @param title - title of the plot
#' #' 
#' #' @return Return a plot of the temporal phenotype variation of the trait
#' 
#' plot_gam <- function(data = data_mcmc_masse_marmotton,
#'                      n.knots = 2, 
#'                      periode = NULL,
#'                      vi = "annee_emergence",
#'                      trait = "tibia",
#'                      xlab = "Cohort",
#'                      ylab = "Tibia",
#'                      title = "Pups",
#'                      ylim = NULL){
#'   
#'   # SELECT PERIOD OF YEAR
#'   if(!is.null(periode)){
#'     data <- data[data$annee_emergence %in% periode,]
#'   }
#'   
#'   data$annee_emergence <- drop.levels(data$annee_emergence)  # drop unused levels
#'   
#'   # CONVERT VI TO NUMERIC
#'   data[,vi] <- as.numeric(as.character(data[,vi]))
#'   
#'   
#'   # NORMALIZATION
#'   data[,"vi_gam"]  = data[,vi] - min(data[,vi])
#'   data[,"vi_gam"] <- data[,"vi_gam"] / max(data[,"vi_gam"])
#'   
#'   
#'   
#'   # GENERATION OF THE CUBIC MODEL
#'   X <- cubic_spline(covariable = data[,"vi_gam"], n.knots = n.knots)$X
#'   
#'   # DATASET FOR MODEL
#'   data_mod <- data.frame(trait = data[,trait], X)  
#'   
#'   # Model creation to get the estimates associated with the b_j(x) basis functions
#'   mod.1 <- lm(trait ~ . - 1, data = data_mod)  # - 1 to remove intercept
#'   coef(mod.1)
#'   
#'   x_prediction <- 0:100/100  # values to predict
#'   
#'   # prediction
#'   new <- cubic_spline(covariable = data[,"vi_gam"],  # needed to calculate knots position 
#'                       n.knots = n.knots,
#'                       prediction = T,
#'                       prediction_data = x_prediction)$X
#'   
#'   new <- as.data.frame(new)  # format for predict() function
#'   names(new) <-  colnames(model.matrix(mod.1))  # Matching new_data names with model's predictors names
#'   
#'   
#'   # CONFIDENCE INTERVAL
#'   pred_ic <- predict(object = mod.1,
#'                      se.fit = T,  # get standard error associated with fit
#'                      newdata = new,
#'                      interval = c("confidence"),
#'                      level = 0.95  # levle of confidence
#'   )
#'   
#'   # PREDICTION INTERVAL
#'   pred_prediction <- predict(object = mod.1,
#'                              se.fit = T,  # get standard error associated with fit
#'                              newdata = new,
#'                              interval = c("prediction"),
#'                              level = 0.95  # levle of confidence
#'   )
#'   
#'   # DATAFRAME WITH PREDICTION AND INTERVALS
#'   data_predict <-
#'     data.frame(
#'       x_prediction = x_prediction,  # data used for prediction
#'       ic_neg = pred_ic$fit[, "lwr"],  # lower IC
#'       ic_pos = pred_ic$fit[, "upr"],  # upper IC
#'       fitted = pred_ic$fit[, "fit"],  # fitted values
#'       se = pred_ic$se.fit,  # stadard error associated with fitted values
#'       pred_neg = pred_prediction$fit[, "lwr"],
#'       pred_pos = pred_prediction$fit[, "upr"]
#'     )
#'   
#'   
#'   # POSITION OF KNOTS 
#'   # X POSITION 
#'   X_knot_location <- cubic_spline(covariable = data[,"vi_gam"],  # needed to calculate knots position 
#'                                   n.knots = n.knots,
#'                                   prediction = T,
#'                                   prediction_data = x_prediction)$knots
#'   
#'   
#'   X_knots <- cubic_spline(  # matrix of prediction
#'     covariable = data[,"vi_gam"],
#'     n.knots = n.knots,
#'     prediction = T,
#'     prediction_data = c(min(data[,"vi_gam"]), X_knot_location, max(data[,"vi_gam"]))
#'   )
#'   
#'   Y_knots <- as.vector(t(X_knots$X %*% coef(mod.1)))  # y position of knots
#'   
#'   
#'   
#'   
#'   sum(!data_mcmc_masse_marmotton$individual_id %in% data_mcmc_tibia_marmotton$individual_id)
#'   
#'   # PLOT
#'   # NAMES OF LABELS
#'   lab_lm <- TeX("$\\hat{y} = \\beta_0 + \\beta_1 x$")
#'   lab_gam <- TeX("$\\hat{y} = \\sum^q_{i=1} f(x_i)$")
#'   
#'   
#'   
#'   # NAMES FOR THE X AXIS (FROM 0 / 1 TO REAL VALUES OF AGE)
#'   # --- Fonction pour retrouver les valeurs d'origines
#'   back_norm <- function(x_norm, xmin = min(data[,vi], na.rm = T), xmax = max(data[,vi], na.rm = T)){
#'     
#'     x <- x_norm * (xmax - xmin) + xmin
#'     return(x)
#'   }
#'   
#'   
#'   # COMPUTE MEAN FOR EACH YEAR
#'   data_vi <- 
#'     data %>% 
#'     group_by(!!sym(vi)) %>% 
#'     summarise(m = mean(!!sym(trait)), x = mean(!!sym(vi)),  sd = sd(!!sym(trait), na.rm = T))
#'   
#'   
#'   # --- Plot model
#'   plot <- ggplot() +
#'       # Data weight ~ age
#'       geom_point(
#'         data = data,
#'         mapping = aes(x = !!sym(vi), y = !!sym(trait)),
#'         alpha = 0.1,
#'         size = .7
#'       ) +
#'       
#'       # # Set custom x-axis tick labels
#'       # scale_x_continuous(
#'       #   breaks = seq(0, 16, 2),
#'       #   labels = seq(0, 16, 2),
#'       #   expand = c(0.05, 0.05)
#'       # ) +
#'       
#'       # MEAN PER AGE YEAR --------------
#'     geom_point(
#'       data = data_vi,
#'       mapping = aes(x = !!sym(vi), y = m),
#'       alpha = 1,
#'       color = "#4A4A4A",
#'       size = 2
#'     ) +
#'       
#'       # geom_errorbar(data = mean_trait_year,
#'       #               mapping = aes(x = age_year, ymin = m - se, ymax = m + se),
#'       #               width = .2) +
#'       
#'       # LM CLASSIC -----------------
#'     geom_smooth(
#'       data = data,
#'       mapping = aes(x =as.numeric(!!sym(vi)), y = !!sym(trait)), method = "lm", se = F,
#'       linewidth = 1.2, col = "darkred"
#'     ) +
#'       
#'       # geom_label(aes(x = 8, y = 6.5), label = lab_lm, fill = "darkred", parse = T,
#'       #            hjust = 0, vjust = 0, size = 6, color = "white", family = "mono") +
#'       
#'       # GAM -------------------------
#'     geom_line(
#'       data = data_predict,
#'       mapping = aes(x = back_norm(x_prediction), y = fitted),
#'       linewidth = 1.2, col = "black"
#'     ) +
#'       
#'       # # IC interval
#'       # geom_ribbon(
#'       #   data = data_predict,
#'       #   mapping = aes(x = seq(0, 16, length.out = 101), ymin = pred_neg, ymax = pred_pos),
#'       #   alpha = 0.3,
#'       #   col = "black",
#'       #   linetype = "dashed",
#'       #   linewidth = 1,
#'       #   fill = NA
#'       # ) +
#'     
#'     # geom_label(aes(x = 12.5, y = 3.8), label = lab_gam, fill = "black", parse = T,
#'     #            hjust = 0, vjust = 0, size = 6, color = "white", family = "mono") +
#'     
#'     # KNOTS --------------------------
#'     geom_point(
#'       aes(x = c(min(data[,vi], na.rm = T), back_norm(X_knots$knots_position), max(data[,vi], na.rm = T)), y = Y_knots),
#'       size = c(4, rep(4.5, length(X_knots$knots_position)), 2), col = "black"
#'     ) +
#'       
#'       geom_point(
#'         aes(x = c(min(data[,vi], na.rm = T), back_norm(X_knots$knots_position), max(data[,vi], na.rm = T)), y = Y_knots),
#'         size = c(2.5, rep(2.5, length(X_knots$knots_position)), 4), col = c("black", rep("white", length(X_knots$knots_position)), "black")
#'       ) +
#'       
#'       theme(
#'         panel.grid = element_line(color = "#EEEEE0", linetype = "dashed"),  # grid of plot
#'         # background color
#'         panel.background = element_rect(fill = "white"),
#'         # axis.line = element_line(color = "black"),
#'         # box around plot
#'         panel.border = element_rect(
#'           color = "black",
#'           fill = NA,
#'           linewidth = 1
#'         ),
#'         
#'         # Increase space between axis and labels
#'         axis.text.x = element_text(margin = margin(t = 6)),
#'         axis.text.y = element_text(margin = margin(r = 6)),
#'         
#'         # Text size
#'         text = element_text(size = 20, family = "sans")
#'       ) +
#'     
#'     # Label xlab, ylab and title
#'       labs(x = xlab,
#'            y = ylab,
#'            title = title) 
#'        # theme(aspect.ratio=1)  # shape plot as squared
#'   
#' 
#'   if(!is.null(ylim)){
#'     plot <- plot +
#'       ylim(c(ylim[1], ylim[2]))
#'   }
#'   
#'   return(plot)
#' }
#' 







# COMPUTE GENETIC PARAMETERS FOR MULTIVARIATE MODELS  -------------------------------------------------------
#' @param mod_result - MCMCglmm model result
#' @param trait - trait from which calculate parameters
#' @param decimal - precision wanted for genetic parameters report
#' @param non_us_re - # random variable(s) with no covariance us() 
#' @param position - "median" or "mean", used to compute genetic parameters

#'@WARNINGS : trait have to be entered in the same way as in the model 

#' @return a dataframe containing genetic parameters statistics : median and IC 

compute_stats_re_bivariate <- function(mod_result = mod_result,
                                       decimal = 2,
                                       trait = c("masse", "tibia"),
                                       position = "mean")
{
  
  # GET RE NAMES
  re <- colnames(mod_result$VCV)
  
  # IS FITNESS IN MODEL ?
  fitness_bool <- any(grepl(pattern = "LRS", x = re))
  
  # EXTRACT RELEVANT VARIABLE FOR DENOMINATOR TRAIT IN GENETIC PARAMETERS COMPUTATION
  for(trait.temp in trait){
    trait_re.temp <- re[grep(pattern = paste0(paste0("trait", trait.temp, ":trait", trait.temp), collapse = "|"), x = re)]
    assign(paste0(trait.temp, "_re"), trait_re.temp)
  }
  
  # Get cohort variance effect [syntax depend on the number of trait and the presence of fitness as cohort effect is not dependent on cohort]
  if(fitness_bool == T){
    # If more than two non-fitness trait included in model, set covariance name with cohort (trait[i], trait[j])
    if(length(trait)>1){
      for(i in  1:length(trait)){
        add.cohort.re.temp <- c(get(paste0(trait[i], "_re")),
                                paste0("at.level(trait, c(\"",trait[1], "\", \"", trait[2], "\"))", which(trait == trait[i]), ":at.level(trait, c(\"", trait[1], "\", \"", trait[2], "\"))", which(trait == trait[i]),".annee_emergence"))
        
        assign(paste0(trait[i], "_re"), add.cohort.re.temp)
      }
      
      #if only one non-fitness trait, set covariance name with cohort as (trait[i], trait[i])
    } else if(length(trait) == 1){
      add.cohort.re.temp <- c(get(paste0(trait, "_re")),
                              paste0("at.level(trait, c(\"", trait, "\")):at.level(trait, c(\"", trait, "\")).annee_emergence"))
      assign(paste0(trait, "_re"), add.cohort.re.temp)
    }
    
  } 
  
  
  
  # DATAFRAME TO STORE RESULTS
  df <- as.data.frame(matrix(ncol = length(trait) * 3,  # for y / ymin / ymax
                             nrow = length(get(paste0(trait[1], "_re"))) + 1)
  )
  
  # EXTRACT DYNAMICALLY NAME FROM COLNAMES $VCV
  # --- Get names of RE
  # if(!fitness_bool == T){
  
  variance_parameters <- gsub(pattern = paste0("trait", trait[1], ":trait", trait[1], ".|trait", trait[1], "."),
                              replacement = "",
                              x = get(paste0(trait[1], "_re")))
  
  # --- Rename cohort effect as "cohort"
  cohort_name <- variance_parameters[grepl(pattern = "emergence", x = variance_parameters)]
  variance_parameters[which(variance_parameters == cohort_name)] <- "cohort"
  
  # --- Add as rownames
  rownames(df) <- c(variance_parameters, "vf")  # change colnames for parameters values summary
  
  
  #} 
  
  
  
  # GET Vf THE VARIANCE ASSOCIATED WITH FIXED EFFECT (De Villemereuil, 2018)
  compute_vcvpred <- function(beta, design_matrix, ntraits) { 
    
    list(cov(matrix(design_matrix %*% beta, ncol = ntraits))) 
    
  }
  
  X <- mod_result$X  # desgin matrix
  X_mat <- as.matrix(X)  # convert to matrix
  col.names.temp <- colnames(mod_result$Sol)  # colnames of prediction
  # get index associating with last fixed effect
  last_row <- grep(pattern = "animal", col.names.temp)[1] - 1
  last_fixed_effect_index <- ifelse(!is.na(last_row), last_row, length(col.names.temp)) # needs to adjust because fixed effects change according to period modeled
  sol <- mod_result$Sol[,1:last_fixed_effect_index]  # select only fixed effects
  
  # ntrait ~ fitness included or not
  if(fitness_bool == T){ntraits = length(trait)+1}else{ntraits = length(trait)}
  fixed_effects_var <- flatten( apply(sol, 1, compute_vcvpred, design_matrix = X, ntraits = ntraits) )  # get variance-covariance matrix associated with prediction
  
  # extract variance associated with each trait
  for(trait.temp in trait){
    vf.temp <- unlist(lapply(X = fixed_effects_var, FUN = function(x) x[which(trait.temp == trait), which(trait.temp == trait)]))
    assign(paste0("vf_", trait.temp), vf.temp)
  }  
  
  
  
  # --- Change colnames containing median [or mean], lower IC and upper IC
  names.temp <- c()
  for(trait.temp in trait){
    for(stats in c("_y", "_ymin", "_ymax")){
      names.temp <- c(names.temp, paste0(trait.temp, stats))
    }
  }  
  
  colnames(df) <- names.temp
  
  
  
  # COMPUTATION OF GENETIC PARAMETERS
  # --- Function used to compute genetic parameters
  if(position == "mean"){
    position_f <- function(x){mean(x, na.rm = T)}
  } else if (position == "median"){
    position_f <- function(x){median(x, na.rm = T)}
  }
  
  
  # Result for genetic parameter (variance)
  df_variance <- df
  
  # Result for genetic parameter (ratio)
  df_ratio <- df
  
  # COMPUTE MEAN AND IC PER TRAIT FOR EACH RANDOM EFFECT
  # --- Variance
  col.index <- 1  # to ajdust column filling according to the number of trait
  for(trait.temp in trait){  # for each trait
    for(i in 1:length(get(paste0(trait[1], "_re")))){  # for each RE
      df_variance[i, col.index] <- round(position_f(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]]), decimal) 
      df_variance[i, col.index + 1] <- round(HPDinterval(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]]), decimal)[1] 
      df_variance[i, col.index + 2] <- round(HPDinterval(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]]), decimal)[2] 
    }
    col.index <- col.index + 3  # get to next trait
  }
  
  # add Vf
  col.index <- 1  # to ajdust column filling according to the number of trait
  for(trait.temp in trait){
    df_variance[nrow(df_variance), col.index] <- round(position_f(as.mcmc(get(paste0("vf_", trait.temp)))), decimal) 
    df_variance[nrow(df_variance), col.index + 1] <- round(HPDinterval(as.mcmc(get(paste0("vf_", trait.temp)))), decimal)[1] 
    df_variance[nrow(df_variance), col.index + 2] <- round(HPDinterval(as.mcmc(get(paste0("vf_", trait.temp)))), decimal)[2] 
    
    col.index <- col.index + 3  # get to next trait
  }

  
  
  
  # --- Ratio
  # COMPUTE MEAN AND IC FOR RELATIVE VARIANCE PARAMETERS
  col.index <- 1  # to ajdust column filling according to the number of trait
  for(trait.temp in trait){  # for each trait
    for(i in 1:length( get(paste0(trait[1], "_re")))){  # for each RE
      df_ratio[i, col.index] <- round(position_f(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]] / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal)
      df_ratio[i, col.index + 1] <- round(HPDinterval(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]] / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal)[1] 
      df_ratio[i, col.index + 2] <- round(HPDinterval(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]] / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal)[2] 
    }
    col.index <- col.index + 3  # get to next trait
  }
  
  # add Vf
  col.index <- 1  # to ajdust column filling according to the number of trait
  for(trait.temp in trait){
    df_ratio[nrow(df_ratio), col.index] <- round(position_f(as.mcmc(get(paste0("vf_", trait.temp))) / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal) 
    df_ratio[nrow(df_ratio), col.index + 1] <- round(HPDinterval(as.mcmc(get(paste0("vf_", trait.temp))) / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal)[1] 
    df_ratio[nrow(df_ratio), col.index + 2] <- round(HPDinterval(as.mcmc(get(paste0("vf_", trait.temp))) / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal)[2] 
    
    col.index <- col.index + 3  # get to next trait
  }
  
  
  
  # REFORMAT TABLE FOR GGPLOT
  reformat_gg <- function(data){
    data["params"] <- rownames(data)
    
    # --- Get all y from each trait
    df <- data.frame(y = eval(parse(text = paste0("c(", paste0(paste0("data$", trait[1:length(trait)], "_y"), collapse = ", "), ")"))))
    
    # --- Add RE as rownames and trait names
    df["params"] <- rep(rownames(data), length(trait))
    df["trait"] <- rep(x = trait, each = length(get(paste0(trait[1], "_re"))))
    
    # --- Get IC for each trait
    df["ymin"] <-  eval(parse(text = paste0("c(", paste0(paste0("data$", trait[1:length(trait)], "_ymin"), collapse = ","), ")")))
    df["ymax"] <-  eval(parse(text = paste0("c(", paste0(paste0("data$", trait[1:length(trait)], "_ymax"), collapse = ","), ")")))
    
    return(df)
  }
  
  
  
  # RETURN VALUES OF PARAMETERS AND PLOT
  return(list(variance = t(as.data.frame(df_variance)),
              ratio = t(as.data.frame(df_ratio))
  ))
  
}



# COMPUTE GENETIC PARAMETERS FOR TRIVARIATE MODELS  -------------------------------------------------------
#' @param mod_result - MCMCglmm model result
#' @param trait - trait from which calculate parameters
#' @param decimal - precision wanted for genetic parameters report
#' @param non_us_re - # random variable(s) with no covariance us() 
#' @param position - "median" or "mean", used to compute genetic parameters

#'@WARNINGS : trait have to be entered in the same way as in the model 

#' @return a dataframe containing genetic parameters statistics : median and IC 

compute_stats_re_trivariate <- function(mod_result = mod_result,
                                       decimal = 2,
                                       trait = c("masse", "tibia"),
                                       position = "mean")
{
  
  # GET RE NAMES
  re <- colnames(mod_result$VCV)
  
  # IS FITNESS IN MODEL ?
  fitness_bool <- any(grepl(pattern = "LRS", x = re))
  
  # EXTRACT RELEVANT VARIABLE FOR DENOMINATOR TRAIT IN GENETIC PARAMETERS COMPUTATION
  for(trait.temp in trait){
    trait_re.temp <- re[grep(pattern = paste0(paste0("trait", trait.temp, ":trait", trait.temp), collapse = "|"), x = re)]
    assign(paste0(trait.temp, "_re"), trait_re.temp)
  }
  
  # Get cohort variance effect [syntax depend on the number of trait and the presence of fitness as cohort effect is not dependent on cohort]
  if(fitness_bool == T){
    # If more than two non-fitness trait included in model, set covariance name with cohort (trait[i], trait[j])
    if(length(trait)>1){
      for(i in  1:length(trait)){
        add.cohort.re.temp <- c(get(paste0(trait[i], "_re")),
                                paste0("at.level(trait, c(\"",trait[1], "\", \"", trait[2], "\"))", which(trait == trait[i]), ":at.level(trait, c(\"", trait[1], "\", \"", trait[2], "\"))", which(trait == trait[i]),".annee_emergence"))
        
        assign(paste0(trait[i], "_re"), add.cohort.re.temp)
      }
      
      #if only one non-fitness trait, set covariance name with cohort as (trait[i], trait[i])
    } else if(length(trait) == 1){
      add.cohort.re.temp <- c(get(paste0(trait, "_re")),
                              paste0("at.level(trait, c(\"", trait, "\")):at.level(trait, c(\"", trait, "\")).annee_emergence"))
      assign(paste0(trait, "_re"), add.cohort.re.temp)
    }
    
  } 
  
  LRS_relative_re <- LRS_relative_re[-length(LRS_relative_re)]
  
  # DATAFRAME TO STORE RESULTS
  df <- as.data.frame(matrix(ncol = length(trait) * 3,  # for y / ymin / ymax
                             nrow = length(get(paste0(trait[1], "_re"))) + 1)
  )
  
  # EXTRACT DYNAMICALLY NAME FROM COLNAMES $VCV
  # --- Get names of RE
  # if(!fitness_bool == T){
  
  variance_parameters <- gsub(pattern = paste0("trait", trait[1], ":trait", trait[1], ".|trait", trait[1], "."),
                              replacement = "",
                              x = get(paste0(trait[1], "_re")))
  
  # --- Rename cohort effect as "cohort"
  cohort_name <- variance_parameters[grepl(pattern = "emergence", x = variance_parameters)]
  variance_parameters[which(variance_parameters == cohort_name)] <- "cohort"
  
  # --- Add as rownames
  rownames(df) <- c(variance_parameters, "vf")  # change colnames for parameters values summary
  
  
  #} 
  
  
  
  # GET Vf THE VARIANCE ASSOCIATED WITH FIXED EFFECT (De Villemereuil, 2018)
  compute_vcvpred <- function(beta, design_matrix, ntraits) { 
    
    list(cov(matrix(design_matrix %*% beta, ncol = ntraits))) 
    
  }
  
  X <- mod_result$X  # desgin matrix
  X_mat <- as.matrix(X)  # convert to matrix
  col.names.temp <- colnames(mod_result$Sol)  # colnames of prediction
  # get index associating with last fixed effect
  last_row <- grep(pattern = "animal", col.names.temp)[1] - 1
  last_fixed_effect_index <- ifelse(!is.na(last_row), last_row, length(col.names.temp)) # needs to adjust because fixed effects change according to period modeled
  sol <- mod_result$Sol[,1:last_fixed_effect_index]  # select only fixed effects
  
  # ntrait ~ fitness included or not
  if(fitness_bool == T){ntraits = length(trait)}else{ntraits = length(trait)}
  fixed_effects_var <- flatten( apply(sol, 1, compute_vcvpred, design_matrix = X, ntraits = ntraits) )  # get variance-covariance matrix associated with prediction
  
  # extract variance associated with each trait
  for(trait.temp in trait){
    vf.temp <- unlist(lapply(X = fixed_effects_var, FUN = function(x) x[which(trait.temp == trait), which(trait.temp == trait)]))
    assign(paste0("vf_", trait.temp), vf.temp)
  }  
  
  
  
  # --- Change colnames containing median [or mean], lower IC and upper IC
  names.temp <- c()
  for(trait.temp in trait){
    for(stats in c("_y", "_ymin", "_ymax")){
      names.temp <- c(names.temp, paste0(trait.temp, stats))
    }
  }  
  
  colnames(df) <- names.temp
  
  
  
  # COMPUTATION OF GENETIC PARAMETERS
  # --- Function used to compute genetic parameters
  if(position == "mean"){
    position_f <- function(x){mean(x, na.rm = T)}
  } else if (position == "median"){
    position_f <- function(x){median(x, na.rm = T)}
  }
  
  
  # Result for genetic parameter (variance)
  df_variance <- df
  
  # Result for genetic parameter (ratio)
  df_ratio <- df
  
  # COMPUTE MEAN AND IC PER TRAIT FOR EACH RANDOM EFFECT
  # --- Variance
  col.index <- 1  # to ajdust column filling according to the number of trait
  for(trait.temp in trait){  # for each trait
    for(i in 1:length(get(paste0(trait[1], "_re")))){  # for each RE
      df_variance[i, col.index] <- round(position_f(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]]), decimal) 
      df_variance[i, col.index + 1] <- round(HPDinterval(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]]), decimal)[1] 
      df_variance[i, col.index + 2] <- round(HPDinterval(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]]), decimal)[2] 
    }
    col.index <- col.index + 3  # get to next trait
  }
  
  # add Vf
  col.index <- 1  # to ajdust column filling according to the number of trait
  for(trait.temp in trait){
    df_variance[nrow(df_variance), col.index] <- round(position_f(as.mcmc(get(paste0("vf_", trait.temp)))), decimal) 
    df_variance[nrow(df_variance), col.index + 1] <- round(HPDinterval(as.mcmc(get(paste0("vf_", trait.temp)))), decimal)[1] 
    df_variance[nrow(df_variance), col.index + 2] <- round(HPDinterval(as.mcmc(get(paste0("vf_", trait.temp)))), decimal)[2] 
    
    col.index <- col.index + 3  # get to next trait
  }
  
  
  
  
  # --- Ratio
  # COMPUTE MEAN AND IC FOR RELATIVE VARIANCE PARAMETERS
  col.index <- 1  # to ajdust column filling according to the number of trait
  for(trait.temp in trait){  # for each trait
    for(i in 1:length( get(paste0(trait[1], "_re")))){  # for each RE
      df_ratio[i, col.index] <- round(position_f(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]] / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal)
      df_ratio[i, col.index + 1] <- round(HPDinterval(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]] / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal)[1] 
      df_ratio[i, col.index + 2] <- round(HPDinterval(mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]] / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal)[2] 
    }
    col.index <- col.index + 3  # get to next trait
  }
  
  # add Vf
  col.index <- 1  # to ajdust column filling according to the number of trait
  for(trait.temp in trait){
    df_ratio[nrow(df_ratio), col.index] <- round(position_f(as.mcmc(get(paste0("vf_", trait.temp))) / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal) 
    df_ratio[nrow(df_ratio), col.index + 1] <- round(HPDinterval(as.mcmc(get(paste0("vf_", trait.temp))) / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal)[1] 
    df_ratio[nrow(df_ratio), col.index + 2] <- round(HPDinterval(as.mcmc(get(paste0("vf_", trait.temp))) / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))), decimal)[2] 
    
    col.index <- col.index + 3  # get to next trait
  }
  
  
  
  # REFORMAT TABLE FOR GGPLOT
  reformat_gg <- function(data){
    data["params"] <- rownames(data)
    
    # --- Get all y from each trait
    df <- data.frame(y = eval(parse(text = paste0("c(", paste0(paste0("data$", trait[1:length(trait)], "_y"), collapse = ", "), ")"))))
    
    # --- Add RE as rownames and trait names
    df["params"] <- rep(rownames(data), length(trait))
    df["trait"] <- rep(x = trait, each = length(get(paste0(trait[1], "_re"))))
    
    # --- Get IC for each trait
    df["ymin"] <-  eval(parse(text = paste0("c(", paste0(paste0("data$", trait[1:length(trait)], "_ymin"), collapse = ","), ")")))
    df["ymax"] <-  eval(parse(text = paste0("c(", paste0(paste0("data$", trait[1:length(trait)], "_ymax"), collapse = ","), ")")))
    
    return(df)
  }
  
  
  
  # RETURN VALUES OF PARAMETERS AND PLOT
  return(list(variance = t(as.data.frame(df_variance)),
              ratio = t(as.data.frame(df_ratio))
  ))
  
}


# PLOT BREEDING VALUES  -------------------------------------------------------
#' @param mod_result - lme4 model result
#' @param data - data frame used for model 
#' @param trait - trait to compute
#' @param position - "median" or "mean", used to calculate bv 
#' @param periode - a list of equal length than "trait" containg the periode used to compute tmeporal variation of BV
#' @param method_regression - "lm" or "loess", regression method to compute temporal trend 
#' @param main - list containing title of upper plots

#' @return Return mean breeding values according to cohort 

plot_breeding_values <- function(mod_result,
                                           data,
                                           trait = c("tibia", "masse"),
                                           periode = list(1997:2023, 1997:2023),
                                           position = "mean",
                                           method_regression = "lm",
                                           main = list(c("Pups : Tibia"), c("Pups : Mass")),
                                           color_period = c("#4e9785", "#cc3300")
)
{
  
  # COMPUTE BV --------------------------------------------------------------
  # EXTRACTION BV
  bv <- mod_result$Sol
  
  for(trait.temp in trait){
    
    # Extract pattern to remove to extract individual id from bv names from $Sol
    if(length(trait) == 2){
      pattern <- paste0("trait", trait.temp, ".animal.")
    } else if (length(trait) == 1){
      pattern <- paste0("animal.")
    }
    
    # EXTRACTION ROWS CORRESPONDING TO BV (based on names 'animal')
    row_bv_temp <- grep(pattern = pattern, x = dimnames(bv)[[2]])
    assign(x = paste0("row_bv_", trait.temp), row_bv_temp)
    
    # CALCULATION OF BV PER INDIVIDUAL
    # --- Select position indexx computation according to function argument "position"
    if(position == "mean"){
      position_f <- function(x){mean(x, na.rm = T)}
    } else if (position == "median"){
      position_f <- function(x){median(x, na.rm = T)}
    }
    
    # --- EBV per individual
    individual_bv.temp <- apply(bv[, get(paste0("row_bv_", trait.temp))], 2, function(x) position_f(x)) 
    assign(paste0("summary_bv_", trait.temp), individual_bv.temp)  # Dynamically assign object
    
    
    # DATAFRAME CONTAINING BV AND INDIVIDUAL IDs
    df.temp <- data.frame(bv = get(paste0("summary_bv_", trait.temp)),
                          # extract name from cell names
                          individual_id = gsub(pattern = pattern,
                                               replacement = "",
                                               x = names(get(paste0("summary_bv_", trait.temp))))
    )
    
    assign(paste0("df_", trait.temp), df.temp)
    
    
    # ADD YEAR OF OBSERVATION
    merge.temp <- merge(x = get(paste0("df_", trait.temp)),
                        y = data[c("individual_id", "annee_emergence")], by = "individual_id", all.x = T)
    
    assign(paste0("df_summary_bv_", trait.temp), merge.temp)
    
    # SELECT PERIOD OVER WHICH COMPUTE BREEDING VALUES
    df.filter.temp <- 
      get(paste0("df_summary_bv_", trait.temp)) %>% 
      filter(annee_emergence %in% periode[[which(trait == trait.temp)]])  # filter per periode for trait i
    
    df.filter.temp$annee_emergence <- drop.levels(df.filter.temp$annee_emergence)  # drop unuused levels
    df.filter.temp["periode"] <- ifelse(as.numeric(as.character(df.filter.temp$annee_emergence)) <= 2005, "period1", "period2")
    
    assign(paste0("df_summary_bv_", trait.temp), df.filter.temp)
    
  }

    # PLOT --------------------------------------------------------------------
  #' [HISTOGRAMS]
  # HISTOGRAM OF BV
  plot_gg_bv_histo <- function(data){
    
    
    # Center hist on 0
    if(abs(min(data$bv, na.rm = T)) > abs(max(data$bv, na.rm = T))){
      lim <- abs(min(data$bv, na.rm = T))
    } else {
      lim <- abs(max(data$bv, na.rm = T))
    }
    
    # HISTOGRAM
    ggplot(data = data, mapping = aes(x = bv, fill = periode)) +
      geom_histogram(color = "black",
                     # fill = "white",
                     linewidth = 1, alpha = 0.4,
                     bins = 30) +
      theme_classic() +
      scale_fill_manual(values = c(color_period[1], color_period[2])) +
      labs(x = "Predicted breeding values", y = "Frequency") +
      xlim(c(-lim, lim)) +
      geom_vline(xintercept = 0, linewidth = 1.1, linetype = "longdash", col = "black") +
      tidyquant::theme_tq() +
      
      theme(
        # Increase space between axis and labels
        axis.title.y = element_text(size = 18, margin = margin(0, 8, 0, 0)),  # space between ylab and plot
        axis.text.x = element_text(margin = margin(t = 6)),
        axis.text.y = element_text(margin = margin(r = 6)),
        # Text size
        text = element_text(size = 18),
        legend.position = "none")
  }
  
  # --- Histogram per trait
  for(trait.temp in trait){
    plot.hist.temp <- plot_gg_bv_histo(get(paste0("df_summary_bv_", trait.temp)))
    assign(paste0("plot_hist_", trait.temp), plot.hist.temp)
  }

    #' [TEMPORAL TREND PER COHORT]
  # STATS PER COHORT
  compute_stat_per_cohort <- function(data, position = position){
    
    # If we want the mean of BV
    # --- Function to compute IC
    if(position == "mean"){
      
      compute_ci <- function(sd){
        n = length(sd)
        return(qnorm(p = .95) * sd / sqrt(n))
      }
      
      output <- 
        data %>% 
        group_by(annee_emergence) %>% 
        summarise(y = mean(bv),
                  ymin = mean(bv) - compute_ci(sd = sd(bv)),
                  ymax = mean(bv) + compute_ci(sd = sd(bv))) 
      
      # If we want the median
    } else if(position == "median"){
      
      output <- 
        data %>% 
        group_by(annee_emergence) %>% 
        summarise(y = median(bv),
                  ymin = quantile(bv, probs = 0.1),
                  ymax = quantile(bv, probs = 0.9)) 
    }
    
    return(output)
  }
  
  
  # --- Compute stats per trait
  for(trait.temp in trait){
    gg.temp <- compute_stat_per_cohort(data = get(paste0("df_summary_bv_", trait.temp)), position = position)
    gg.temp["periode"] <- ifelse(as.numeric(as.character(gg.temp$annee_emergence)) <= 2005, "period1", "period2")
    assign(paste0("gg_", trait.temp), gg.temp)
  }
  
  

    # PLOT
  plot_gg_bv <- function(data, main, method_regression = method_regression, ylab = ylab){
    # # Reformating cohort modality names
    # labels_cohort <- substr(x = data$annee_emergence, start = 3, stop = 4)
    # labels_cohort <- factor(labels_cohort, levels = labels_cohort)
    # data$annee_emergence <- labels_cohort
    
    ggplot(data = na.omit(data), mapping = aes(x = as.numeric(as.character(annee_emergence)), y = y, col = periode)) +
      geom_point(size = 2) +
      geom_smooth(data = na.omit(data),
                  mapping = aes(x = as.numeric(as.character(annee_emergence)), y = y, col = periode),
                  method = method_regression,
                  se = T) +
      geom_errorbar(data = data, aes(ymin = ymin, ymax = ymax), width= 0.5) +
      theme_classic() +
      scale_color_manual(values = c(color_period[1], color_period[2])) +
      labs(x = "Cohort",
           y = ylab,
           title = main) +
      
      tidyquant::theme_tq() +
      
 theme(
         # Increase space between axis and labels
         axis.title.y = element_text(size = 18, margin = margin(0, 8, 0, 0)),  # space between ylab and plot
         axis.text.x = element_text(margin = margin(t = 6)),
         axis.text.y = element_text(margin = margin(r = 6)),
         # Text size
         text = element_text(size = 18),
         legend.position = "none") +
     
       # # X LABELS OF DATE
      # # Set custom x-axis tick labels
      scale_x_continuous(
        breaks = seq(1997, 2023, 4),
        labels = seq(1997, 2023, 4),
        expand = c(0.05, 0.05),
        limits = c(1997, 2023)
      ) 
  }
  
  # --- Plot temporal pattern pattern of breeding values
  for(trait.temp in trait){
    if(trait.temp == "tibia"){
      ylab <- "PBV tibia length (mm)"
    } else if (trait.temp == "masse"){
      ylab <- "PBV mass (g)"
    }
    plot.cohort.temp <- plot_gg_bv(data = get(paste0("gg_", trait.temp)),
                                   main = main[[which(trait == trait.temp)]],
                                   method_regression = method_regression,
                                   ylab = ylab)  
    assign(paste0("plot_cohort_", trait.temp), plot.cohort.temp)
    
  }
  
  
  
  # RETURN PLOTS
  list_hist <- list()
  list_cohort <- list()
  i = 1  # index for list 
  
  for(trait.temp in trait){
    # --- Histograms
    list_hist[[i]] <- get(paste0("plot_hist_", trait.temp))
    
    # --- Cohort plots
    list_cohort[[i]] <- get(paste0("plot_cohort_", trait.temp))
    # update index
    i = i + 1
  }
  
  # --- Display all four plots
  grid.arrange(grobs = list_cohort, ncol = length(trait), nrow = 1)
  
  
}








#  DRIFT FOR BIVARIATE MODELS -------------------------------------------------------
#' @param mod_result - MCMCglmm model result
#' @param pedigree - pedigree used for model 
#' @param position - "median" or "mean", used to compute breeding values
#' @param data - data used for modelling to remove NAs for trait BV calculation 
#' @param periode1 - first periode to test 
#' @param periode2 - "median" or "mean", used to compute breeding values 
#' @param GenerationTime - generation time to express result per generation
#' @param trait - trait tested for BV ["tibia" and/or "masse"] 
#' @param iteration - number of MCMC iteration 

#' @return Return mean breeding values according to cohort 
 
# mod_result = mod_result
# pedigree = pedigree_data_biv_marmotton
# position = "mean"
# periode_pooled_tibia = 1997:2023
# periode_pooled_masse = 1992:2023
# periode1_tibia = 1990:2005
# periode2_tibia = 2006:2023
# periode1_masse = 1990:2005
# periode2_masse = 2006:2023
# trait = c("tibia", "masse")
# iteration = 100


test_drift <- function(mod_result = mod_result,
                       pedigree = pedigree_data_biv_marmotton,
                       # data = data_mcmc_biv_marmotton,
                       position = "mean",
                       periode_pooled_tibia = 1997:2023,
                       periode_pooled_masse = 1992:2023,
                       periode1_tibia = 1990:2005,
                       periode2_tibia = 2006:2023,
                       periode1_masse = 1990:2005,
                       periode2_masse = 2006:2023,
                       GenerationTime = 5.9,
                       trait = c("tibia", "masse"),
                       iteration = 100,
                       color_period = c("#4e9785", "#cc3300", "#36465d")
                       ){
  
  # Saving results
  post_df <- data.frame(trait = NA, ebv = NA, ebv_drift = NA)
  post_df_periode1 <- data.frame(trait = NA, ebv = NA, ebv_drift = NA)
  post_df_periode2 <- data.frame(trait = NA, ebv = NA, ebv_drift = NA)
  
  row.temp = 0  # inititalize for loop for each trait
  
  
  for(trait.temp in trait){
    
    # FORMATTING AND DATA EXTRACTION ------------------------------------------
    # COLUMN CONTAINING EBV FOR SIZE
    col.name.temp <- colnames(mod_result$Sol)
    col.size.indexes <- grep(pattern = paste0(trait.temp, ".animal"), x = col.name.temp)
    
    # --- Extract $Sol corresponding to size only (for computation time)
    sol_size <- mod_result$Sol[, col.size.indexes]
    
    # --- Renaming col with individual id
    colnames(sol_size) <- gsub(pattern = paste0("trait", trait.temp,".animal.", collapse = ""), replacement = "", x = colnames(sol_size))
    
    
    # CONSTRUCT DATAFRAME WITH INDIVIDUAL ID / COHORT / EBV
    # --- Extract individual ID associated with EBV
    id_ebv <- gsub(pattern = paste0("trait", trait.temp,".animal.", collapse = ""), replacement = "", x = colnames(mod_result$Sol)[col.size.indexes])  # remove extra text to get individual id
    df_bv <- data.frame(individual_id = id_ebv)  # formatting for merging
    
    # --- Extract cohort date associated with id
    df_bv <- merge(x = pedigree_data[c("individual_id", "cohort", "inbreeding")], y = df_bv, by = "individual_id")
    
    # --- Reorder df_bv to match the order of individual id of $Sol
    matching.indexes <- match(x = colnames(sol_size), table = df_bv$individual_id)
    df_bv <- df_bv[matching.indexes,]
    
    
    # REMOVE NA FOR TRAIT
    # --- in df_bv
    # id.na <- data_mcmc_biv_marmotton$individual_id[which(is.na(data_mcmc_biv_marmotton[trait.temp]))]
    # row.id.na <- which(df_bv$individual_id %in% id.na)
    # df_bv <- df_bv[-row.id.na,]
    # 
    # # --- in MCMC $Sol
    # col.id.na <- which(colnames(sol_size) %in% id.na)
    # sol_size <- sol_size[,-col.id.na]
    
    
    # MEAN EBV PER COHORTE AND DRIFT SIMULATION -------------------------------
    # FOR LOOP TO CALCULATE SIGNIFICANCE OF EBV TEMPORAL PATTERN AGAINST DRIFT
    if(position == "mean"){
      position_f <- function(x){mean(x, na.rm = T)}
    } else if (position == "median"){
      position_f <- function(x){median(x, na.rm = T)}
    }
    
    for(i in 1:iteration) { # for each iteration
      
      # Mean ebv per cohort per iteration
      YEAR.MEAN.MCMC <- tapply(X = sol_size[i,],
                               INDEX =  df_bv$cohort,
                               FUN = function(x) position_f(x))
      
      # Simulate breeding values with VA = estimate at iteration i
      repbv <- rbv(pedigree = pedigree[1:3],
                   G = mod_result$VCV[,paste0("trait", trait.temp,":trait", trait.temp, ".animal", collapse = "")][i]  # additive covariance matrix
      )
      
      # --- remove NAs measure for repbv
      # repbv <- repbv[-row.id.na]
      
      # Mean ebv per cohort per iteration of simulated data
      # --- Filter id %in% pedigree
      
      YEAR.MEAN.MCMC.REPLICATE <- tapply(X = repbv,
                                         INDEX = df_bv$cohort[df_bv$individual_id %in% pedigree$individual_id],
                                         FUN = function(x) position_f(x))
      
      
      # --- Global
      # select period over which to compute slope of bv for pooled period
      if(trait.temp == "tibia"){year_periode_pooled <- periode_pooled_tibia}else{year_periode_pooled <- periode_pooled_masse}  # adjust depending to trait test [define in function parameters]
      year_periode_pooled_indexes <- which(names(YEAR.MEAN.MCMC) %in% year_periode_pooled)
      
      # ... for periode 1 
      if(trait.temp == "tibia"){year_periode1 <- periode1_tibia}else{year_periode1 <- periode1_masse}  # adjust depending to trait test [define in function parameters
      year_periode1_indexes <- which(names(YEAR.MEAN.MCMC) %in% year_periode1)
      
      # ... for periode 2
      if(trait.temp == "tibia"){year_periode2 <- periode2_tibia}else{year_periode2 <- periode2_masse}  # adjust depending to trait test [define in function parameters
      year_periode2_indexes <- which(names(YEAR.MEAN.MCMC) %in% year_periode2)
      
      
      # REGRESSION SLOPES -------------------------------------------------------
      # # df with threshold manipulation
      # 
      # if(threshold_model == TRUE){
      #   df.temp <- data.frame(y = YEAR.MEAN.MCMC[year_periode_pooled_indexes],
      #                         ydrift = YEAR.MEAN.MCMC.REPLICATE[year_periode_pooled_indexes],
      #                         x_pooled = as.numeric(as.character(names(YEAR.MEAN.MCMC[year_periode_pooled_indexes]))))
      #   
      #   df.temp$x_periode1 <- df.temp$x_periode2 <- df.temp$x_pooled
      #   df.temp$x_periode1[df.temp$x_periode1 > 2005] <- 2005
      #   df.temp$x_periode2[df.temp$x_periode2 <= 2005] <- 2005
      #   df.temp$x_periode2[df.temp$x_periode2 > 2005] <- df.temp$x_periode2[df.temp$x_periode2 > 2005]
      #   
      #   
      #   # extract slope of bv ~ year
      #   # POOLED
      #   post_df[i + row.temp, "ebv"] <- summary(lm(y ~ x_pooled, data = df.temp))$coeff[2]  
      #   post_df[i + row.temp, "ebv_drift"] <- summary(lm(ydrift ~ x_pooled, data = df.temp))$coeff[2]
      #   
      #   # PERIODES
      #   mod <- summary(lm(y ~ x_periode1 + x_periode2, data = df.temp))  # Threshold model
      #   mod.drift <- summary(lm(ydrift ~ x_periode1 + x_periode2, data = df.temp))  # Threshold model
      #   
      #   # --- Periode 1
      #   post_df_periode1[i + row.temp, "ebv"] <- mod$coeff[2]  # periode 1 [we take the first coefficient]
      #   post_df_periode1[i + row.temp, "ebv_drift"] <- mod.drift$coeff[2]  # periode 1 [we take the first coefficient]
      #   
      #   # --- Periode 2
      #   post_df_periode2[i + row.temp, "ebv"] <- mod$coeff[3]  # periode 2 [we take the second coefficient]
      #   post_df_periode2[i + row.temp, "ebv_drift"] <- mod.drift$coeff[3]  # periode 2 [we take the second coefficient]
      #   
      #   
      # } else {
        
      
      
      # LINEAR REGRESSION if EBV calculated with mean
        if(position == "mean"){
          # extract slope of bv ~ year
          post_df[i + row.temp, "ebv"] <- summary(lm(YEAR.MEAN.MCMC[year_periode_pooled_indexes] ~
                                                       as.numeric(as.character(names(YEAR.MEAN.MCMC[year_periode_pooled_indexes])))))$coef[2]*GenerationTime
          
          post_df[i + row.temp, "ebv_drift"] <- summary(lm(YEAR.MEAN.MCMC.REPLICATE[year_periode_pooled_indexes] ~
                                                             as.numeric(as.character(names(YEAR.MEAN.MCMC.REPLICATE[year_periode_pooled_indexes])))))$coef[2]*GenerationTime
          
          # --- Periode 1 : phneotype decreasing
          # extract slope of bv ~ year
          post_df_periode1[i + row.temp, "ebv"] <- summary(lm(YEAR.MEAN.MCMC[year_periode1_indexes] ~
                                                                as.numeric(as.character(names(YEAR.MEAN.MCMC[year_periode1_indexes])))))$coef[2]*GenerationTime
          
          post_df_periode1[i + row.temp, "ebv_drift"] <- summary(lm(YEAR.MEAN.MCMC.REPLICATE[year_periode1_indexes] ~
                                                                      as.numeric(as.character(names(YEAR.MEAN.MCMC.REPLICATE[year_periode1_indexes])))))$coef[2]*GenerationTime
          
          # Periode 2 : phenotype increasing
          
          # extract slope of bv ~ year
          post_df_periode2[i + row.temp, "ebv"] <- summary(lm(YEAR.MEAN.MCMC[year_periode2_indexes] ~
                                                                as.numeric(as.character(names(YEAR.MEAN.MCMC[year_periode2_indexes])))))$coef[2]*GenerationTime
          
          post_df_periode2[i + row.temp, "ebv_drift"] <- summary(lm(YEAR.MEAN.MCMC.REPLICATE[year_periode2_indexes] ~
                                                                      as.numeric(as.character(names(YEAR.MEAN.MCMC.REPLICATE[year_periode2_indexes])))))$coef[2]*GenerationTime
          
          
          
          # QUANTILE REGRESSION if position == "median"
        } else if(position == "median"){
          # extract slope of bv ~ year
          # --- quantreg::rq function need a dataframe to work
          df.temp <- data.frame(y = YEAR.MEAN.MCMC,
                                x = as.numeric(names(YEAR.MEAN.MCMC)),
                                y_drift = YEAR.MEAN.MCMC.REPLICATE)
          
          # --- Compute slope 
          post_df[i + row.temp, "ebv"] <- summary(rq(y ~ x, data = df.temp[year_periode_pooled_indexes,]))$coef[2]  # observed
          post_df[i + row.temp, "ebv_drift"] <- summary(rq(y_drift ~ x, data = df.temp[year_periode_pooled_indexes,]))$coef[2]  # drift
          
          # --- Periode 1 : phneotype decreasing
          # extract slope of bv ~ year
          post_df_periode1[i + row.temp, "ebv"] <- summary(rq(y ~ x, data = df.temp[year_periode1_indexes,]))$coef[2]  # observed
          post_df_periode1[i + row.temp, "ebv_drift"] <- summary(rq(y_drift ~ x, data = df.temp[year_periode1_indexes,]))$coef[2]  # drift
          
          
          # Periode 2 : phenotype increasing
          # extract slope of bv ~ year
          post_df_periode2[i + row.temp, "ebv"] <- summary(rq(y ~ x, data = df.temp[year_periode2_indexes,]))$coef[2]  # observed
          post_df_periode2[i + row.temp, "ebv_drift"] <- summary(rq(y_drift ~ x, data = df.temp[year_periode2_indexes,]))$coef[2]  # drift
        }
      
      
      
      
    }
    
    # Add trait corresponding to result
    post_df[(1 + row.temp):(iteration + row.temp), "trait"] <- rep(trait.temp, iteration)
    post_df_periode1[(1 + row.temp):(iteration + row.temp), "trait"] <- rep(trait.temp, iteration)
    post_df_periode2[(1 + row.temp):(iteration + row.temp), "trait"] <- rep(trait.temp, iteration)
    
    # Row index for next trait
    row.temp = iteration
  }
  
  
  # COMPUTE PSEUDO P VALUES PER PERIODE
  # --- function to calculate number of time ebv > drift
  compute_p <- function(data, test = c("drift", "sign_evol"), trait, iteration = iteration){
    df.temp <- data[data$trait == trait,]
    
    if(test == "drift"){
      return(sum(df.temp$ebv < df.temp$ebv_drift) / iteration)
    } else if(test == "sign_evol"){
      return(sum(df.temp$ebv < 0) / iteration)
    }
  }
  
  # --- Compute results per periode dataset per trait
  dataset <- c("post_df", "post_df_periode1", "post_df_periode2")
  trait.test <- trait[1:length(trait)]
  
  # --- DRIFT (PROBA EBV > DRIFT)
  # matrix for result
  p_result_drift <- as.data.frame(matrix(NA, ncol = 2, nrow = 3,
                                         dimnames = list(dataset[1:length(dataset)],
                                                         trait[1:length(trait)])
  ))
  
  # --- Compute results per periode dataset per trait
  for(i in 1:length(dataset)){
    for(j in 1:length(trait.test)){
      p_result_drift[i, j] <- compute_p(data = get(dataset[i]),
                                        trait = trait.test[j],
                                        test = "drift",
                                        iteration = iteration)
    }
  }
  
  # --- SIGN OF EVOLUTION (PROBA EBV > 0)
  # matrix for result
  p_result_evol <- as.data.frame(matrix(NA, ncol = 2, nrow = 3,
                                        dimnames = list(dataset[1:length(dataset)],
                                                        trait[1:length(trait)])
  ))
  
  # --- Compute results per periode dataset per trait
  for(i in 1:length(dataset)){
    for(j in 1:length(trait.test)){
      p_result_evol[i, j] <- compute_p(data = get(dataset[i]),
                                       trait = trait.test[j],
                                       test = "sign_evol",
                                       iteration = iteration)
    }
  }
  
  
  # PLOT --------------------------------------------------------------------
  # ... To compare ebv to drift, as there is negative number, we take the absolute difference  
  # ... between ebv and drift, and then multiply by the sign of the maximum
  #Formatting dataframe for ggplot
  # --- Compute ebv - drift for each periode
  # pooled
  post_df["delta"] <- abs(post_df$ebv) - abs(post_df$ebv_drift)
  sign_max <- apply(X = post_df[c("ebv", "ebv_drift")], MARGIN = 1, FUN = function(x) sign(max(x[1], x[2])))  # extract sign of the max per row
  post_df["delta"] <- post_df["delta"] * sign_max
  post_df["period"] <- rep("Pooled", iteration)
  
  # --- P1
  post_df_periode1["delta"] <- abs(post_df_periode1$ebv) - abs(post_df_periode1$ebv_drift)
  sign_max <- apply(X = post_df_periode1[c("ebv", "ebv_drift")], MARGIN = 1, FUN = function(x) sign(max(x[1], x[2])))  # extract sign of the max per row
  post_df_periode1["delta"] <- post_df_periode1["delta"] * sign_max
  post_df_periode1["period"] <- rep("Period 1", iteration)
  
  
  post_df_periode2["delta"] <- abs(post_df_periode2$ebv) - abs(post_df_periode2$ebv_drift)
  sign_max <- apply(X = post_df_periode2[c("ebv", "ebv_drift")], MARGIN = 1, FUN = function(x) sign(max(x[1], x[2])))  # extract sign of the max per row
  post_df_periode2["delta"] <- post_df_periode2["delta"] * sign_max
  post_df_periode2["period"] <- rep("Period 2", iteration)
  

  

  
  # --- compute mean and IC of slopes
  # compute_mean_ic <- function(vector, decimal = decimal){
  #   moyenne <- round(mean(vector), decimal)
  #   ic_lwr <-  round(qnorm(p = 0.025, mean = moyenne, sd = sd(vector) / sqrt(length(vector))), decimal)
  #   ic_upr <-  round(qnorm(p = 0.975, mean = moyenne, sd = sd(vector) / sqrt(length(vector))), decimal)
  #   return(paste0(moyenne, " [", ic_lwr, "; ", ic_upr, "]"))
  # }
  
  compute_mean_ic <- function(vector, decimal = decimal){
    moyenne <- round(median(vector), decimal)
    ic_lwr <-  round(quantile(vector, probs = 0.05), decimal)
    ic_upr <-  round(quantile(vector, probs = 0.95), decimal)
    return(paste0(moyenne, " [", ic_lwr, "; ", ic_upr, "]"))
  }
  
  # --- for loop for all trait
  col.names.temp <- c()
  for(trait.temp in trait){
    col.names.temp <- c(col.names.temp, paste0(c("ebv_", "ebv_drift_", "delta_"), trait.temp))
  }
  
  df_stats <- as.data.frame(matrix(NA, ncol = length(trait)*3, nrow = 3, dimnames = list(c("Pooled", "P1", "P2"), col.names.temp)))  #df to store result
  
  col.index <- 1  # index for column                          
  for(trait.temp in trait){
    df_stats[1, col.index] <- compute_mean_ic(vector = post_df$ebv[post_df$trait == trait.temp], decimal = 3)
    df_stats[1, col.index+1] <- compute_mean_ic(vector = post_df$ebv_drift[post_df$trait == trait.temp], decimal = 3)
    df_stats[1, col.index+2] <- compute_mean_ic(vector = post_df$delta[post_df$trait == trait.temp], decimal = 3)
    
    df_stats[2, col.index] <- compute_mean_ic(vector = post_df_periode1$ebv[post_df_periode1$trait == trait.temp], decimal = 3)
    df_stats[2, col.index+1] <- compute_mean_ic(vector = post_df_periode1$ebv_drift[post_df_periode1$trait == trait.temp], decimal = 3)
    df_stats[2, col.index+2] <- compute_mean_ic(vector = post_df_periode1$delta[post_df_periode1$trait == trait.temp], decimal = 3)
    
    df_stats[3, col.index] <- compute_mean_ic(vector = post_df_periode2$ebv[post_df_periode2$trait == trait.temp], decimal = 3)
    df_stats[3, col.index+1] <- compute_mean_ic(vector = post_df_periode2$ebv_drift[post_df_periode2$trait == trait.temp], decimal = 3)
    df_stats[3, col.index+2] <- compute_mean_ic(vector = post_df_periode2$delta[post_df_periode2$trait == trait.temp], decimal = 3)
    
    col.index = col.index + 3
  }                          

  
  
  # --- Create dataframe gathering results per periode
  df_tot <- rbind(post_df, post_df_periode1, post_df_periode2)
  df_tot$period <- factor(df_tot$period, levels = rev(c("Pooled", "Period 1", "Period 2")))
  df_tot$trait <- factor(df_tot$trait, levels = c("masse", "tibia"))
  
  # --- Formatting for ggplot
  df_tot_gg <- pivot_longer(data = df_tot, cols = c("ebv", "delta"), values_to = "y", names_to = "condition")
  df_tot_gg$condition <- factor(df_tot_gg$condition, levels = rev(c("delta", "ebv")))
  
  # GOOG TUTO FOR ggdist::stat_halfeye
  #'[https://cran.r-project.org/web/packages/ggdist/vignettes/slabinterval.html]
  
  # # FACETS ylim TO BE CENTERED ON 0
  # # --- tibia
  # if(abs(min(df_tot$delta[df_tot$trait == "tibia"])) > abs(max(df_tot$delta[df_tot$trait == "tibia"]))){
  #   lim_tibia.temp <- abs(min(df_tot$delta[df_tot$trait == "tibia"]))
  # } else {
  #   lim_tibia.temp <- abs(max(df_tot$delta[df_tot$trait == "tibia"]))
  # }
  # 
  # # --- masse
  # if(abs(min(df_tot$delta[df_tot$trait == "masse"])) > abs(max(df_tot$delta[df_tot$trait == "masse"]))){
  #   lim_masse.temp <- abs(min(df_tot$delta[df_tot$trait == "masse"]))
  # } else {
  #   lim_masse.temp <- abs(max(df_tot$delta[df_tot$trait == "masse"]))
  # }
  # 
  # # --- set limits to object for ggplot
  # lim_tibia <- c(-lim_tibia.temp, lim_tibia.temp)
  # lim_masse <- c(-lim_masse.temp, lim_masse.temp)
  
  
  
  # GGPLOT
 plot_gg <-    ggplot(data = df_tot_gg, aes(x = period, y = y, fill = period)) +
   stat_halfeye(
     mapping = aes(fill_ramp = after_stat(y > 0)),
     justification = -.2,  # Move geom to the right to have space for boxplot
     .width = 0,   # remove default boxplot
     point_colour = NA,  # remove default boxplot
     slab_color = "black",  # line around density curves
     slab_size = 0.6,  # thickness of line around density curves
     alpha = 0.8,
     expand = F,  # expand lines below plot to min and max
     adjust = 0.8,  # smothness of curves (max to 1)
     normalize = "panels"  # density of each plot varying between 0 and 1
   ) +

   # Colour for the different periodes
   scale_fill_manual(values = c(color_period[3], color_period[2], color_period[1])) +
   # Colour per condition : bv < 0 or bv > 0
   scale_fill_ramp_discrete(na.translate = F) +
   # Pic at 0
   # stat_spike(at = 0, linewidth = 0.5) +

   geom_boxplot(width = .12,
                linewidth = 0.5,  # linewidth of boxplot
                outlier.color = NA,  # if NA, do not show outliers
                alpha = 0.7) +
   coord_flip() +
   # Facets and parameters
   facet_grid(condition ~ trait,
              scales = "free",  # allow each facet to have its own y and x axis
              labeller = labeller(trait = c("masse" = "Mass", "tibia" = "Tibia length"),
                                  condition = c(ebv = "EBV", delta = "EBV - DRIFT"))  # name for facets panels
   ) +


   geom_hline(yintercept = 0, size = 0.5) +  # add black line of 0

   # Set ylab for each facets with <ggh4x> according above the ggplot script
   # ggh4x::facetted_pos_scales(y = list(trait == "tibia" ~ scale_y_continuous(limits = lim_tibia),
   #                                     trait == "masse" ~ scale_y_continuous(limits = lim_masse))) +

   # labs(y = latex2exp::TeX("$\\Delta (BV_{predicted} \\ , \\ BV_{drift})$"),
   #     x = "Density") +
   labs(y = "Slope estimates",
        x = "Probability density") +

   # OVERALL THEME
   tidyquant::theme_tq() +

   # Annotate p-values

   # THEME
   theme(text = element_text(size = 20),
         legend.position = "none",
         strip.background.x = element_rect(fill = "#4F4F4F", color = "#4F4F4F", linewidth = 0.5),  # fill of panels x
         strip.background.y = element_rect(fill = "#4F4F4F", color = "#4F4F4F", linewidth = 0.5),  # fill of panels x
         strip.text = element_text(color = "white", size = 15, face = "bold"),  # text color of panels
         strip.text.x = element_text(margin = margin (3, 3, 5, 3)),  # space around x panels text (t, r, b, l)
         strip.text.y = element_text(margin = margin (3, 4, 3, 4)),  # space around y panels text
         axis.ticks.y = element_blank(),
         axis.text.x = element_text(margin = margin(t = 5)), # Adjust padding at the bottom of x-axis labels
         axis.title.x = element_text(margin = margin(t = 10)), # Adjust padding at the bottom of x-axis title
         axis.title.y = element_text(size = 22, margin = margin(0, 9, 0, 0)),  # space between ylab and plot

   )

  # DISPLAY PLOT
  gridExtra::grid.arrange(plot_gg)

  
  # RETURN P VALUES
  return(list(bv = df_tot, sign_evolution_p_values = p_result_evol, drift_p_values = p_result_drift, df_stats = df_stats))
  
}





# COMPUTE SELECTION  ------------------------------------------
# COMPARE OBSERVED VERSUS DRIFT FOR BIVARIATE MODELS -------------------------------------------------------
#' @param mod_result - MCMCglmm model result
#' @param trait - trait for which to compute selection
#' @param fitness - name of fitness variable
#' @param decimal - precision of parameters

#' @WARNINGS : trait must be entered in the same order as in model
#' @return Return a dataframe containing Va(w), cov(w, trait) for each RE, selection differential and gradient of selection for phenotype and genetic

compute_selection <- function(
    mod_result = mod_result,
    trait = c("tibia", "masse"),
    fitness = c("LRS_relative"),
    decimal = 2
) {
  # --- store result
  df <- data.frame(matrix(NA, ncol = length(trait),
                          nrow = 8,
                          dimnames = list(c("Va(w)", "cov_Va(trait, w)", "cov_Vm(trait, w)", "cov_Vr(trait, w)", "S_pheno", "B_pheno", "S_gen", "B_gen"), trait[1:length(trait)])))
  
  # mod result for random effects
  vcv <- mod_result$VCV
  
  # GET VARIANCE ASSOCIATED WITH FIXED EFFECT (De Villemereuil, 2018)
  compute_vcvpred <- function(beta, design_matrix, ntraits) { 
    
    list(cov(matrix(design_matrix %*% beta, ncol = ntraits))) 
    
  }

  X <- mod_result$X  # desgin matrix
  X_mat <- as.matrix(X)  # convert to matrix
  col.names.temp <- colnames(mod_result$Sol)  # colnames of prediction
  
  # get index associating with last fixed effect
  last_row <- grep(pattern = "animal", col.names.temp)[1] - 1
  last_fixed_effect_index <- ifelse(!is.na(last_row), last_row, length(col.names.temp)) # needs to adjust because fixed effects change according to period modeled
  sol <- mod_result$Sol[,1:last_fixed_effect_index]  # select only fixed effects
  
  fixed_effects_var <- flatten( apply(sol, 1, compute_vcvpred, design_matrix = X, ntraits = length(trait) +1) )  # get variance-covariance matrix associated with prediction
  
  for(trait.temp in trait){
    vf.temp <- unlist(lapply(X = fixed_effects_var, FUN = function(x) x[which(trait.temp == trait), which(trait.temp == trait)]))
    assign(paste0("vf_", trait.temp), vf.temp)
  }  # extract variance associated with each trait
  
  
  # COMPUTE SELECTION PARAMETERS
  for(trait.temp in trait){
    
    
    # RESULTS OUTPUT 
    # --- Additive variance of fitness
    va_w <- vcv[,paste0("trait", fitness, ":trait", fitness ,".animal")]  
    
    # --- additive genetic covariance (trait, fitness)
    gen_cov_w <- vcv[,paste0("trait", fitness, ":trait", trait.temp ,".animal")]  
    
    # --- maternal covariance (trait, fitness)
    mother_cov_w <- vcv[,paste0("trait", fitness, ":trait", trait.temp ,".mother")]  
    
    # --- environmental covariance (trait, fitness)
    env_cov_w <- vcv[,paste0("trait", fitness, ":trait", trait.temp ,".units")]  
    
    # --- Numerator phenotypic differential selection
    S.phen <- gen_cov_w + env_cov_w + mother_cov_w
  
    # --- Numerator phenitypic differential selection
    S.gen <- gen_cov_w
    
    # --- Total phenotypic variance per trait
      if(length(trait) > 1){
        vp <- vcv[,paste0("trait", trait.temp, ":trait", trait.temp ,".animal")] +   # additive variance
          vcv[,paste0("trait", trait.temp, ":trait", trait.temp ,".mother")] +  # maternal effect variance
          vcv[,paste0("at.level(trait, c(\"",trait[1], "\", \"", trait[2], "\"))", which(trait == trait[1]), ":at.level(trait, c(\"", trait[1], "\", \"", trait[2], "\"))", which(trait == trait[1]),".annee_emergence")] +  # cohort variance
          vcv[,paste0("trait", trait.temp, ":trait", trait.temp ,".units")] +  # residuals variance
          get(paste0("vf_", trait.temp))  # variance from fixed effect
        
        
      } else if (length(trait) == 1){
        vp <- vcv[,paste0("trait", trait.temp, ":trait", trait.temp ,".animal")] +   # additive variance
          vcv[,paste0("trait", trait.temp, ":trait", trait.temp ,".mother")] +  # maternal effect variance
          vcv[,paste0("at.level(trait, c(\"",trait, "\")):at.level(trait, c(\"", trait, "\")).annee_emergence")] +  # cohort variance
          vcv[,paste0("trait", trait.temp, ":trait", trait.temp ,".units")] +  # residuals variance
          get(paste0("vf_", trait.temp))  # variance from fixed effect
        
      }
      

    # COMPUTE PARAMETERS
    compute_mean_ic <- function(vector, decimal = decimal){
        moyenne <- round(mean(vector), decimal)
        ic_lwr <-  round(HPDinterval(vector)[1], decimal)
        ic_upr <- round(HPDinterval(vector)[2], decimal)
      return(paste0(moyenne, " [", ic_lwr, "; ", ic_upr, "]"))
    }
    
    # Phenotypic differential selection
    df["Va(w)", which(trait == trait.temp)] <- compute_mean_ic(va_w, decimal = decimal)
    df["cov_Va(trait, w)", which(trait == trait.temp)] <- compute_mean_ic(gen_cov_w, decimal = decimal)
    df["cov_Vm(trait, w)", which(trait == trait.temp)] <- compute_mean_ic(mother_cov_w, decimal = decimal)
    df["cov_Vr(trait, w)", which(trait == trait.temp)] <- compute_mean_ic(env_cov_w, decimal = decimal)
    df["S_pheno", which(trait == trait.temp)] <- compute_mean_ic(S.phen, decimal = decimal)
    df["B_pheno", which(trait == trait.temp)] <- compute_mean_ic(S.phen / vp, decimal = decimal)
    df["S_gen", which(trait == trait.temp)] <- compute_mean_ic(S.gen, decimal = decimal)
    df["B_gen", which(trait == trait.temp)] <- compute_mean_ic(S.gen / vp, decimal = decimal)
    
  }
  return(df)
}



# COMPUTE ADDITIVE COVARIANCE  ------------------------------------------
#' @param mod_result - MCMCglmm FITNESS model result
#' @param trait - trait for which to compute selection
#' @param fitness - name of fitness variable
#' @param decimal - precision of parameters

#' @WARNINGS : trait must be entered in the same order as in model
#' @return Return a dataframe containing the (co)variance additive genetic between mass, tibia and relative fitness from a trivariate model

compute_additive_cov <- function(
    mod_result = mod_result_fitness,
    trait = c("tibia", "masse"),
    fitness = c("LRS_relative"),
    decimal = 2
) {
  # mod result for random effects
  vcv <- mod_result$VCV
  sol <- mod_result$Sol
  
  
  
  #  RECONSTRUCITNG VCV OBJECT AND SEPARATE PERIOD 1 AND PERIOD 2
  # column containing result of RE per iteration per individual
  col.name.temp <- colnames(sol)
  col.name.vcv <- colnames(vcv)
  
  # Filter $Sol for breedinvcv# Filter $Sol for breeding values only
  sol <- sol[,grep(pattern = "animal", x = col.name.temp)]
  
  # Filter $VCV for .animal only
  vcv <- vcv[,grep(pattern = "animal", x = col.name.vcv)]
  
  # reconstruct vcv object
  df_vcv <- as.data.frame(matrix(NA, nrow = nrow(vcv), ncol = ncol(vcv), dimnames = dimnames(vcv)))
  
  
  # ASSOCIATE COHORT PER INDIVIDUAL
  # --- Extract individual ID associated with EBV
  col.index.temp <-  grep(pattern = paste0(trait[1], ".animal"), x = col.name.temp)
  id_ebv <- gsub(pattern = paste0("trait", trait[1],".animal.", collapse = ""), replacement = "", x = col.name.temp[col.index.temp])  # remove extra text to get individual id
  df_bv <- data.frame(individual_id = id_ebv)
  
  # --- Extract cohort date associated with id
  df_bv <- merge(x = pedigree_data[c("individual_id", "cohort")], y = df_bv, by = "individual_id")
  
  # GET ID PER PERIOD
  id_pooled <- df_bv[df_bv$cohort %in% 1997:2016, "individual_id"]
  id_P1 <- df_bv[df_bv$cohort %in% 1997:2005, "individual_id"]
  id_P2 <- df_bv[df_bv$cohort %in% 2006:2016, "individual_id"]  # 2016 because of the fitness age threshold
  
  
  # RECONSTRUCTING VARIANCE / CO-VARIANCE EFFECT
  # --- ADDITIVE AND MATERNAL VARIANCE
  get_additive_cov <- function(period_id = id_P2, df_vcv = df_vcv){
    
    for(trait.temp1 in c(trait, fitness)){
      for(trait.temp2 in c(trait, fitness)){  # include a second trait to calculate covariance as well as variance (as cov(x, x) = var(x))
        
        # store filtered result per trait
        df1 <- as.data.frame(sol[,grep(pattern = trait.temp1, x = colnames(sol))])
        df2 <- as.data.frame(sol[,grep(pattern = trait.temp2, x = colnames(sol))])
        
        # renaming column names with individual ID
        colnames(df1) <- gsub(pattern = paste0("trait", trait.temp1,".animal.", collapse = ""), replacement = "", x = colnames(df1))
        colnames(df2) <- gsub(pattern = paste0("trait", trait.temp2,".animal.", collapse = ""), replacement = "", x = colnames(df1))
        
        # Filter by period
        df1 <- df1[,colnames(df1) %in% period_id]
        df2 <- df2[,colnames(df2) %in% period_id]
        
        # compute covariance between paire of row
        res.temp <- 
          apply(cbind(df1, df2), 1, function(x){
            cov(x[1:length(df1)],
                x[(length(df2)+1):(length(df1)*2)])  # *2 to go at the end of dataframe
          }
          )
        
        df_vcv[,paste0("trait", trait.temp1, ":trait", trait.temp2, ".animal")] <- res.temp
        
      }
    }
    
    return(df_vcv)
  }
  
  vcv_pooled <- get_additive_cov(period_id = id_pooled, df_vcv = df_vcv)
  vcv_P1 <- get_additive_cov(period_id = id_P1, df_vcv = df_vcv)
  vcv_P2 <- get_additive_cov(period_id = id_P2, df_vcv = df_vcv)
  
  return(list(vcv_pooled = vcv_pooled, vcv_P1 = vcv_P1, vcv_P2 = vcv_P2))
  
  # # --- COHORT EFFECT
  # col.cohort <- grep(pattern = "emergence", col.name.temp)
  # 
  #     
  #     # Get name associated with cohort effects
  #     cohort_name <- function(x, trait = trait){
  #       name.var.temp <- paste0(which(trait == x), ".annee_emergence")  # the "1" or "2" yielded by which(trait == x) corresponds to the index for trait in cohort [see summary]
  #       return(name.var.temp)
  #     }
  #     
  #     
  #     for(trait.temp1 in trait){
  #       for(trait.temp2 in trait){
  #         # temporary df for apply() function formatting
  #         col.size.indexes1 <- grep(pattern = cohort_name(x = trait.temp1, trait = c("tibia", "masse")), x = col.name.temp)  # get sol associared with trait 1
  #         col.size.indexes2 <- grep(pattern = cohort_name(x = trait.temp2, trait = c("tibia", "masse")), x = col.name.temp)  # get sol associated with trait 2
  #         
  #         df1 <- as.data.frame(sol[,col.size.indexes1])
  #         df2 <- as.data.frame(sol[,col.size.indexes2])
  #         
  #         # compute covariance between paire of row
  #         res.temp <- 
  #           apply(cbind(df1, df2), 1, function(x){
  #             cov(x[1:length(df1)],
  #                 x[(length(df2)+1):(length(df1)*2)])  # *2 to go at the ned of dataframe
  #           })
  #             # fill vcv
  #             name_cohort.temp <- name.var.temp <- paste0("at.level(trait, c(\"",trait[1], "\", \"", trait[2], "\"))", which(trait == trait.temp1), ":at.level(trait, c(\"", trait[1], "\", \"", trait[2], "\"))", which(trait == trait.temp2),".annee_emergence")
  # 
  #             df_vcv[, name_cohort.temp] <- res.temp
  #           
  #       }
  #     }
  
  
}











# COMPUTE MBE  ------------------------------------------
#' @param GenerationTime - Generation time of the species to compute response to selection
#' @param mod_result - MCMCglmm model result for bivariate model to obtain variance covariance matrix of both trait
#' @param mod_result_fitness - MCMCglmm model result from trivariate model (mass, tibia, relative fitness)
#' @param mod_result_fitness_pheno - MCMCglmm model result from trivariate model (mass, tibia, relative fitness) with only residuals as random effect
#' @param additive_covariance_matrix - additive_covariance_matrix reconstructed from trivariate mod_result_fitness
#' @param decimal - precision of estimates
#' @param study_period_tibia - length of the study period for tibia length
#' @param study_period_masse - length of the study period for masse

#' @WARNINGS : can be disfunctioning according to trait order specification in model (to be corrected)
#' @return Return a dataframe containing selection differential, selection gradient normalized by variance beta_sigma or variance AND mean beta_mu
#' ... and predicted response from multivariate breeder equation accouting for genetic correlation (CG) or not (no_CG)
#' 
compute_selection_MBE <- function(GenerationTime = 5.9,
                                  mod_result = mod_result,
                                  mod_result_fitness = mod_result_fitness,
                                  mod_result_fitness_pheno = mod_result_fitness_pheno,
                                  additive_covariance_matrix = compute_additive_cov(mod_result = mod_result_fitness)$vcv_pooled,
                                  decimal = 3,
                                  study_period_tibia = study_period_tibia,
                                  study_period_masse = study_period_masse
)
{
  
  # GET VARIANCE ASSOCIATED WITH FIXED EFFECT (De Villemereuil, 2018) TO ADD IT TO PHENOTYPE VARIANCE
  compute_vcvpred <- function(beta, design_matrix, ntraits) { 
    
    list(cov(matrix(design_matrix %*% beta, ncol = ntraits))) 
    
  }
  
  X <- mod_result_fitness_pheno$X  # desgin matrix
  X_mat <- as.matrix(X)  # convert to matrix
  col.names.temp <- colnames(mod_result_fitness_pheno$Sol)  # colnames of prediction
  
  # get index associating with last fixed effect
  last_row <- grep(pattern = "animal", col.names.temp)[1] - 1
  last_fixed_effect_index <- ifelse(!is.na(last_row), last_row, length(col.names.temp)) # needs to adjust because fixed effects change according to period modeled
  sol <- mod_result_fitness_pheno$Sol[,1:last_fixed_effect_index]  # select only fixed effects
  
  fixed_effects_var <- flatten( apply(sol, 1, compute_vcvpred, design_matrix = X, ntraits = 3))  # get variance-covariance matrix associated with prediction
  
  # compute vector of vf for each trait
  vf_masse <- unlist(lapply(X = fixed_effects_var, FUN = function(x) x[2, 2]))
  vf_tibia <- unlist(lapply(X = fixed_effects_var, FUN = function(x) x[1, 1]))
  vf_cov <- unlist(lapply(X = fixed_effects_var, FUN = function(x) x[1, 2]))
  
  
  
  multivarBreederLBS <- list()  # storage of response according to multivariate breeder equation accounting for genetic correlation
  multivarBreederLBSnoCG <- list()  # with no genetic correlation
  Gmat <- list()  # variance covariance additive genetic matrix
  Pmat <- list()  # phenotype variance covariancre matrix
  STS_list <- list()
  cov_mother_list <- list()
  cov_residuals_list <- list()
  Slist <- list()  # selection differential
  S_normlist <- list()  # selection differential normalized
  Gradlist <- list()  # selection gradient beta 
  Grad_normlist <- list()  # selection gradient beta_sigma normalised by total phenotype variance
  
  
  
  
  for (i in 1:length(mod_result_fitness$VCV[,1])){
    
    # ADDITIVE GENETIC VARIANCE COV MATRIX
    G <- matrix(c(mod_result$VCV[i,"traittibia:traittibia.animal"], # tibia tibia
                  mod_result$VCV[i,"traitmasse:traittibia.animal"], # tibia masse
                  mod_result$VCV[i,"traitmasse:traittibia.animal"], # masse tibia
                  mod_result$VCV[i,"traitmasse:traitmasse.animal"]), nrow = 2) # masse masse
    
    
    # SELECTION GRADIENT VARIANCE COV MATRIX
    # ... [Somme de toutes les covariances traits / fitness]
    S <- c(
      # mod_result_fitness$VCV[i,"traittibia:traitLRS_relative.mother"] + 
      #   mod_result_fitness$VCV[i,"traittibia:traitLRS_relative.units"] +
      #   additive_covariance_matrix[i,"traittibia:traitLRS_relative.animal"] +
      mod_result_fitness_pheno$VCV[i,"traittibia:traitLRS_relative.units"],
      
      # mod_result_fitness$VCV[i,"traitmasse:traitLRS_relative.mother"] + 
      #   mod_result_fitness$VCV[i,"traitmasse:traitLRS_relative.units"] +
      #   additive_covariance_matrix[i,"traitmasse:traitLRS_relative.animal"] +
      mod_result_fitness_pheno$VCV[i,"traitmasse:traitLRS_relative.units"]
    )
    
    
    # PHENOTYPIC VARIANCE COVARIANCE MATRIX - Somme de toutes les variance de chaque trait en diagonale et off-diag somme des covariance
    P <- matrix(c(
      # somme des variance tibia
      #mod_result$VCV[i,"traittibia:traittibia.animal"]+
      #mod_result$VCV[i,"traittibia:traittibia.mother"]+
      #mod_result$VCV[i,"traittibia:traittibia.annee_emergence"]+
      vf_tibia[i] +  # variance from fixed effect
      mod_result_fitness_pheno$VCV[i,"traittibia:traittibia.units"],
      #mod_result$VCV[i,"at.level(trait, c(\"tibia\", \"masse\"))1:at.level(trait, c(\"tibia\", \"masse\"))1.annee_emergence"] +
      #mod_result$VCV[i,"traittibia:traittibia.units"],
      
      # somme des covariance tibia / masse
      #mod_result$VCV[i,"traittibia:traitmasse.animal"]+
      # mod_result$VCV[i,"traittibia:traitmasse.mother"]+
      #mod_result$VCV[i,"traittibia:traitmasse.annee_emergence"] +
      vf_cov[i] + # covariance from fix effect
      mod_result_fitness_pheno$VCV[i,"traittibia:traitmasse.units"],
      # mod_result$VCV[i,"at.level(trait, c(\"tibia\", \"masse\"))2:at.level(trait, c(\"tibia\", \"masse\"))1.annee_emergence"] +
      #mod_result$VCV[i,"traittibia:traitmasse.units"],
      
      #mod_result$VCV[i,"traittibia:traitmasse.animal"]+
      # mod_result$VCV[i,"traittibia:traitmasse.mother"]+
      #mod_result$VCV[i,"traittibia:traitmasse.annee_emergence"]+
      vf_cov[i] + # covariance from fix effect
      mod_result_fitness_pheno$VCV[i,"traittibia:traitmasse.units"],
      
      # mod_result$VCV[i,"at.level(trait, c(\"tibia\", \"masse\"))2:at.level(trait, c(\"tibia\", \"masse\"))1.annee_emergence"] +
      #mod_result$VCV[i,"traittibia:traitmasse.units"],
      
      # somme des variance masse
      #mod_result$VCV[i,"traitmasse:traitmasse.animal"]+
      # mod_result$VCV[i,"traitmasse:traitmasse.mother"]+
      #mod_result$VCV[i,"traitmasse:traitmasse.annee_emergence"]+
      vf_masse[i] +  # variance from fixed effect
      mod_result_fitness_pheno$VCV[i,"traitmasse:traitmasse.units"]
      # mod_result$VCV[i,"at.level(trait, c(\"tibia\", \"masse\"))2:at.level(trait, c(\"tibia\", \"masse\"))2.annee_emergence"] +
      # mod_result$VCV[i,"traitmasse:traitmasse.units"]),
    ),
    nrow = 2)
    
    
    
    # Price equation (cov_A(trait, w))
    STS <- matrix(c(additive_covariance_matrix[i,"traittibia:traitLRS_relative.animal"],
                    NA,
                    NA, 
                    additive_covariance_matrix[i,"traitmasse:traitLRS_relative.animal"]),
                  byrow = T, 
                  nrow = 2)
    
    # Covariance with other RE
    # --- maternal effect
    cov_mother <- matrix(c(mod_result_fitness$VCV[i,"traittibia:traitLRS_relative.mother"],
                           NA,
                           NA, 
                           mod_result_fitness$VCV[i,"traitmasse:traitLRS_relative.mother"]),
                         byrow = T, 
                         nrow = 2)
    
    # --- residuals
    cov_residuals <- matrix(c(mod_result_fitness$VCV[i,"traittibia:traitLRS_relative.units"],
                              NA,
                              NA, 
                              mod_result_fitness$VCV[i,"traitmasse:traitLRS_relative.units"]),
                            byrow = T, 
                            nrow = 2)
    
    
    # gradient computation
    Grad_tibia_masse <- S %*% solve(P)  # cov(trait, w) / vP
    Grad_tibia_masse_norm <- Grad_tibia_masse %*%  matrix(c(sqrt(P[1,1]),0, 0, sqrt(P[2,2])), ncol = 2) # beta * sd
    
    Gmat[[i]] <- G
    Pmat[[i]] <- P
    Slist[[i]] <- S
    S_normlist[[i]] <- S  / c(sqrt(P[1,1]), sqrt(P[2,2])) # selection differential standardized
    Gradlist[[i]] <- Grad_tibia_masse
    Grad_normlist[[i]] <- Grad_tibia_masse_norm  #beta standardized
    STS_list[[i]] <- STS
    cov_mother_list[[i]] <- cov_mother
    cov_residuals_list[[i]] <- cov_residuals
    multivarBreederLBS[[i]] <- Grad_tibia_masse %*% G  # with genetic correlation
    multivarBreederLBSnoCG[[i]] <- Grad_tibia_masse %*% (G * diag(2))  # without genetic correlation
  }
  
  # 
  # # P MATRIX
  Pmode <- matrix(c(posterior.mode(as.mcmc(unlist(lapply(Pmat, function(x){x[1,1]})))),
                    posterior.mode(as.mcmc(unlist(lapply(Pmat, function(x){x[1,2]})))),
                    posterior.mode(as.mcmc(unlist(lapply(Pmat, function(x){x[2,1]})))),
                    posterior.mode(as.mcmc(unlist(lapply(Pmat, function(x){x[2,2]}))))),
                  nrow=2)
  # Pmode
  # 
  # #phenotypic correlation
  # posterior.mode(as.mcmc(unlist(lapply(Pmat, function(x){
  #   x[1,2]/sqrt(x[1,1]*x[2,2])
  #   
  # })) ) )
  # 
  # HPDinterval(as.mcmc(unlist(lapply(Pmat, function(x){
  #   x[1,2]/sqrt(x[1,1]*x[2,2])
  #   
  # })) ) )
  # 
  # # G MATRIX
  Gmode <- matrix(c(mean(as.mcmc(unlist(lapply(Gmat, function(x){x[1,1]})))),
                    mean(as.mcmc(unlist(lapply(Gmat, function(x){x[1,2]})))),
                    mean(as.mcmc(unlist(lapply(Gmat, function(x){x[2,1]})))),
                    mean(as.mcmc(unlist(lapply(Gmat, function(x){x[2,2]}))))),
                  nrow=2)
  
  # S MATRIX
  Smode <- matrix(c(mean(as.mcmc(unlist(lapply(Slist, function(x){x[1]})))),
                    mean(as.mcmc(unlist(lapply(Slist, function(x){x[2]}))))),
                  nrow=1)
  
  # Gmode
  # 
  # # genetic covariance
  # mean(as.mcmc(unlist(lapply(Gmat, function(x){x[1,2]})) ) )
  # HPDinterval(as.mcmc(unlist(lapply(Gmat, function(x){x[1,2]})) ) )
  # 
  # # genetic correlation
  # mean(as.mcmc(unlist(lapply(Gmat, function(x){x[1,2]/sqrt(x[1,1]*x[2,2])})) ) )
  # HPDinterval(as.mcmc(unlist(lapply(Gmat, function(x){x[1,2]/sqrt(x[1,1]*x[2,2])})) ) )
  # mean(as.mcmc(unlist(lapply(Gmat, function(x){x[1,2]/sqrt(x[1,1]*x[2,2])})) ) >0)
  
  ####################
  # CORRELATION 
  ###################
  PM_cor_pheno <- round(mean(as.mcmc(unlist(lapply(Pmat, function(x){x[2,1] / sqrt(x[1,1] * x[2,2])})))) , decimal)  # phenotype
  HPD_cor_pheno <- round(HPDinterval(as.mcmc(unlist(lapply(Pmat, function(x){x[2,1] / sqrt(x[1,1] * x[2,2])})))) , decimal)  # phenotpe
  
  PM_cor_gen <- round(mean(as.mcmc(unlist(lapply(Gmat, function(x){x[2,1] / sqrt(x[1,1] * x[2,2])})))) , decimal)  # additive genetic
  HPD_cor_gen <- round(HPDinterval(as.mcmc(unlist(lapply(Gmat, function(x){x[2,1] / sqrt(x[1,1] * x[2,2])})))) , decimal)  # additive genetic
  
  
  
  #####################
  # VA AND HERITABILITY
  #####################
  # additive genetic variance in tibia
  PM_VA_tibia <- round(mean(as.mcmc(unlist(lapply(Gmat, function(x){x[1,1]})) ) ), decimal)
  HPD_VA_tibia <- round(HPDinterval(as.mcmc(unlist(lapply(Gmat, function(x){x[1,1]})) ) ), decimal)
  
  
  # heritability of tibia
  PM_H_tibia <- round(mean(as.mcmc(unlist(lapply(Gmat, function(x){x[1,1]})) ) /
                             as.mcmc(unlist(lapply(Pmat, function(x){x[1,1]})) )), decimal)
  
  HPD_H_tibia <-   round(HPDinterval(as.mcmc(unlist(lapply(Gmat, function(x){x[1,1]})) ) /
                                       as.mcmc(unlist(lapply(Pmat, function(x){x[1,1]})) )), decimal)
  
  
  
  # additive genetic variance in mass
  PM_VA_masse <-  round(mean(as.mcmc(unlist(lapply(Gmat, function(x){x[2,2]})) ) ), decimal)
  HPD_VA_masse <- round(HPDinterval(as.mcmc(unlist(lapply(Gmat, function(x){x[2,2]})) ) ), decimal)
  
  
  # heritability of masse
  PM_H_masse <- round(mean(as.mcmc(unlist(lapply(Gmat, function(x){x[2,2]})) ) /
                             as.mcmc(unlist(lapply(Pmat, function(x){x[2,2]})) )), decimal)
  
  HPD_H_masse <-   round(HPDinterval(as.mcmc(unlist(lapply(Gmat, function(x){x[2,2]})) ) /
                                       as.mcmc(unlist(lapply(Pmat, function(x){x[2,2]})) )), decimal)
  
  
  
  ##################### ##############
  # STS AND OTHER COV(TRAIT, FITNESS)
  ##################### ##############
  # TIBIA
  # sts
  STS_tibia <- as.mcmc(unlist(lapply(STS_list, function(x){x[1,1]})) )
  PM_STS_tibia <- round(mean(STS_tibia ), decimal)
  HPD_STS_tibia <- round(HPDinterval(STS_tibia ), decimal)
  Prob_STS_tibia_pos <- mean(STS_tibia<0)
  
  # mother
  cov_mother_tibia <- as.mcmc(unlist(lapply(cov_mother_list, function(x){x[1,1]})) )
  PM_cov_mother_tibia <- round(mean(cov_mother_tibia ), decimal)
  HPD_cov_mother_tibia <- round(HPDinterval(cov_mother_tibia ), decimal)
  Prob_cov_mother_tibia_pos <- mean(cov_mother_tibia<0)
  
  # residuals
  cov_residuals_tibia <- as.mcmc(unlist(lapply(cov_residuals_list, function(x){x[1,1]})) )
  PM_cov_residuals_tibia <- round(mean(cov_residuals_tibia ), decimal)
  HPD_cov_residuals_tibia <- round(HPDinterval(cov_residuals_tibia ), decimal)
  Prob_cov_residuals_tibia_pos <- mean(cov_residuals_tibia<0)
  
  
  # MASS
  # sts
  STS_masse <- as.mcmc(unlist(lapply(STS_list, function(x){x[2,2]})) )
  PM_STS_masse <- round(mean(STS_masse), decimal)
  HPD_STS_masse <- round(HPDinterval(STS_masse), decimal)
  Prob_STS_masse_pos <- mean(STS_masse<0)
  
  # mother
  cov_mother_masse <- as.mcmc(unlist(lapply(cov_mother_list, function(x){x[2,2]})) )
  PM_cov_mother_masse <- round(mean(cov_mother_masse), decimal)
  HPD_cov_mother_masse <- round(HPDinterval(cov_mother_masse), decimal)
  Prob_cov_mother_masse_pos <- mean(cov_mother_masse<0)
  
  # residuals
  cov_residuals_masse <- as.mcmc(unlist(lapply(cov_residuals_list, function(x){x[2,2]})) )
  PM_cov_residuals_masse <- round(mean(cov_residuals_masse), decimal)
  HPD_cov_residuals_masse <- round(HPDinterval(cov_residuals_masse), decimal)
  Prob_cov_residuals_masse_pos <- mean(cov_residuals_masse<0)
  
  
  ##########################################
  # FUNDAMENTAL THEOREM OF NATURAL SELECTION
  ##########################################
  PM_var_additive_fitness <- round(mean(additive_covariance_matrix$`traitLRS_relative:traitLRS_relative.animal`), 2)
  HPD_var_additive_fitness <- round(HPDinterval(as.mcmc(additive_covariance_matrix$`traitLRS_relative:traitLRS_relative.animal`)), 2)
  
  
  ##########################################################
  # SELECTION DIFFERENTIAL
  ##########################################################
  # PER YEAR
  # --- tibia
  S_tibia_year <- as.mcmc(as.numeric(unlist(lapply(Slist, function(x){x[1]})))) / GenerationTime
  PM_S_tibia_year <- round(mean(S_tibia_year), decimal)
  HPD_S_tibia_year <- round(HPDinterval(S_tibia_year), decimal) 
  Prob_S_tibia_pos_year <- mean(S_tibia_year<0)  # probability S > 0
  
  # --- masse
  S_masse_year <- as.mcmc(as.numeric(unlist(lapply(Slist, function(x){x[2]})))) / GenerationTime
  PM_S_masse_year <- round(mean(S_masse_year), decimal)
  HPD_S_masse_year <- round(HPDinterval(S_masse_year), decimal) 
  Prob_S_masse_pos_year <- mean(S_masse_year<0)  # probability S > 0
  
  
  # PER GENERATION
  # --- tibia
  S_tibia_generation <- as.mcmc(as.numeric(unlist(lapply(Slist, function(x){x[1]} ))))
  PM_S_tibia_generation <- round(mean(S_tibia_generation), decimal)
  HPD_S_tibia_generation <- round(HPDinterval(S_tibia_generation), decimal) 
  Prob_S_tibia_pos_generation <- mean(S_tibia_generation<0)  # probability S > 0
  
  # --- masse
  S_masse_generation <- as.mcmc(as.numeric(unlist(lapply(Slist, function(x){x[2]} ))))
  PM_S_masse_generation <- round(mean(S_masse_generation), decimal)
  HPD_S_masse_generation <- round(HPDinterval(S_masse_generation), decimal) 
  Prob_S_masse_pos_generation <- mean(S_masse_generation<0)  # probability S > 0
  
  
  ##########################################################
  # SELECTION DIFFERENTIAL (normalized sd)
  ##########################################################
  # PER YEAR
  # --- tibia
  S_norm_tibia_year <- as.mcmc(as.numeric(unlist(lapply(S_normlist, function(x){x[1]})))) / GenerationTime
  PM_S_norm_tibia_year <- round(mean(S_norm_tibia_year), decimal)
  HPD_S_norm_tibia_year <- round(HPDinterval(S_norm_tibia_year), decimal) 
  Prob_S_norm_tibia_pos_year <- mean(S_norm_tibia_year<0)  # probability S > 0
  
  # --- masse
  S_norm_masse_year <- as.mcmc(as.numeric(unlist(lapply(S_normlist, function(x){x[2]})))) / GenerationTime
  PM_S_norm_masse_year <- round(mean(S_norm_masse_year), decimal)
  HPD_S_norm_masse_year <- round(HPDinterval(S_norm_masse_year), decimal) 
  Prob_S_norm_masse_pos_year <- mean(S_norm_masse_year<0)  # probability S > 0
  
  
  # PER GENERATION
  # --- tibia
  S_norm_tibia_generation <- as.mcmc(as.numeric(unlist(lapply(S_normlist, function(x){x[1]} ))))
  PM_S_norm_tibia_generation <- round(mean(S_norm_tibia_generation), decimal)
  HPD_S_norm_tibia_generation <- round(HPDinterval(S_norm_tibia_generation), decimal) 
  
  # --- masse
  S_norm_masse_generation <- as.mcmc(as.numeric(unlist(lapply(S_normlist, function(x){x[2]} ))))
  PM_S_norm_masse_generation <- round(mean(S_norm_masse_generation), decimal)
  HPD_S_norm_masse_generation <- round(HPDinterval(S_norm_masse_generation), decimal) 
  
  
  
  ##########################################################
  # SELECTION GRADIENTS (beta)
  ##########################################################
  # PER YEAR
  # --- tibia
  grad_tibia_year <- as.mcmc(as.numeric(unlist(lapply(Gradlist, function(x){x[1]} )))) /
    (GenerationTime)
  PM_beta_tibia_year <- round(mean(grad_tibia_year), decimal)
  HPD_beta_tibia_year <- round(HPDinterval(grad_tibia_year), decimal) 
  Prob_beta_tibia_pos_year <- mean(grad_tibia_year<0)  # probability S > 0
  
  
  # --- masse
  grad_masse_year <- as.mcmc(as.numeric(unlist(lapply(Gradlist, function(x){x[2]} )))) /
    (GenerationTime)
  PM_beta_masse_year <- round(mean(grad_masse_year), decimal)
  HPD_beta_masse_year <- round(HPDinterval(grad_masse_year), decimal) 
  Prob_beta_masse_pos_year <- mean(grad_masse_year<0)  # probability S > 0
  
  
  # PER GENERATION
  # --- tibia
  grad_tibia_generation <- as.mcmc(as.numeric(unlist(lapply(Gradlist, function(x){x[1]} ))))
  PM_beta_tibia_generation <- round(mean(grad_tibia_generation), decimal)
  HPD_beta_tibia_generation <- round(HPDinterval(grad_tibia_generation), decimal) 
  Prob_beta_tibia_pos_generation <- mean(grad_tibia_generation<0)  # probability S > 0
  
  
  # --- masse
  grad_masse_generation <- as.mcmc(as.numeric(unlist(lapply(Gradlist, function(x){x[2]} ))))
  PM_beta_masse_generation <- round(mean(grad_masse_generation), decimal)
  HPD_beta_masse_generation <- round(HPDinterval(grad_masse_generation), decimal) 
  Prob_beta_masse_pos_generation <- mean(grad_masse_generation<0)  # probability S > 0
  
  
  
  
  
  
  ##########################################################
  # SELECTION GRADIENTS BY MEAN (beta_sd)
  ##########################################################
  # PER YEAR
  # --- tibia
  grad_sd_tibia_year <- as.mcmc(as.numeric(unlist(lapply(Grad_normlist, function(x){x[1]} )))) /
    (GenerationTime)
  PM_beta_sd_tibia_year <- round(mean(grad_sd_tibia_year), decimal)
  HPD_beta_sd_tibia_year <- round(HPDinterval(grad_sd_tibia_year), decimal) 
  Prob_beta_sd_tibia_pos_year <- mean(grad_sd_tibia_year<0)  # probability S > 0
  
  
  # --- masse
  grad_sd_masse_year <- as.mcmc(as.numeric(unlist(lapply(Grad_normlist, function(x){x[2]} )))) /
    (GenerationTime)
  PM_beta_sd_masse_year <- round(mean(grad_sd_masse_year), decimal)
  HPD_beta_sd_masse_year <- round(HPDinterval(grad_sd_masse_year), decimal) 
  Prob_beta_sd_masse_pos_year <- mean(grad_sd_masse_year<0)  # probability S > 0
  
  
  # PER GENERATION
  # --- tibia
  grad_sd_tibia_generation <- as.mcmc(as.numeric(unlist(lapply(Grad_normlist, function(x){x[1]} ))))
  PM_beta_sd_tibia_generation <- round(mean(grad_sd_tibia_generation), decimal)
  HPD_beta_sd_tibia_generation <- round(HPDinterval(grad_sd_tibia_generation), decimal) 
  Prob_beta_sd_tibia_pos_generation <- mean(grad_sd_tibia_generation<0)  # probability S > 0
  
  
  # --- masse
  grad_sd_masse_generation <- as.mcmc(as.numeric(unlist(lapply(Grad_normlist, function(x){x[2]} ))))
  PM_beta_sd_masse_generation <- round(mean(grad_sd_masse_generation), decimal)
  HPD_beta_sd_masse_generation <- round(HPDinterval(grad_sd_masse_generation), decimal) 
  Prob_beta_sd_masse_pos_generation <- mean(grad_sd_masse_generation<0)  # probability S > 0
  
  
  
  
  ##########################################################
  # MULTIVARIATE BREEDER EQUATION PREDICTION [predicted response]
  # ... Rmv for "Response multivariate"
  ##########################################################
  
  # EN CONSIDERANT LA CORRELATION GENETIQUE
  # PER YEAR ====================
  # --- tibia
  Rmv_tibia_year <- as.mcmc(unlist(lapply(multivarBreederLBS, FUN = function(x){x[1]}))) /
    (GenerationTime)
  PM_rmv_tibia_year <- round(mean(Rmv_tibia_year, adjust = 1), decimal)
  HPD_rmv_tibia_year <- round(HPDinterval(Rmv_tibia_year), decimal)
  Prob_rmv_tibia_pos_year <- mean(Rmv_tibia_year<0)
  
  # --- masse
  Rmv_masse_year <- as.mcmc(unlist(lapply(multivarBreederLBS, FUN = function(x){x[2]}))) /
    (GenerationTime)
  PM_rmv_masse_year <- round(mean(Rmv_masse_year, adjust = 1), decimal)
  HPD_rmv_masse_year <- round(HPDinterval(Rmv_masse_year), decimal)
  Prob_rmv_masse_pos_year <- mean(Rmv_masse_year<0)
  
  
  
  # PER GENERATION ==========
  Rmv_tibia_generation <- as.mcmc(unlist(lapply(multivarBreederLBS, FUN = function(x){x[1]})))
  PM_rmv_tibia_generation <- round(mean(Rmv_tibia_generation, adjust = 1), decimal)
  HPD_rmv_tibia_generation <- round(HPDinterval(Rmv_tibia_generation), decimal)
  Prob_rmv_tibia_pos_generation <- mean(Rmv_tibia_generation<0)
  
  # --- masse
  Rmv_masse_generation <- as.mcmc(unlist(lapply(multivarBreederLBS, FUN = function(x){x[2]})))
  PM_rmv_masse_generation <- round(mean(Rmv_masse_generation, adjust = 1), decimal)
  HPD_rmv_masse_generation <- round(HPDinterval(Rmv_masse_generation), decimal)
  Prob_rmv_masse_pos_generation <- mean(Rmv_masse_generation<0)
  
  
  
  # OVERALL STUDY ===============
  Rmv_tibia_overall <- as.mcmc(unlist(lapply(multivarBreederLBS, FUN = function(x){x[1]}))) /
    (GenerationTime)
  PM_rmv_tibia_overall <- round(mean(study_period_tibia*Rmv_tibia_overall, adjust = 1), decimal)
  HPD_rmv_tibia_overall <- round(HPDinterval(study_period_tibia*Rmv_tibia_overall), decimal)
  Prob_rmv_tibia_pos_overall <- mean(Rmv_tibia_overall<0)
  
  # --- masse
  Rmv_masse_overall <- as.mcmc(unlist(lapply(multivarBreederLBS, FUN = function(x){x[2]}))) /
    (GenerationTime)
  PM_rmv_masse_overall <- round(mean(study_period_masse*Rmv_masse_overall, adjust = 1), decimal)
  HPD_rmv_masse_overall <- round(HPDinterval(study_period_masse*Rmv_masse_overall), decimal)
  Prob_rmv_masse_pos_overall <- mean(Rmv_masse_overall<0)
  
  
  
  # SANS CONSIDERER LA CORRELATION GENETIQUE
  # PER YEAR ====================
  # --- tibia
  Rmv_noCG_tibia_year <- as.mcmc(unlist(lapply(multivarBreederLBSnoCG, FUN = function(x){x[1]}))) /
    (GenerationTime)
  PM_Rmv_noCG_tibia_year <- round(mean(Rmv_noCG_tibia_year, adjust = 1), decimal)
  HPD_Rmv_noCG_tibia_year <- round(HPDinterval(Rmv_noCG_tibia_year), decimal)
  Prob_Rmv_noCG_tibia_pos_year <- mean(Rmv_noCG_tibia_year<0)
  
  # --- masse
  Rmv_noCG_masse_year <- as.mcmc(unlist(lapply(multivarBreederLBSnoCG, FUN = function(x){x[2]}))) /
    (GenerationTime)
  PM_Rmv_noCG_masse_year <- round(mean(Rmv_noCG_masse_year, adjust = 1), decimal)
  HPD_Rmv_noCG_masse_year <- round(HPDinterval(Rmv_noCG_masse_year), decimal)
  Prob_Rmv_noCG_masse_pos_year <- mean(Rmv_noCG_masse_year<0)
  
  
  
  # PER GENERATION ==========
  Rmv_noCG_tibia_generation <- as.mcmc(unlist(lapply(multivarBreederLBSnoCG, FUN = function(x){x[1]})))
  PM_Rmv_noCG_tibia_generation <- round(mean(Rmv_noCG_tibia_generation, adjust = 1), decimal)
  HPD_Rmv_noCG_tibia_generation <- round(HPDinterval(Rmv_noCG_tibia_generation), decimal)
  Prob_Rmv_noCG_tibia_pos_generation <- mean(Rmv_noCG_tibia_generation<0)
  
  # --- masse
  Rmv_noCG_masse_generation <- as.mcmc(unlist(lapply(multivarBreederLBSnoCG, FUN = function(x){x[2]})))
  PM_Rmv_noCG_masse_generation <- round(mean(Rmv_noCG_masse_generation, adjust = 1), decimal)
  HPD_Rmv_noCG_masse_generation <- round(HPDinterval(Rmv_noCG_masse_generation), decimal)
  Prob_Rmv_noCG_masse_pos_generation <- mean(Rmv_noCG_masse_generation<0)
  
  
  
  # OVERALL STUDY ===============
  # --- TIBIA
  Rmv_noCG_tibia_overall <- as.mcmc(unlist(lapply(multivarBreederLBSnoCG, FUN = function(x){x[1]}))) /
    (GenerationTime)
  PM_Rmv_noCG_tibia_overall <- round(mean(study_period_tibia*Rmv_noCG_tibia_overall, adjust = 1), decimal)
  HPD_Rmv_noCG_tibia_overall <- round(HPDinterval(study_period_tibia*Rmv_noCG_tibia_overall), decimal)
  Prob_Rmv_noCG_tibia_pos_overall <- mean(Rmv_noCG_tibia_overall<0)
  
  # --- masse
  Rmv_noCG_masse_overall <- as.mcmc(unlist(lapply(multivarBreederLBSnoCG, FUN = function(x){x[2]}))) /
    (GenerationTime)
  PM_Rmv_noCG_masse_overall <- round(mean(study_period_masse*Rmv_noCG_masse_overall, adjust = 1), decimal)
  HPD_Rmv_noCG_masse_overall <- round(HPDinterval(study_period_masse*Rmv_noCG_masse_overall), decimal)
  Prob_Rmv_noCG_masse_pos_overall <- mean(Rmv_noCG_masse_overall<0)
  
  
  ################################################
  # CONTRIBUTION OF GENETIC CORRELATION TO RESPONSE
  ################################################
  # --- tibia  
  cor_tibia <- cor.test(Rmv_tibia_year, Rmv_noCG_tibia_year)
  cor_tibia_estimate <- round(cor_tibia$estimate, decimal)
  cor_tibia_ic <- round(cor_tibia$conf.int, decimal)
  # 
  # contrib_Rmv_tibia_year <- (Rmv_tibia_year -  Rmv_noCG_tibia_year)^2 /   (Rmv_noCG_tibia_year + Rmv_tibia_year)^2
  # PM_contrib_rmv_tibia_year <- round(mean(contrib_Rmv_tibia_year), decimal)
  # HPD_contrib_rmv_tibia_year <- round(HPDinterval(contrib_Rmv_tibia_year), decimal)
  
  # --- masse
  cor_masse <- cor.test(Rmv_masse_year, Rmv_noCG_masse_year)
  cor_masse_estimate <- round(cor_masse$estimate, decimal)
  cor_masse_ic <- round(cor_masse$conf.int, decimal)
  # contrib_Rmv_masse_year <- (Rmv_masse_year -  Rmv_noCG_masse_year)^2  /   (Rmv_masse_year)^2
  # PM_contrib_rmv_masse_year <- round(mean(contrib_Rmv_masse_year), decimal)
  # HPD_contrib_rmv_masse_year <- round(HPDinterval(contrib_Rmv_masse_year), decimal)
  
  
  # STORE RESULT 
  output <- as.data.frame(matrix(NA, ncol = 2, nrow = 30, dimnames = list(
    c("selection_differential_year",
      "selection_differential_generation",
      "selection_differential_norm_year",
      "selection_differential_norm_generation",
      "p",
      "beta_year",
      "beta_generation",
      "beta_norm_year",
      "beta_norm_generation",
      "p",
      "R_MVB_CG_year",
      "R_MVB_CG_generation",
      "R_MVB_CG_overall",
      "p",
      "R_MVB_NO_CG_year",
      "R_MVB_NO_CG_generation",
      "R_MVB_NO_CG_overall",
      "p",
      "cor_CG_noCG",
      "va",
      "h2",
      "STS",
      "p",
      "cov_mother_w",
      "p",
      "cov_residuals_w",
      "p",
      "var_additive_fitness",
      "cor_pheno",
      "cor_gen"),
    c("masse", "tibia")
  )))
  
  # masse
  output[1,1] <- paste0(PM_S_masse_year, " [", HPD_S_masse_year[1], "; ", HPD_S_masse_year[2], "]")
  output[2,1] <- paste0(PM_S_masse_generation, " [", HPD_S_masse_generation[1], "; ", HPD_S_masse_generation[2], "]")
  output[3,1] <- paste0(PM_S_norm_masse_year, " [", HPD_S_norm_masse_year[1], "; ", HPD_S_norm_masse_year[2], "]")
  output[4,1] <- paste0(PM_S_norm_masse_generation, " [", HPD_S_norm_masse_generation[1], "; ", HPD_S_norm_masse_generation[2], "]")
  output[5,1] <- Prob_S_masse_pos_generation
  
  output[6,1] <- paste0(PM_beta_masse_year, " [", HPD_beta_masse_year[1], "; ", HPD_beta_masse_year[2], "]")
  output[7,1] <- paste0(PM_beta_masse_generation, " [", HPD_beta_masse_generation[1], "; ", HPD_beta_masse_generation[2], "]")
  
  output[8,1] <- paste0(PM_beta_sd_masse_year, " [", HPD_beta_sd_masse_year[1], "; ", HPD_beta_sd_masse_year[2], "]")
  output[9,1] <- paste0(PM_beta_sd_masse_generation, " [", HPD_beta_sd_masse_generation[1], "; ", HPD_beta_sd_masse_generation[2], "]")
  output[10,1] <- Prob_beta_sd_masse_pos_generation
  
  output[11,1] <- paste0(PM_rmv_masse_year, " [", HPD_rmv_masse_year[1], "; ", HPD_rmv_masse_year[2], "]")
  output[12,1] <- paste0(PM_rmv_masse_generation, " [", HPD_rmv_masse_generation[1], "; ", HPD_rmv_masse_generation[2], "]")
  output[13,1] <- paste0(PM_rmv_masse_overall, " [", HPD_rmv_masse_overall[1], "; ", HPD_rmv_masse_overall[2], "]")
  output[14,1] <- Prob_rmv_masse_pos_overall
  
  output[15,1] <- paste0(PM_Rmv_noCG_masse_year, " [", HPD_Rmv_noCG_masse_year[1], "; ", HPD_Rmv_noCG_masse_year[2], "]")
  output[16,1] <- paste0(PM_Rmv_noCG_masse_generation, " [", HPD_Rmv_noCG_masse_generation[1], "; ", HPD_Rmv_noCG_masse_generation[2], "]")
  output[17,1] <- paste0(PM_Rmv_noCG_masse_overall, " [", HPD_Rmv_noCG_masse_overall[1], "; ", HPD_Rmv_noCG_masse_overall[2], "]")
  output[18,1] <- Prob_Rmv_noCG_masse_pos_overall
  
  output[19,1] <- paste0(cor_masse_estimate, " [", cor_masse_ic[1], "; ", cor_masse_ic[2], "]")
  
  output[20,1] <- paste0(PM_VA_masse, " [", HPD_VA_masse[1], "; ", HPD_VA_masse[2], "]")
  output[21,1] <- paste0(PM_H_masse, " [", HPD_H_masse[1], "; ", HPD_H_masse[2], "]")
  
  output[22,1] <- paste0(PM_STS_masse, " [", HPD_STS_masse[1], "; ", HPD_STS_masse[2], "]")
  output[23,1] <- Prob_STS_masse_pos

  output[24,1] <- paste0(PM_cov_mother_masse, " [", HPD_cov_mother_masse[1], "; ", HPD_cov_mother_masse[2], "]")
  output[25,1] <- Prob_cov_mother_masse_pos

  output[26,1] <- paste0(PM_cov_residuals_masse, " [", HPD_cov_residuals_masse[1], "; ", HPD_cov_residuals_masse[2], "]")
  output[27,1] <- Prob_cov_residuals_masse_pos

  output[28,1] <- paste0(PM_var_additive_fitness, " [", HPD_var_additive_fitness[1], "; ", HPD_var_additive_fitness[2], "]")
  
  output[29,1] <- paste0(PM_cor_pheno, " [", HPD_cor_pheno[1], "; ", HPD_cor_pheno[2], "]")
  output[30,1] <- paste0(PM_cor_gen, " [", HPD_cor_gen[1], "; ", HPD_cor_gen[2], "]")
  
  
  # tibia
  output[1,2] <- paste0(PM_S_tibia_year, " [", HPD_S_tibia_year[1], "; ", HPD_S_tibia_year[2], "]")
  output[2,2] <- paste0(PM_S_tibia_generation, " [", HPD_S_tibia_generation[1], "; ", HPD_S_tibia_generation[2], "]")
  output[3,2] <- paste0(PM_S_norm_tibia_year, " [", HPD_S_norm_tibia_year[1], "; ", HPD_S_norm_tibia_year[2], "]")
  output[4,2] <- paste0(PM_S_norm_tibia_generation, " [", HPD_S_norm_tibia_generation[1], "; ", HPD_S_norm_tibia_generation[2], "]")
  output[5,2] <- Prob_S_tibia_pos_generation
  output[6,2] <- paste0(PM_beta_tibia_year, " [", HPD_beta_tibia_year[1], "; ", HPD_beta_tibia_year[2], "]")
  output[7,2] <- paste0(PM_beta_tibia_generation, " [", HPD_beta_tibia_generation[1], "; ", HPD_beta_tibia_generation[2], "]")
  output[8,2] <- paste0(PM_beta_sd_tibia_year, " [", HPD_beta_sd_tibia_year[1], "; ", HPD_beta_sd_tibia_year[2], "]")
  output[9,2] <- paste0(PM_beta_sd_tibia_generation, " [", HPD_beta_sd_tibia_generation[1], "; ", HPD_beta_sd_tibia_generation[2], "]")
  output[10,2] <- Prob_beta_sd_tibia_pos_generation
  output[11,2] <- paste0(PM_rmv_tibia_year, " [", HPD_rmv_tibia_year[1], "; ", HPD_rmv_tibia_year[2], "]")
  output[12,2] <- paste0(PM_rmv_tibia_generation, " [", HPD_rmv_tibia_generation[1], "; ", HPD_rmv_tibia_generation[2], "]")
  output[13,2] <- paste0(PM_rmv_tibia_overall, " [", HPD_rmv_tibia_overall[1], "; ", HPD_rmv_tibia_overall[2], "]")
  output[14,2] <- Prob_rmv_tibia_pos_overall
  output[15,2] <- paste0(PM_Rmv_noCG_tibia_year, " [", HPD_Rmv_noCG_tibia_year[1], "; ", HPD_Rmv_noCG_tibia_year[2], "]")
  output[16,2] <- paste0(PM_Rmv_noCG_tibia_generation, " [", HPD_Rmv_noCG_tibia_generation[1], "; ", HPD_Rmv_noCG_tibia_generation[2], "]")
  output[17,2] <- paste0(PM_Rmv_noCG_tibia_overall, " [", HPD_Rmv_noCG_tibia_overall[1], "; ", HPD_Rmv_noCG_tibia_overall[2], "]")
  output[18,2] <- Prob_Rmv_noCG_tibia_pos_overall
  output[19,2] <- paste0(cor_tibia_estimate, " [", cor_tibia_ic[1], "; ", cor_tibia_ic[2], "]")
  output[20,2] <- paste0(PM_VA_tibia, " [", HPD_VA_tibia[1], "; ", HPD_VA_tibia[2], "]")
  output[21,2] <- paste0(PM_H_tibia, " [", HPD_H_tibia[1], "; ", HPD_H_tibia[2], "]")
  output[22,2] <- paste0(PM_STS_tibia, " [", HPD_STS_tibia[1], "; ", HPD_STS_tibia[2], "]")
  output[23,2] <- Prob_STS_tibia_pos
  output[24,2] <- paste0(PM_cov_mother_tibia, " [", HPD_cov_mother_tibia[1], "; ", HPD_cov_mother_tibia[2], "]")
  output[25,2] <- Prob_cov_mother_tibia_pos
  output[26,2] <- paste0(PM_cov_residuals_tibia, " [", HPD_cov_residuals_tibia[1], "; ", HPD_cov_residuals_tibia[2], "]")
  output[27,2] <- Prob_cov_residuals_tibia_pos
  output[28,2] <- paste0(PM_var_additive_fitness, " [", HPD_var_additive_fitness[1], "; ", HPD_var_additive_fitness[2], "]")
  output[29,2] <- paste0(PM_cor_pheno, " [", HPD_cor_pheno[1], "; ", HPD_cor_pheno[2], "]")
  output[30,2] <- paste0(PM_cor_gen, " [", HPD_cor_gen[1], "; ", HPD_cor_gen[2], "]")
  
  
  return(list(beta_sd = cbind(masse = grad_sd_masse_generation, tibia = grad_sd_tibia_generation),
              s_differential_norm = cbind(masse = S_norm_masse_generation, tibia = S_norm_tibia_generation),
              MBE_CG = cbind(masse = Rmv_masse_generation, tibia = Rmv_tibia_generation),
              MBE_NO_CG = cbind(masse = Rmv_noCG_masse_generation, tibia = Rmv_noCG_tibia_generation),
              STS = cbind(masse = STS_masse, tibia = STS_tibia),  # by generation
              cov_m = cbind(masse = cov_mother_masse, tibia = cov_mother_tibia),  # by generation
              cov_residuals = cbind(masse = cov_residuals_masse, tibia = cov_residuals_tibia),  # by generation
              G = Gmode,
              S = Smode,
              P = Pmode,
              summary = output))
}







# VIOLIN PLOT -------------------------------------------------------
#' @param mod_result - MCMCglmm model result
#' @param trait - trait from which calculate parameters
#' @param type - should result be visualized as ratio or in orginal variance scale

#'@WARNINGS : trait have to be entered in the same way as in the model (with LRS_relative at last)

#' @return a dataframe containing genetic parameters statistics : median and IC 

compute_stats_re_bivariate_violin_plot <- function(mod_result = mod_result,
                                                   type = c("variance", "ratio"),
                                                   trait = c("tibia", "masse"),
                                                   ylim = c(0, 0.5)
)
{
  
  # GET RE NAMES
  re <- colnames(mod_result$VCV)
  
  # IS FITNESS IN MODEL ?
  fitness_bool <- any(grepl(pattern = "LRS", x = re))
  
  # EXTRACT RELEVANT VARIABLE FOR DENOMINATOR TRAIT IN GENETIC PARAMETERS COMPUTATION
  for(trait.temp in trait){
    trait_re.temp <- re[grep(pattern = paste0(paste0("trait", trait.temp, ":trait", trait.temp), collapse = "|"), x = re)]
    assign(paste0(trait.temp, "_re"), trait_re.temp)
  }
  
  # Get cohort variance effect [syntax depend on the number of trait and the presence of fitness as cohort effect is not dependent on cohort]
  if(fitness_bool == T){
    # If more than two non-fitness trait included in model, set covariance name with cohort (trait[i], trait[j])
    if(length(trait)>1){
      for(i in  1:length(trait)){
        add.cohort.re.temp <- c(get(paste0(trait[i], "_re")),
                                paste0("at.level(trait, c(\"",trait[1], "\", \"", trait[2], "\"))", which(trait == trait[i]), ":at.level(trait, c(\"", trait[1], "\", \"", trait[2], "\"))", which(trait == trait[i]),".annee_emergence"))
        
        assign(paste0(trait[i], "_re"), add.cohort.re.temp)
      }
      
      #if only one non-fitness trait, set covariance name with cohort as (trait[i], trait[i])
    } else if(length(trait) == 1){
      add.cohort.re.temp <- c(get(paste0(trait, "_re")),
                              paste0("at.level(trait, c(\"", trait, "\")):at.level(trait, c(\"", trait, "\")).annee_emergence"))
      assign(paste0(trait, "_re"), add.cohort.re.temp)
    }
    
  } 
  
  
  
  # DATAFRAME TO STORE RESULTS
  df <- as.data.frame(matrix(ncol = length(get(paste0(trait[1], "_re"))) * length(trait),  # for y / ymin / ymax
                             nrow = nrow(mod_result$VCV))
  )
  
  # EXTRACT DYNAMICALLY NAME FROM COLNAMES $VCV
  variance_parameters <- gsub(pattern = paste0("trait", trait[1], ":trait", trait[1], ".|trait", trait[1], "."),
                              replacement = "",
                              x = get(paste0(trait[1], "_re")))
  
  # --- Rename cohort effect as "cohort"
  cohort_name <- variance_parameters[grepl(pattern = "emergence", x = variance_parameters)]
  variance_parameters[which(variance_parameters == cohort_name)] <- "cohort"
  
  # --- Add as colnames
  col.names.temp <- c()
  for(trait.temp in trait){
    col.names.temp <- c(col.names.temp, paste0(variance_parameters, '_', trait.temp))
  }
  
  colnames(df) <- col.names.temp  # change colnames for parameters values summary
  
  
  #} 
  
  
  
  # COMPUTE MEAN AND IC PER TRAIT FOR EACH RANDOM EFFECT
  # GET VARIANCE ASSOCIATED WITH FIXED EFFECT (De Villemereuil, 2018)
  compute_vcvpred <- function(beta, design_matrix, ntraits) { 
    
    list(cov(matrix(design_matrix %*% beta, ncol = ntraits))) 
    
  }
  
  X <- mod_result$X  # desgin matrix
  X_mat <- as.matrix(X)  # convert to matrix
  col.names.temp <- colnames(mod_result$Sol)  # colnames of prediction
  # get index associating with last fixed effect
  last_row <- grep(pattern = "animal", col.names.temp)[1] - 1
  last_fixed_effect_index <- ifelse(!is.na(last_row), last_row, length(col.names.temp)) # needs to adjust because fixed effects change according to period modeled
  sol <- mod_result$Sol[,1:last_fixed_effect_index]  # select only fixed effects
  fixed_effects_var <- flatten( apply(sol, 1, compute_vcvpred, design_matrix = X, ntraits = length(trait)) )  # get variance-covariance matrix associated with prediction
  
  
  # extract variance associated with each trait
  for(trait.temp in trait){
    vf.temp <- unlist(lapply(X = fixed_effects_var, FUN = function(x) x[which(trait.temp == trait), which(trait.temp == trait)]))
    assign(paste0("vf_", trait.temp), vf.temp)
  }  
  
  
 
  
  if(type == "ratio"){
    col.index <- 0  # to ajdust column filling according to the number of trait
    for(trait.temp in trait){  # for each trait
      for(i in 1:length(get(paste0(trait[1], "_re")))){
        
        df[, i + col.index] <- mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]] / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))
        
      }
      col.index <- length(get(paste0(trait[1], "_re")))  # to ajdust column filling according to the number of trait and random effects
      
      # Add Vf
      df[paste0("vf_", trait.temp)] <- get(paste0("vf_", trait.temp))  / (rowSums(mod_result$VCV[,get(paste0(trait.temp, "_re"))]) + get(paste0("vf_", trait.temp)))
      
      }
    

    
  } else if (type == "variance"){
    col.index <- 0  # to ajdust column filling according to the number of trait
    for(trait.temp in trait){  # for each trait
      for(i in 1:length(get(paste0(trait[1], "_re")))){
        
        df[, i + col.index] <- mod_result$VCV[,get(paste0(trait.temp, "_re"))[i]]
        
      }
      col.index <- length(get(paste0(trait[1], "_re")))  # to adjust column filling according to the number of trait and random effects
     
       # Add Vf
      df[paste0("vf_", trait.temp)] <- get(paste0("vf_", trait.temp))
      }
  }
 
  
  
  
  # REFORMAT TABLE FOR GGPLOT
  variance_parameters <- c(variance_parameters, "vf")  # add vf to variance parameters
  params <- rep(rep(variance_parameters, each = nrow(mod_result$VCV)), length(trait))
  df_gg.temp <- pivot_longer(data = df, cols = 1:ncol(df), names_to = "params", values_to = "y")
  df_gg <- separate(data = df_gg.temp, col = params, into = c("params", "trait"), sep = "_", remove = FALSE)  # create two columns "params" and "trait" by separating strings
  df_gg$trait <- factor(df_gg$trait, levels = c("masse", "tibia"))
  df_gg$params <- factor(df_gg$params, levels = c("animal", "mother", "cohort", "vf", "units"))
  
  
  # --- x axis labels
  # labels_x_axis <- 
  #   c(TeX("$\\sigma^2 _{a}$"),
  #     TeX("$\\sigma^2 _{m}$"),
  #     TeX("$\\sigma^2 _{c}$"),
  #     TeX("$\\sigma^2 _{f}$"),
  #     TeX("$\\sigma^2 _{r}$")
  #   )
  
  labels_x_axis <- 
    c(TeX("a", bold = T),
      TeX("m", bold = T),
      TeX("c", bold = T),
      TeX("f", bold = T),
      TeX("r", bold = T)
    )
  
  # Add unity if y in original scale
  if(type == "variance"){
    df_gg$trait <- factor(df_gg$trait, levels = c("masse", "tibia"), labels = c(TeX("Mass $(g^2)$", bold = T), TeX("Tibia length $(mm^2)$", bold = T)))
  } else{
    df_gg$trait <- factor(df_gg$trait, levels = c("masse", "tibia"), labels = c(TeX("Mass", bold = T), TeX("Tibia length$", bold = T)))
  }
  
  for(trait.temp in trait){
    # --- ylab ~ type
    if(type == "ratio"){
      ylab = TeX("$\\sigma^2_i / \\sigma^2_P$")
    } else if (type == "variance"){
      ylab = TeX("$\\sigma^2$", bold = T)
    }
    
    
    plot.temp1 <- ggplot(df_gg[df_gg$trait ==  "bold(`Tibia length `(mm^{\n    2\n}))",], aes(x = params, y = y, fill = params)) +
      geom_violin(trim = FALSE) +
      geom_boxplot(width = 0.1, size = 0.7, fill = "white", outliers = F) +
      facet_wrap(.~trait, scale = "free", labeller = label_parsed, nrow = 1, ncol = 2) +
      labs(title = NULL, x = NULL, y = ylab) +
      tidyquant::theme_tq() +
      scale_x_discrete(
        labels = labels_x_axis) +  # define above
      scale_fill_brewer(palette="Greys", direction = -1) +
      
      xlab("Variance components") +
      
      # ggh4x::facetted_pos_scales(y = list(trait == "tibia" ~ scale_y_continuous(limits = c(0, tibia_y_scale),
      #                                                                           breaks = seq(0, tibia_y_scale, by = 5)
      #                                                                           
      #                                                                           ),
      #                                     trait == "masse" ~ scale_y_continuous(limits = c(0, max(df_gg$y[df_gg$trait == "masse"])*1.1),
      #                                                                            breaks = seq(0, tibia_y_scale, by = 5)
      #                                                                           
      #                                                                           )
    #                                     )
    #                            ) +
    
    # scale_y_continuous(
    #   breaks = seq(0, max(y), length.out = 4),
    #   labels = seq(0, max(y), length.out = 4),
    #   limits = c(0, max(y))
    # ) +
    theme(text = element_text(size = 20),
          legend.position = "none",
          axis.text.x = element_text(size = 23),
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(size = 25, margin = margin(0, 9, 0, 0)),  # space between ylab and plot
          axis.title.x = element_text(size = 20, face = "bold"),  # space between ylab and plot
          
          strip.background.x = element_rect(fill = "#4F4F4F", color = "#4F4F4F", linewidth = 0.5),  # fill of panels x
          strip.background.y = element_rect(fill = "#4F4F4F", color = "#4F4F4F", linewidth = 0.5),  # fill of panels x
          strip.text = element_text(color = "white", size = 17, face = "bold"),  # text color of panels
          strip.text.x = element_text(margin = margin (3, 3, 5, 3)),  # space around x panels text (t, r, b, l)
          strip.text.y = element_text(margin = margin (3, 4, 3, 4)),  # space around y panels text
          axis.ticks.y = element_blank(),
    ) 
    
    if(!is.null(ylim)){
      plot.temp <- plot.temp +
        ylim(ylim) 
    }
  }
  
  
  gridExtra::grid.arrange(plot.temp, plot.temp1, nrow = 2)
  
}







# PLOT RESPONSE TO SELECTION ---------------------------------------------
plot_response_selection <- function(output_pooled = output_pooled,
                                    output_P1 = output_P1,
                                    output_P2 = output_P2,
                                    bv_data_drift = drift$bv,
                                    color_period = rev(c("#36465d", "#cc3300", "#4e9785"))){
  
  # reformatin col "ebv" of data drift to match those of other df
  names(bv_data_drift)[which(colnames(bv_data_drift) == "ebv")] <- "response"
  
  # Get EBV per generation
  bv_data_drift$response <- bv_data_drift$response
  
  # formating df to store response to selection
  # --- pooled
  df.temp1 <- pivot_longer(data = as.data.frame(output_pooled$MBE_NO_CG), cols = 1:2, names_to = "trait", values_to = "response")  # MBE without genetic correlation
  df.temp2 <- pivot_longer(data = as.data.frame(output_pooled$MBE_CG), cols = 1:2, names_to = "trait", values_to = "response")  # MBE
  df.temp3 <- pivot_longer(data = as.data.frame(output_pooled$STS), cols = 1:2, names_to = "trait", values_to = "response")  # STS
  df.temp4 <- bv_data_drift[bv_data_drift$period == "Pooled", c("trait", "response")]   # Breeding values
  
  df_pooled <- rbind(df.temp1, df.temp2, df.temp3, df.temp4)
  df_pooled["equation"] <- rep(c("MBE_NO_CG", "MBE_CG", "STS", "BV"), each = nrow(output_pooled$MBE_CG)*2)  # *2 == for each trait
  df_pooled["period"] <- "Overall study"
  
  # --- P1
  df.temp1 <- pivot_longer(data = as.data.frame(output_P1$MBE_NO_CG), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp2 <- pivot_longer(data = as.data.frame(output_P1$MBE_CG), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp3 <- pivot_longer(data = as.data.frame(output_P1$STS), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp4 <- bv_data_drift[bv_data_drift$period == "Period 1", c("trait", "response")]   # Breeding values
  
  df_P1 <- rbind(df.temp1, df.temp2, df.temp3, df.temp4)
  df_P1["equation"] <- rep(c("MBE_NO_CG", "MBE_CG", "STS", "BV"), each = nrow(output_P1$MBE_CG)*2)
  df_P1["period"] <- "1997-2005"
  
  
  # --- P2
  df.temp1 <- pivot_longer(data = as.data.frame(output_P2$MBE_NO_CG), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp2 <- pivot_longer(data = as.data.frame(output_P2$MBE_CG), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp3 <- pivot_longer(data = as.data.frame(output_P2$STS), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp4 <- bv_data_drift[bv_data_drift$period == "Period 2", c("trait", "response")]   # Breeding values
  
  df_P2 <- rbind(df.temp1, df.temp2, df.temp3, df.temp4)
  df_P2["equation"] <- rep(c("MBE_NO_CG", "MBE_CG", "STS", "BV"), each = nrow(output_P2$MBE_CG)*2)
  df_P2["period"] <- "2006-2023"
  
  
  # Merge dataframe for visualisation
  df <- rbind(df_pooled, df_P1, df_P2)
  
  
  # Reorder factor
  df$period <- factor(df$period, levels = rev(c("Overall study", "1997-2005", "2006-2023")))
  df$trait <- factor(df$trait, levels = c("masse", "tibia"), labels = c(TeX("\\textbf{Mass (g)}"), TeX("\\textbf{Tibia length (mm)}")))
  df$equation <- factor(df$equation, levels = rev(c("BV", "MBE_NO_CG", "MBE_CG", "STS")),
                        labels = rev(c(TeX("PBV's slopes", bold = T),
                                       TeX("$MBE_{\\rho=0}$", bold = T),
                                       TeX("MBE", bold = T),
                                       TeX("STS", bold = T))))
  
  
  # VISUALISATION  
  ggplot(data = df, aes(x = period, y = response, fill = period)) +
    stat_halfeye(
      mapping = aes(fill_ramp = after_stat(y > 0)),  
      justification = -.2,  # Move geom to the right to have space for boxplot
      .width = 0,   # remove default boxplot
      point_colour = NA,  # remove default boxplot
      slab_color = "black",  # line around density curves
      slab_size = 0.6,  # thickness of line around density curves
      alpha = 0.8,
      expand = F,  # expand lines below plot to min and max
      adjust = 0.8,  # smothness of curves (max to 1)
      normalize = "panels"  # density of each plot varying between 0 and 1
    ) +
    
    # Colour for the different periodes
    scale_fill_manual(values = c(color_period[3], color_period[2], color_period[1]), 
                      limits = rev(levels(df$period)), 
                      guide = guide_legend(title = "period")) +
   
     # Colour per condition : bv < 0 or bv > 0
    scale_fill_ramp_discrete(na.translate = F,
                             guide = "none") +  # remove it from legend
    # Pic at 0
    # stat_spike(at = 0, linewidth = 0.5) +
    
    geom_boxplot(width = .12,
                 linewidth = 0.5,  # linewidth of boxplot
                 outlier.color = NA,  # if NA, do not show outliers
                 alpha = 0.7, show.legend = F) +
    coord_flip() +
    # Facets and parameters
    facet_grid(equation ~ trait,  
               scales = "free",  # allow each facet to have its own y and x axis
               labeller = label_parsed  # name for facets panels
    ) +
    
    
    geom_hline(yintercept = 0, size = 0.5, show.legend = NA) +  # add black line of 0
    
    # ylab per facet
    ggh4x::facetted_pos_scales(y = list(trait == "tibia" ~ scale_y_continuous(breaks = seq(-2, 3, 1),
                                                                              labels = seq(-2, 3, 1),
                                                                              limits = c(-2, 3)
                                                                              ),
                                        
                                        trait == "masse" ~ scale_y_continuous(breaks = seq(-40, 60, 20),
                                                                              labels = seq(-40, 60, 20),
                                                                              limits = c(-41, 61)
                                                                              )
                                        )
                               ) +

    # labs(y = latex2exp::TeX("$\\Delta (BV_{predicted} \\ , \\ BV_{drift})$"),
    #     x = "Density") +
    labs(y = "Predicted response per generation",
         x = "Probability density") +
    
    # OVERALL THEME
    tidyquant::theme_tq() +
    
    # Annotate p-values
    
    # THEME
    theme(text = element_text(size = 20),
          legend.position = "none",
          # legend.title = element_blank(),
          # legend.text = element_text(size = 13),
          legend.box.background = element_rect(color = "black", size = 1),
          strip.background.x = element_rect(fill = "#4F4F4F", color = "#4F4F4F", linewidth = 0.5),  # fill of panels x
          strip.background.y = element_rect(fill = "#4F4F4F", color = "#4F4F4F", linewidth = 0.5),  # fill of panels x
          strip.text = element_text(color = "white", size = 15, face = "bold"),  # text color of panels
          strip.text.x = element_text(margin = margin (3, 3, 5, 3), face = "bold"),  # space around x panels text (t, r, b, l)
          strip.text.y = element_text(margin = margin (3, 4, 3, 4), face = "bold"),  # space around y panels text
          axis.ticks.y = element_blank(),
          axis.text.y=element_text(size = 15),
          # axis.text.x = element_text(margin = margin(t = 5)), # Adjust padding at the bottom of x-axis labels
          axis.title.x = element_text(margin = margin(t = 10), face = "bold"), # Adjust padding at the bottom of x-axis title
          axis.title.y = element_text(size = 22, margin = margin(0, 9, 0, 0), face = "bold"),  # space between ylab and plot
          
    ) 
    
   # guides(fill = guide_legend(override.aes = list(linetype = 0, shape = 0, color = "black")))  # keep only fill as legend color and style
}










# PLOT GRADIENT AND SELECTION ---------------------------------------------
plot_selection_gradient <- function(output_pooled = output_pooled,
                                    output_P1 = output_P1,
                                    output_P2 = output_P2,
                                    color_period = c("#36465d", "#cc3300", "#4e9785")){
  
  # formating df to store selection parameters
  # --- pooled
  df.temp1 <- pivot_longer(data = as.data.frame(output_pooled$s_differential_norm), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp2 <- pivot_longer(data = as.data.frame(output_pooled$beta_sd), cols = 1:2, names_to = "trait", values_to = "response")
  
  df_pooled <- rbind(df.temp1, df.temp2)
  df_pooled["equation"] <- rep(c("s", "beta"), each = nrow(output_pooled$MBE_CG)*2)
  df_pooled["period"] <- "Overall study"
  
  # --- P1
  df.temp1 <- pivot_longer(data = as.data.frame(output_P1$s_differential_norm), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp2 <- pivot_longer(data = as.data.frame(output_P1$beta_sd), cols = 1:2, names_to = "trait", values_to = "response")
  
  df_P1 <- rbind(df.temp1, df.temp2)
  df_P1["equation"] <- rep(c("s", "beta"), each = nrow(output_P1$MBE_CG)*2)
  df_P1["period"] <- "1997-2005"
  
  
  # --- P2
  df.temp1 <- pivot_longer(data = as.data.frame(output_P2$s_differential_norm), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp2 <- pivot_longer(data = as.data.frame(output_P2$beta_sd), cols = 1:2, names_to = "trait", values_to = "response")
  
  df_P2 <- rbind(df.temp1, df.temp2)
  df_P2["equation"] <- rep(c("s", "beta"), each = nrow(output_P2$MBE_CG)*2)
  df_P2["period"] <- "2006-2023"
  
  
  # Merge dataframe for visualisation
  df <- rbind(df_pooled, df_P1, df_P2)
  
  
  # Reorder factor
  df$period <- factor(df$period, levels = rev(c("Overall study", "1997-2005", "2006-2023")))
  df$trait <- factor(df$trait, levels = c("masse", "tibia"))
  df$equation <- factor(df$equation, levels = c("s", "beta"))
  
  
  # VISUALISATION  
  ggplot(data = df, aes(x = period, y = response, fill = period)) +
    stat_halfeye(
      mapping = aes(fill_ramp = after_stat(y > 0)), 
      justification = -.2,  # Move geom to the right to have space for boxplot
      .width = 0,   # remove default boxplot
      point_colour = NA,  # remove default boxplot
      slab_color = "black",  # line around density curves
      slab_size = 0.6,  # thickness of line around density curves
      alpha = 0.8,
      expand = F,  # expand lines below plot to min and max
      adjust = 0.8,  # smothness of curves (max to 1)
      normalize = "panels"  # density of each plot varying between 0 and 1
    ) +
    
    # Colour for the different periodes
    scale_fill_manual(values = c(color_period[3], color_period[2], color_period[1])) +
    # Colour per condition : bv < 0 or bv > 0
    scale_fill_ramp_discrete(na.translate = F) +  
    # Pic at 0
    # stat_spike(at = 0, linewidth = 0.5) +
    
    geom_boxplot(width = .12,
                 linewidth = 0.5,  # linewidth of boxplot
                 outlier.color = NA,  # if NA, do not show outliers
                 alpha = 0.7) +
    coord_flip() +
    # Facets and parameters
    facet_grid(equation ~ trait,  
               scales = "free",  # allow each facet to have its own y and x axis
               labeller = labeller(trait = c("masse" = "Mass", "tibia" = "Tibia length"),
                                   equation = c(s = "S", beta = "Beta"))  # name for facets panels
    ) +
    
    
    geom_hline(yintercept = 0, size = 0.5) +  # add black line of 0
    
    # Set ylab for each facets with <ggh4x> according above the ggplot script
    # ggh4x::facetted_pos_scales(y = list(trait == "tibia" ~ scale_y_continuous(limits = lim_tibia),
    #                                     trait == "masse" ~ scale_y_continuous(limits = lim_masse))) +
    
    # labs(y = latex2exp::TeX("$\\Delta (BV_{predicted} \\ , \\ BV_{drift})$"),
    #     x = "Density") +
    labs(y = "Parameter values",
         x = "Probability density") +
    
    # OVERALL THEME
    tidyquant::theme_tq() +
    
    # Annotate p-values
    
    # THEME
    theme(text = element_text(size = 20),
          legend.position = "none",
          strip.background.x = element_rect(fill = "#4F4F4F", color = "#4F4F4F", linewidth = 0.5),  # fill of panels x
          strip.background.y = element_rect(fill = "#4F4F4F", color = "#4F4F4F", linewidth = 0.5),  # fill of panels x
          strip.text = element_text(color = "white", size = 18, face = "bold"),  # text color of panels
          strip.text.x = element_text(margin = margin (3, 3, 5, 3)),  # space around x panels text (t, r, b, l)
          strip.text.y = element_text(margin = margin (3, 4, 3, 4)),  # space around y panels text
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(margin = margin(t = 5)), # Adjust padding at the bottom of x-axis labels
          axis.title.x = element_text(margin = margin(t = 10)), # Adjust padding at the bottom of x-axis title
          axis.title.y = element_text(size = 22, margin = margin(0, 9, 0, 0)),  # space between ylab and plot
          
    ) 
}




# PLOT COVARIANCE FITNESS ---------------------------------------------
plot_cov_fitness <- function(output_pooled = output_pooled,
                                    output_P1 = output_P1,
                                    output_P2 = output_P2,
                                    color_period = c("#36465d", "#cc3300", "#4e9785")){
  
  # additive genetic covariance trait, fitness
  # cova_pooled <- additive_covariance_matrix$vcv_pooled
  # cova_pooled <- data.frame(masse = cova_pooled$`traitLRS_relative:traitmasse.animal`, tibia = cova_pooled$`traitLRS_relative:traittibia.animal`)
  # cova_P1 <- additive_covariance_matrix$vcv_P1
  # cova_P1 <- data.frame(masse = cova_P1$`traitLRS_relative:traitmasse.animal`, tibia = cova_P1$`traitLRS_relative:traittibia.animal`)
  # cova_P2 <- additive_covariance_matrix$vcv_P2
  # cova_P2 <- data.frame(masse = cova_P2$`traitLRS_relative:traitmasse.animal`, tibia = cova_P2$`traitLRS_relative:traittibia.animal`)
  
  # formating df to store response to selection
  # --- pooled
  df.temp1 <- pivot_longer(data = as.data.frame(output_pooled$STS), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp2 <- pivot_longer(data = as.data.frame(output_pooled$cov_m), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp3 <- pivot_longer(data = as.data.frame(output_pooled$cov_residuals), cols = 1:2, names_to = "trait", values_to = "response")
  
  df_pooled <- rbind(df.temp1, df.temp2, df.temp3)
  df_pooled["equation"] <- rep(c("sts", "cov_m", "cov_r"), each = nrow(output_pooled$STS)*2)
  df_pooled["period"] <- "Overall study"
  
  # --- P1
  df.temp1 <- pivot_longer(data = as.data.frame(output_P1$STS), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp2 <- pivot_longer(data = as.data.frame(output_P1$cov_m), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp3 <- pivot_longer(data = as.data.frame(output_P1$cov_residuals), cols = 1:2, names_to = "trait", values_to = "response")
  
  df_P1 <- rbind(df.temp1, df.temp2, df.temp3)
  df_P1["equation"] <- rep(c("sts", "cov_m", "cov_r"), each = nrow(output_P1$STS)*2)
  df_P1["period"] <- "1997-2005"
  
  
  # --- P2
  df.temp1 <- pivot_longer(data = as.data.frame(output_P2$STS), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp2 <- pivot_longer(data = as.data.frame(output_P2$cov_m), cols = 1:2, names_to = "trait", values_to = "response")
  df.temp3 <- pivot_longer(data = as.data.frame(output_P2$cov_residuals), cols = 1:2, names_to = "trait", values_to = "response")
  
  df_P2 <- rbind(df.temp1, df.temp2, df.temp3)
  df_P2["equation"] <- rep(c("sts", "cov_m", "cov_r"), each = nrow(output_P2$STS)*2)
  df_P2["period"] <- "2006-2023"
  
  
  # Merge dataframe for visualisation
  df <- rbind(df_pooled, df_P1, df_P2)
  
  
  # Reorder factor
  df$period <- factor(df$period, levels = c("Overall study", "1997-2005", "2006-2023"))
  df$trait <- factor(df$trait, levels = c("masse", "tibia"))
  df$equation <- factor(df$equation, levels = c("sts", "cov_m", "cov_r"),
                        labels = c("Genetic", "Mother", "Environment"))
  
  
  # VISUALISATION  
  ggplot(data = df, aes(x = equation, y = response, fill = period)) +
    geom_hline(yintercept = 0, size = 0.5, linetype = 2) +  # add black line of 0
    
    geom_boxplot(width = 0.5,
                 linewidth = 0.8,  # linewidth of boxplot
                 outlier.color = NA,  # if NA, do not show outliers
                 alpha = 0.8
    ) +
    
    # Colour for the different periodes
    scale_fill_manual(values = c(color_period[1], color_period[2], color_period[3])) +
    # Colour per condition : bv < 0 or bv > 0
    # Pic at 0
    # stat_spike(at = 0, linewidth = 0.5) +
    
    # Facets and parameters
    facet_grid(trait ~ .,  
               scales = "free", labeller = labeller(trait = c(masse = "Mass (g)", tibia = "Tibia length (mm)"))
    ) +
    
    

    
    #Set ylab for each facets with <ggh4x> according above the ggplot script
    ggh4x::facetted_pos_scales(y = list(trait == "tibia" ~ scale_y_continuous(limits = c(-1.5,3)),
                                        trait == "masse" ~ scale_y_continuous(limits = c(-30, 60)))) +
    
    # labs(y = latex2exp::TeX("$\\Delta (BV_{predicted} \\ , \\ BV_{drift})$"),
    #     x = "Density") +
    labs(y = "Selection differential",
         x = "Source of covariation with fitness") +
    
    # OVERALL THEME
    tidyquant::theme_tq() +
    
    # Annotate p-values
    
    # THEME
    theme(
          legend.position = "none",
          
          strip.background.y =  element_blank(),  # fill of panels x
          strip.text = element_text(color = "black", size = 22),  # text color of panels
          strip.text.y = element_text(margin = margin (t = 3, r = 4, b = 3, l = 9)),  # space around y panels text
          axis.ticks.y = element_blank(),
          axis.text.y = element_text( size = 18), # Adjust padding at the bottom of x-axis labels
          axis.text.x = element_text(margin = margin(t = 5), size = 18), # Adjust padding at the bottom of x-axis labels
          axis.title.x = element_text(margin = margin(t = 10), face = "bold"), # Adjust padding at the bottom of x-axis title
          axis.title.y = element_text(size = 22, margin = margin(0, 9, 0, 0), face = "bold"),  # space between ylab and plot
          text = element_text(size = 20)
          
    ) 
}



# PLOT BREEDING VALUES POSTERIOR  -------------------------------------------------------
#' @param mod_result - lme4 model result
#' @param data - data frame used for model 
#' @param trait - trait to compute
#' @param position - "median" or "mean", used to calculate bv 
#' @param periode - a list of equal length than "trait" containg the periode used to compute tmeporal variation of BV
#' @param method_regression - "lm" or "loess", regression method to compute temporal trend 
#' @param main - list containing title of upper plots

#' @return Return mean breeding values according to cohort 

plot_breeding_values_posterior <- function(mod_result,
                                           data = data,
                                           trait = c("masse", "tibia"))
{
  
  for(trait.temp in trait){
    # COMPUTE BV --------------------------------------------------------------
    # EXTRACTION BV
    sol <- mod_result$Sol
    
    # Column name of sol
    col.names.sol <- colnames(mod_result$Sol)
    
    # Filter sol for BV 
    sol <- sol[, grep(pattern = paste0("trait", trait.temp, ".animal"), col.names.sol) ]
    col.names.sol <- colnames(sol)
    
    # Extract id
    id <- gsub(pattern = paste0("trait", trait.temp, ".animal."), replacement = "", x =  col.names.sol)
    colnames(sol) <- id
    
    # Filter with individuals in df
    sol <- sol[, id %in% data$individual_id]
    
    
    # Associated cohort year
    index <- match(x = colnames(sol), table = pedigree_data$individual_id)
    cohort <- pedigree_data[index, "cohort"]
    
    # x for plot
    cohort_level <- 1997:2023
    
    # Overall mean
    mean_bv_posterior <- 
      t(apply(X = sol, MARGIN = 1, FUN = function(x){
        tapply(X = x, INDEX = cohort, FUN = mean)
      }))
    
    overall_mean <- apply(X = mean_bv_posterior, MARGIN = 2, FUN = mean)
    
    # Create model
    # mod_pooled <- lm(overall_mean ~ cohort_level)
    # mod_P1 <- lm(overall_mean[1:9] ~ cohort_level[1:9])
    # mod_P2 <- lm(overall_mean[10:length(cohort_level)] ~ cohort_level[10:length(cohort_level)])
    # 
    # # IC for prediction
    # ic_pooled <- as.data.frame(predict(mod_pooled, interval = "confidence", se.fit = T))
    # ic_pooled["annee_emergence"] <- rownames(ic_pooled)
    # ic_P1 <- as.data.frame(predict(mod_P1, interval = "confidence"))
    # ic_P2 <- as.data.frame(predict(mod_P2, interval = "confidence"))
    
    
    
    
    # MEAN BREEDING VALUES PER COHORT FOR EACH ITERATION MCMC
    mean_bv_iteration <-  apply(X = sol, MARGIN = 1, FUN = function(x){
      predict(loess(tapply(x, cohort, mean) ~ cohort_level, span = 0.7))
    })
    
    # --- reformating for gggplot
    mean_bv_iteration_raw <- t(as.data.frame(mean_bv_iteration))
    colnames(mean_bv_iteration_raw) <- cohort_level
    mean_bv_iteration <- as.data.frame(mean_bv_iteration_raw)
    mean_bv_iteration["iteration"] <- 1:1000  # add iteration to plot each sample
    mean_bv_iteration <- pivot_longer(data = mean_bv_iteration, cols = 1:27, values_to = "bv", names_to = "annee_emergence")
    
    # --- add period
    mean_bv_iteration["period"] <- ifelse(as.numeric(as.character(mean_bv_iteration$annee_emergence)) <= 2005, "period 1", "period 2")
    
    
    # IC INTERVAL OF BV
    ic_lwr <- apply(mean_bv_posterior, MARGIN = 2, FUN = function(x) quantile(x, probs = 0.025))
    ic_upr <- apply(mean_bv_posterior, MARGIN = 2, FUN = function(x) quantile(x, probs = 0.975))
    ic_df <- data.frame(ic_lwr, ic_upr, annee_emergence= cohort_level)
    
    
    
    # PLOT
    # --- Parametring labs
    if(trait.temp == "tibia"){
      ylab <- "Tibia's PBV (mm)"
      title <-  ""
    } else if(trait.temp == "masse"){
      ylab <- "Mass PBV (g)"
      title <- "(b)"
    }

    # --- plot
    plot.temp <- ggplot() +
      # HIGHLIGHT 0
      geom_hline(yintercept = 0, col = "#5C5C5C", linetype = 1, linewidth = 0.05) +
      
      
      # QUANTILE 95% OF BV FOR EACH YEAR
      # --- LOWER
      geom_smooth(
        data = ic_df,
        mapping = aes(x = as.numeric(as.character(annee_emergence)), y = ic_lwr),
        linewidth = 1,
        linetype = 3,
        se = F,
        col = "#545454"
      ) +
      
      # --- UPPER
      geom_smooth(
        data = ic_df,
        mapping = aes(x = as.numeric(as.character(annee_emergence)), y = ic_upr),
        linewidth = 1,
        linetype = 3,
        se = F,
        col = "#545454"
      ) +
      
 
      
      # TRAJECTORIES
      geom_line(
        data = mean_bv_iteration,
        mapping = aes(x = as.numeric(as.character(annee_emergence)), y = bv, group = iteration),
        linewidth = 0.05,
        alpha = 0.1,
        col = "#545454"
      ) +
      
      # POOLED
      # geom_smooth(data = mean_bv_iteration,
      #             mapping = aes(x = as.numeric(as.character(annee_emergence)), y = bv),
      #             linewidth = 1,
      #             linetype = 1,
      #             se = F,
      #             col = "#104E8B",
      #             method = "lm") +
      
      # SLOPE FOR EACH PERIOD
      geom_smooth(data = mean_bv_iteration,
                  mapping = aes(x = as.numeric(as.character(annee_emergence)), y = bv, col = period),
                  linewidth = 1.5, 
                  se = F,
                  method = "lm") +
      scale_color_manual(values = c("#cc3300", "#4e9785")) +  # c("#cc3300", "#4e9785")
      
      
      
      
      
      # # X LABELS OF DATE
      # # Set custom x-axis tick labels
      scale_x_continuous(
        breaks = seq(1997, 2023, 4),
        labels = seq(1997, 2023, 4),
        expand = c(0.05, 0.05),
        limits = c(1997, 2023)
      ) +
      #   scale_x_discrete(
      #     breaks = seq(1997, 2023, 4),
      #     labels = c("1997", "2001", "2005", "2009", "2013", "2017", "2021")
      #   ) +
      
      tidyquant::theme_tq() +
      
      # Label xlab, ylab and title
      labs(x = "Cohort",
           y = ylab,
           title = title) +
    
      theme(
        # panel.grid = element_line(color = "#EEEEE0", linetype = "dashed"),  # grid of plot
        # background color
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        plot.title = element_text(size = 18, margin = margin(b = 6)),
        # axis.line = element_line(color = "black"),
        # box around plot
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 1
        ),
        
        
        # Increase space between axis and labels
        axis.title.y = element_text(size = 18, margin = margin(0, 9, 0, 0), face = "bold"),  # space between ylab and plot
        axis.title.x = element_text(face = "bold"),
        axis.text.x = element_text(margin = margin(t = 6)),
        axis.text.y = element_text(margin = margin(r = 6)),
        
        # Text size
        text = element_text(size = 18)
      ) 

    
    assign(x = paste0("plot_", trait.temp), plot.temp)
  }
  
  gridExtra::grid.arrange(plot_masse, plot_tibia, ncol = 2)
  
  return(list(plot_masse = plot_masse, plot_tibia = plot_tibia))
}
