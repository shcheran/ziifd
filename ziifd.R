# ------------------------------------------------------------------------
# Implements the main function for running ZIIFD tool
# ------------------------------------------------------------------------
#
# Args:
#   pathfilelist1/2: path to a folder containg a list of files for each chr 
#                  for the first/second track of data
#   formulaZ1/2:   formula for modelling zero-infaleded component, track 1/2
#   formulaB1/2:   formula for modelling background component, track 1/2
#   formulaE1/2:   formula for modelling enrichment component, track 1/2
#   outputPath:    path to output folder
#   filename1/2:   prefix used to denote the output files for the firts/second
#                  track that are created by ziifd
#   threshold:     threshold of posterior probability for selecting peaks
#                  must be between 0 and 1 (default 0.95)
#   consensus:     boolean - specifies whether the consensus track output 
#                  should be created or not (default FALSE)
#   numCores:      a desirable number of cores to be used for the run 
#                  (default is 4)
#
# Returns:
#   Writes the output peaks files to the output folder
#rsem_hepg2_de_transcripts_upreg_by_snd1_pathways
# Packages required:
#   zinba, MASS, matrixStats, BSgenome, multicore, doMC
#
# ------------------------------------------------------------------------



ziifd <- function(pathfilelist1, pathfilelist2, formulaZ1 = NULL, 
                  formulaB1 = NULL, formulaE1 = NULL, formulaZ2 = NULL, 
                  formulaB2 = NULL, formulaE2 = NULL, outputPath = NULL, 
                  filename1 = NULL, filename2 = NULL, 
                  threshold = 0.95, consensus = FALSE,
                  numCores = 4) {
  
  
  # load the required libraries
  suppressMessages(library(matrixStats))
  suppressMessages(library(MASS))
  suppressMessages(library(BSgenome))
  suppressMessages(library(multicore))
  suppressMessages(library(doMC))
  suppressMessages(library(R.utils))
  
  
  # check whether all the input parameters were given
  if(!hasArg(pathfilelist1) && !file.exists(pathfilelist1)) 
    stop("Provide valid paths to input folders")
  if(!hasArg(pathfilelist2) && !file.exists(pathfilelist2)) 
    stop("Provide valid paths to input folders") 
  if(!hasArg(outputPath) && !file.exists(outputPath)) 
    stop("Provide a valid path to an output folder\n")
  
  
  # read files from the specified folders
  filelist1 <- list.files(pathfilelist1)
  filelist2 <- list.files(pathfilelist2)
  filesNumb <- length(filelist1)
  time <- Sys.time()
  cat("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n")
  cat(paste("ANALYSIS STARTED AT ", strsplit(as.character(time), " ")[[1]][2], " ",
            strsplit(as.character(time), " ")[[1]][1], "\n", sep=""))
  cat("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n")
  
  
  # access multicores packge:
  registerDoMC(cores = numCores)
  mcoptions <- list(preschedule = FALSE, set.seed = FALSE)
  cat(paste("Using ", numCores, " core(s) \n\n", sep=""))
  tmpfile <- paste0(outputPath, "/_log.file")
  write("0", file = tmpfile, append = T)
  
  
  # run simultaneous analysis for the specified files:
  # create temporary file for monitoring total numer of iterations
  ptm <- proc.time()
  pkList <- foreach(i = 1:filesNumb, .inorder = FALSE, .options.multicore = mcoptions) %dopar% 
  {
    
    # read the file names prepared for the analysis
    file1 <- file.path(pathfilelist1, filelist1[i])
    file2 <- file.path(pathfilelist2, filelist2[i])
    time <- paste0("[", strsplit(as.character(Sys.time()), " ", fixed = T)[[1]][2], "]")
    message <- paste(time, " Simultaneous analysis of tracks: \n", 
                     "      ", tail(strsplit(file1, "/")[[1]], 1), 
                     " AND\n      ", tail(strsplit(file2, "/")[[1]], 1), 
                     "\n\n", sep="")
    cat(message)
    
    
    # error handling: reading data table for track 1
    # check whether files are corrupted or not: stop execution if yes
    data1 <- tryCatch(
      read.table(file1, header = T, stringsAsFactors = F), 
      error = function(e) {
        cat(paste("File ", tail(strsplit(file1, "/")[[1]], 1), " is corrupted\n", sep=""))
        cat(paste("Run completed at ", Sys.time(), "\n", sep="")); return(0)}, 
      warning=function(w) {
        cat(paste("File ", tail(strsplit(file1, "/")[[1]], 1), " is corrupted\n",sep=""))
        cat(paste("Run completed at ", Sys.time(), "\n", sep="")); return(0)})
    
    
    # error handling: reading data table for track 2
    # check whether files are corrupted or not: stop execution if yes
    data2 <- tryCatch(
      read.table(file2, header = T, stringsAsFactors = F), 
      error = function(e) {
        cat(paste("File ", tail(strsplit(file2, "/")[[1]], 1), " is corrupted\n", sep=""))
        cat(paste("Run completed at ", Sys.time(), "\n", sep="")); return(0)}, 
      warning=function(w) {
        cat(paste("File ", tail(strsplit(file2, "/")[[1]], 1), " is corrupted\n",sep=""))
        cat(paste("Run completed at ", Sys.time(), "\n", sep="")); return(0)})
    
    if(is.null(formulaZ1) && is.null(formulaZ2) &&
       is.null(formulaB1) && is.null(formulaB2) &&
       is.null(formulaE1) && is.null(formulaE2)){
      cat("Provide the formulas for modelling ZINB components! Terminating. \n")
      return(0)
    }
    
    
    # if files are not corrupted, generate the working data frames    
    covsTrack1 <- list(formulaZ1, formulaB1, formulaE1)
    covsTrack2 <- list(formulaZ2, formulaB2, formulaE2)    
    covsTrack1 <- lapply(1:3, function(k) as.formula(covsTrack1[[k]]))
    covsTrack2 <- lapply(1:3, function(k) as.formula(covsTrack2[[k]]))
    
    
    # create model frames according to the specified formulas
    mf1 <- lapply(1:3, function(j) model.frame(formula = covsTrack1[[j]], data = data1))
    mf2 <- lapply(1:3, function(j) model.frame(formula = covsTrack2[[j]], data = data2))
    y1 <- model.extract(mf1[[1]], "response")
    y2 <- model.extract(mf2[[1]], "response")
    data1 <- cbind(data1, qVal = (y1 > quantile(y1, 0)) ^ 2)
    data2 <- cbind(data2, qVal = (y2 > quantile(y2, 0)) ^ 2)
    
    
    # get the variable names fro constructing the single model frame
    names1 <- lapply(1:3, function(j) attributes(mf1[[j]])$names)
    names2 <- lapply(1:3, function(j) attributes(mf2[[j]])$names)
    dfnames1 <- unique(unlist(names1))
    dfnames2 <- unique(unlist(names2))
    
    
    # generate data frame including column with post.prob for zero-inflated component
    df1 <- data.frame(data1[dfnames1], exp_count_zi = 1 * (y1 <= 0), stringsAsFactors = F)
    df2 <- data.frame(data2[dfnames2], exp_count_zi = 1 * (y2 <= 0), stringsAsFactors = F)
    
    
    # create covariate set for the ZI component for the first and second tracks
    terms1 <- strsplit(formulaZ1, "~")[[1]][2]
    formZ1 <- paste("exp_count_zi ~", terms1)
    terms2 <- strsplit(formulaZ2, "~")[[1]][2]
    formZ2 <- paste("exp_count_zi ~", terms2)
    
    
    # Run em-algorithm
    chrname <- data1[2,1]
    output <- em(df1, df2, formZ1, formulaB1, formulaE1, formZ2, 
                 formulaB2, formulaE2, chrname, tmpfile)
    
    
    # if algorithm didn't converge, program terminates
    if (output$convergence == 0) {
      cat("Algorithm didn't converge. Provide another set of formulas\n")
      cat(paste("Run completed at ", Sys.time(), "\n", sep=""))
      return(0)
    }
    
    # postprocessing
    mui1 <- output$glmFitBg1$fitted
    mui2 <- output$glmFitBg2$fitted
    data1$qVal[y1 - mui1 < 0] <- 0
    data2$qVal[y2 - mui2 < 0] <- 0
    
    
    # merge peaks
    peaksList <- mergepeaks(output, data1, data2, threshold, outputPath, 
                            filename1, filename2, cons=consensus, tmpfile)
    
    time <- paste0("[", strsplit(as.character(Sys.time()), " ", fixed = T)[[1]][2], "]")
    cat(paste0(time, " Run completed for chromosome ", chrname, "\n"))
    rm(output, data1, data2, mf1, mf2, y1, y2)
    return(peaksList)
  }
  
  
  # generate the output
  peaksTotal1 <- rbind()
  peaksTotal2 <- rbind()
  peaksCons <- rbind()
  for(j in 1:length(pkList)){
    print(length(pkList))
    chrListPeaks <- pkList[[j]]
    peaks1 <- chrListPeaks$pk1
    peaks2 <- chrListPeaks$pk2
    peaksTotal1 <- rbind(peaksTotal1, peaks1)
    peaksTotal2 <- rbind(peaksTotal2, peaks2)
    if(!is.null(chrListPeaks$cons)){
      peaksCons <- rbind(peaksCons, chrListPeaks$cons)
    }
  }
  
  
  # write output files
  if(!dir.exists(paste0(outputPath, "/peaks"))){
    dir.create(paste0(outputPath, "/peaks"))
  }
  fname1 <- paste("ziifd_", filename1, "_peaks.bed", sep="")
  fname2 <- paste("ziifd_", filename2, "_peaks.bed", sep="")
  write.table(peaksTotal1, file=file.path(paste0(outputPath, "/peaks/"), fname1), 
              sep="\t", row.names = F, col.names = F, quote = F)
  write.table(peaksTotal2, file=file.path(paste0(outputPath, "/peaks/"), fname2), 
              sep="\t", row.names = F, col.names = F, quote = F)
  if(consensus){
    fname <- paste("ziifd_consensus_", filename1, "_", filename2, "_peaks.bed", sep="")
    write.table(peaksCons, file=file.path(paste0(outputPath, "/peaks/"), fname), 
                sep="\t", row.names = F, col.names = F, quote = F)
  }
  tmpremove <- file.remove(tmpfile)
  
  
  # print the total time used to the screen
  tm <- proc.time() - ptm
  
  cat("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n")
  cat(paste("ANALYSIS ENDED AT ", strsplit(as.character(time), " ")[[1]][2], " ",
            strsplit(as.character(time), " ")[[1]][1], "\n", sep=""))
  if (tm[[3]] / (3600 * 24) >= 1) {  del <- 3600 * 24; time<- " days"} 
  else if (tm[[3]] / (3600) >= 1) { del <- 3600; time<- " hours"}
  else { del <- 60; time<- " minutes"}
  cat(paste("Processing time of files: ", round(tm[[3]] / del, d = 4),
            time, "\n"))
  cat("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n")
}









# ------------------------------------------------------------------------
# Executes extended EM - algorithm for two tracks data
# ------------------------------------------------------------------------
#
# Args:
#   data1/2:       data frames for the first and secod track of data containing 
#                  readcount along with quantified covariates for each window
#   formulaZ1/2:   formula for modelling zero-infaleded component, track 1/2
#   formulaB1/2:   formula for modelling background component, track 1/2
#   formulaE1/2:   formula for modelling enrichment component, track 1/2
#   chrname:       chromosome name
#
# Returns a list contating the following fields:
#   probi1/2       matrix containing three columns, each of which represents 
#                  posterior probability values for a one of model components
#   convergence:   value indicating whether the algorithm converged (1) or not (0)
#   glmFitBg1/2:   output from glm fit for background compomemt for track 1/2
#
# ------------------------------------------------------------------------


em <- function(data1 = NULL, data2 = NULL, formulaZ1 = NULL, 
               formulaB1 = NULL, formulaE1 = NULL, formulaZ2 = NULL, 
               formulaB2 = NULL, formulaE2 = NULL, chrname = NULL, 
               tmpfile = NULL){
  
  
  # extract response from the model
  y1 <- model.extract(model.frame(as.formula(formulaE1), data = data1), "response")
  y2 <- model.extract(model.frame(as.formula(formulaE2), data = data2), "response")

  
  # test multiple partitions for enrichment component and choose one that provides best fit
  # cat("    Calculating initial partitions for enrichment component for ", chrname, "\n")
  enrichProp1 <- suppressWarnings(getenrichprop(data1, formulaZ1, formulaB1, formulaE1, tmpfile))
  enrichProp2 <- suppressWarnings(getenrichprop(data2, formulaZ2, formulaB2, formulaE2, tmpfile))
  if (enrichProp1 == 0 || enrichProp2 == 0) { return(list(convergence = 0)) }
  n <- dim(data1)[1]
  cat(". ")
  

  # initialization step performed for two tracks separately
  init1 <- initialize(data1, enrichProp1, formulaZ1, formulaB1, formulaE1)
  init2 <- initialize(data2, enrichProp2, formulaZ2, formulaB2, formulaE2)
  cat(". ")
  
  
  # for console output
  globIters <- as.numeric(readLines(tmpfile))
  lastIter <- globIters[length(globIters)]
  write(as.character(lastIter + 2), file = tmpfile, append = T)
  
  
  # perform M-step for initial parameters 
  mStep1 <- list(glmFitZi = init1$glmFitZi, glmFitBg = init1$glmFitBg, glmFitE = init1$glmFitE)
  mStep2 <- list(glmFitZi = init2$glmFitZi, glmFitBg = init2$glmFitBg, glmFitE = init2$glmFitE)
  
  
  # initialize variables from initialization steps
  data1$exp_count_zi <- init1$t[,1]
  data2$exp_count_zi <- init2$t[,1]
  probi2 <- init2$t
  probi1 <- init1$t
  
  
  # Update mixture proportion values of background and enrichment components
  prior1 <- getpriors(probi1)
  prior2 <- getpriors(probi2)
  
  
  # calculate kog-likelihood
  llnew <- init1$L + init2$L    
  L <- llold <- 2 * llnew
  
  
  # Perform an iterative process until the algorithm converges
  iter <- 2
  while (abs((llold - llnew) / llold) > 1e-05) {
    
    
    # Estimation of the correlation term, i.e conditional probabiliy matrices
    conditpr <- getconditionalprob(probi1, probi2)
    
    
    # E-step
    probiTmp1 <- estep(y1, conditpr$conditProb1, mStep1$glmFitBg, mStep1$glmFitE,
                       mStep2$glmFitZi$fitted, prior2$pr1, prior2$pr2, probi2)
    probiTmp2 <- estep(y2, conditpr$conditProb2, mStep2$glmFitBg,  mStep2$glmFitE,
                       mStep1$glmFitZi$fitted, prior1$pr1, prior1$pr2, probi1) 
    probi1 <- probiTmp1
    probi2 <- probiTmp2
    data1$exp_count_zi <- probi1[,1]
    data2$exp_count_zi <- probi2[,1]
    
    
    # M-step
    mStep1 <- mstep(data1, probi1, mStep1$glmFitBg, mStep1$glmFitE, formulaZ1,formulaB1, formulaE1)
    mStep2 <- mstep(data2, probi2, mStep2$glmFitBg, mStep2$glmFitE, formulaZ2, formulaB2, formulaE2)
    
    
    # Update mixture proportion values of background and enrichment components for two tracks
    prior1 <- getpriors(probi1)
    prior2 <- getpriors(probi2)
    if (prior1$pr1 <= .5 || prior2$pr1 <= .5) { 
      cat("The estimated proportion of enrichment component has exceeded 0.5.\n")
      return(list(convergence = 0)) 
    }
    
    
    # Compute complete log-likelihood of a model
    llold <- L[max(1, iter - 10)]
    
    loglik1 <- loglikelihood(y1, mStep1$glmFitBg$fitted, mStep1$glmFitE$fitted,
                             mStep1$glmFitZi$fitted, prior1$pr1, prior1$pr2, 
                             mStep1$glmFitBg$theta, mStep1$glmFitE$theta)
    
    loglik2 <- loglikelihood(y2, mStep2$glmFitBg$fitted, mStep2$glmFitE$fitted,
                             mStep2$glmFitZi$fitted, prior2$pr1, prior2$pr2, 
                             mStep2$glmFitBg$theta, mStep2$glmFitE$theta)
    
    llnew <-  loglik1 + loglik2
    if (iter > 300) { 
      return(list(convergence = 0)) 
    }
    
    
    # for console printing in multiple processors mode
    globIters <- as.numeric(readLines(tmpfile))
    lastIter <- globIters[length(globIters)]
    
    if(lastIter %% 50 == 0){
      lastIter <- 0
      cat('\n')
    }
    
    write(as.character(lastIter + 1), file = tmpfile, append = T)
    
    L <- c(L, llnew)
    iter <- iter + 1
    cat(". ")

  }
  cat("\n")
  convergence <- 1
  
  
  # return the list
  list(probi1 = probi1, probi2 = probi2, convergence = convergence,
       glmFitBg1 = mStep1$glmFitBg, glmFitBg2 = mStep2$glmFitBg)
}









# ------------------------------------------------------------------------
# Gets initial partition of inrichment component
# ------------------------------------------------------------------------
#
# Args:
#   df:       Data frame contating the data including read couns and 
#             covariates quantified for each window
#   formulaZ: Formula for modelling zero-inflated component
#   formulaB: Formula for modelling background component
#   formulaE: Formula for modelling enrichment component
#
# Returns:
#   Proportion of enrichment windowns in the analysed track of data
#
# ------------------------------------------------------------------------

getenrichprop <- function(df, formulaZ, formulaB, formulaE, tmpfile=NULL){
  
  
  # iitialize variables
  maxit <- 26
  probs <- data.frame(val = c(0.15, 0.11275, 0.0755, 0.03825, 0.001),
                      flag = rep(1,5), stringsAsFactors = F)
  
  formulaZ <- as.formula(formulaZ)
  formulaB <- as.formula(formulaB)
  formulaE <- as.formula(formulaE)
  
  
  # create data frames
  mfE <- model.frame(formulaE, data = df)
  mfB <- model.frame(formulaB, data = df)
  y <- as.numeric(model.extract(mfE, "response")) 
  n <- length(y)
  
  
  # test multiple proportions of enrichment component
  loglike <- c()
  for (i in 1:length(probs[,1])) {
    
    
    # initialize variables
    t <- matrix(0, n, 3)

    
    # renew the variables for the zero-inflated component for checking each partitions
    df$exp_count_zi <- 1 * (y <= 0)
    mfZ <- model.frame(formulaZ, data = df)
    
    
    # proportion of enriched windows
    pr2 <- probs[i,]$val
    
    
    # proportion of background windows
    pr1 <- 1 - pr2
    n1  <- round(length(y) * (1 - pr2))
    priorWeight <- rep(1 - 10^-10, length(y))
    sortedInd <- sort(y, index.return = T)$ix
    priorWeight[sortedInd[1:n1]] <- 10^-10    
    
    
    # model zero-inflated component using logistic regression
    mustartzi <- (rep(1,n) * (y <=0 ) + 0.5) / (rep(1,n) + 1)
    
    glmFitZi <- suppressWarnings(
      glm(as.formula(formulaZ), control = glm.control(maxit = maxit), 
          family = "binomial", mustart = mustartzi, data = mfZ)
      )
    pr0 <- glmFitZi$fitted 
    
    
    # model background with poisson regression to get starting values   
    glmFitBg <- glm(as.formula(formulaB), weights = 1 - priorWeight, 
                     family = poisson, control = glm.control(maxit = maxit),
                     mustart = y+(y == 0) / 6, data = mfB)
  
    
    # model enrichment with poisson regression to get starting values     
    glmFitE <- glm(as.formula(formulaE), weights = priorWeight, family = poisson, 
              control = glm.control(maxit = maxit), 
              mustart = y + (y == 0) / 6, data = mfE)  
    
    
    # calculate log-likelihood of a model
    ll <- loglikelihood(y, glmFitBg$fitted, glmFitE$fitted, pr0, pr1, pr2, 1, 1)
    
    
    # e-step: calculate posterior probabilities
    f1 <- dnbinom(y, size = 1, mu = glmFitBg$fitted)
    f2 <- dnbinom(y, size = 1, mu = glmFitE$fitted)
    
    t[, 1] <- pr0 / (pr0 + (1 - pr0) * (pr1 * f1 + pr2 * f2))
    t[, 2] <- ((1 - pr0) * pr1 * f1) / (pr0*(y <= 0) + (1 - pr0) * (pr1 * f1 + pr2 * f2))
    t[, 3] <- ((1 - pr0) * pr2 * f2) / (pr0*(y <= 0) + (1 - pr0) * (pr1 * f1 + pr2 * f2))
    t[y > 0, 1] <- 0
    
    
    # chech NAs in calculated resulted values
    NAs<-which(t[, 2] == 'NaN'| t[, 3] == 'NaN')      
    if(length(NAs > 0)){
      t[, 2][NAs] <- 0
      t[, 3][NAs] <- 1
    } 
    
    
    # update mixture proportion values
    pr1 <- sum(t[, 2]) / sum(t[, 2] + t[, 3])
    pr2 <- sum(t[, 3]) / sum(t[, 2] + t[, 3])
    if (pr1 <= .5) {  probs[i,]$flag <- 0 }
    
    
    # update values for zero-inflated component
    df$exp_count_zi <- t[, 1]
    mfZ <- model.frame(formulaZ, data = df)
    mustartzi <- (rep(1,n) * t[, 1] + 0.5) / (rep(1,n) + 1)
    
    
    # temporary variables
    glmFitBg_ <- glmFitBg 
    glmFitE_ <- glmFitE 
    
    
    # glm fit for zero-inflated component
    glmFitZi <- suppressWarnings(
      glm(formulaZ, family = "binomial", control = glm.control(maxit = maxit), 
          mustart = mustartzi, data=mfZ)
      )
    pr0 <- glmFitZi$fitted
    
    
    # glm fit for background component
    glmFitBg <- tryCatch(glm.nb(formulaB, weights = t[, 2], 
                          control = glm.control(maxit = maxit),
                          init.theta = 1, mustart = glmFitBg$fitted, data = mfB),
                   error=function(e) e, warning=function(w) w)
    
    
    # glm fit for enrichment component
    glmFitE <- tryCatch(glm.nb(formulaE, weights = t[, 3], 
                          control = glm.control(maxit = maxit), 
                          init.theta = 1, mustart = glmFitE$fitted, data = mfE),
                  error=function(e) e, warning=function(w) w)
    
    
    # possible warning messages
    w1 <- "step size truncated due to divergence"
    w2 <- "NA/NaN/Inf in 'x'"
    w3 <- "glm.fit: algorithm did not converge"
    w4 <- "alternation limit reached"
    
    
    # check the warnings, errors for the background component
    if(is(glmFitBg,"warning") || is(glmFitBg,"error")) { 
      
      cat("\nBackground component warning occured for proportion", probs[i, 1],
          "Fit simplified negbin regression.\n")
      
      # for console printing in multiple processors mode
      write('0', file = tmpfile, append = T)
      

    if((glmFitBg$message == w1 || glmFitBg$message == w2) || glmFitBg$message == w3 ){
      
      glmFitBg <- glm.nb(formulaB, control = glm.control(maxit = maxit), 
                         init.theta = 1, mustart = glmFitBg_$fitted, data = mfB)
      
      } else if(glmFitBg$message == w4) { 
        
        # do nothing if algorithm did not converge
        glmFitBg <- glm.nb(formulaB, weights = t[, 2], init.theta = 1, 
                           control = glm.control(maxit = maxit), 
                           mustart = glmFitBg_$fitted, data = mfB)
      } else {
        
        # or "iteration limit reached"
        glmFitBg <- glm.nb(formulaB, weights = t[, 2], init.theta = 1, 
                           mustart = glmFitBg_$fitted, data = mfB)
      }
    }
    
    
    # check the warnings, errors for the enrichment component
    if(is(glmFitE,"warning") || is(glmFitE,"error")) { 
      
      cat("\nEnrichment component warning occured for proportion", probs[i, 1], 
          "Fit simplified negbin regression.\n")
      
      # for console printing in multiple processors mode
      write('0', file = tmpfile, append = T)
      
      if ((glmFitE$message == w1 || glmFitE$message == w2) || glmFitE$message == w3) {
        
        glmFitE <- glm.nb(formulaE, control = glm.control(maxit = maxit), 
                    init.theta = 1, mustart = glmFitE_$fitted, data = mfE)
        
      } else if  (glmFitE$message == w4){ 

        # do nothing if algorithm did not converge
        glmFitE <- glm.nb(formulaE, weights = t[, 3], init.theta = 1, 
                          control = glm.control(maxit = maxit), 
                          mustart = glmFitE_$fitted, data = mfE)
      } else {
        # or "iteration limit reached"
        glmFitE <- glm.nb(formulaE, weights = t[, 3], init.theta = 1, 
                    mustart = glmFitE_$fitted, data = mfE)
      }
    }
    
    
    # Compute complete log-likelihood of a model
    loglike[i] <- loglikelihood(y, glmFitBg$fitted, glmFitE$fitted, pr0, pr1, 
                                pr2, glmFitBg$theta, glmFitE$theta)
  }
  
  loglikeSumm <- data.frame(loglik = loglike, val = probs$val, flag = probs$flag, 
                            stringsAsFactors = F)
  loglikeSumm <- loglikeSumm[order(loglikeSumm$loglik, decreasing = T), ] 
  probsSorted <- loglikeSumm[,2:3]
  
  
  # format the output
  if (dim(probsSorted[which(probsSorted$flag == 1),])[1] != 0) {
    prop <- probsSorted[which(probsSorted$flag == 1),][1,]$val
  } else { 
    prop <- 0
    }
  
  
  # return
  prop
}









# ------------------------------------------------------------------------
# Gets initialization parameters for running EM - algorithm 
# ------------------------------------------------------------------------
#
# Args:
#   df:       data frame contating the data including read counts and 
#             covariates quantified for each window
#   formulaZ: formula for modelling zero-inflated component
#   formulaB: formula for modelling background component
#   formulaE: formula for modelling enrichment component
#
# Returns a list of arguments including:
#   t:        posterior probability matrix 
#   glmFitZi: glm fit for zero- inflated component
#   glmFitBg: glm fit for background component
#   glmFitE:  glm fit for enrichmnet component
#   L:        estimated model likelihood
#
# ------------------------------------------------------------------------


initialize <- function(df, prop, formulaZ, formulaB, formulaE){
  
  
  maxit <- 26
  
  formulaZ <- as.formula(formulaZ)
  formulaB <- as.formula(formulaB)
  formulaE <- as.formula(formulaE)
  
  mfE <- model.frame(formulaE, data = df)
  mfB <- model.frame(formulaB, data = df)
  
  y <- model.extract(mfE, "response")
  y0 <- (y == 0)
  n <- length(y)
  t <- matrix(0, n, 3)
  
  
  # zero-inflated component
  df$exp_count_zi <- 1 * (y <= 0)
  mfZ <- model.frame(formulaZ, data = df)
  
  
  # initialize mixture proportion variables according to best partition
  pr2 <- prop
  pr1 <- 1 - pr2
  n1  <- round(length(y) * (1 - pr2))
  priorWeight <- rep(1 - 10 ^ -10, length(y))
  sorted.ind <- sort(y, index.return = T)$ix
  priorWeight[sorted.ind[1:n1]] <- 10 ^ -10
  
  
  # glm fit for zero-inflated component
  mustartzi <- (rep(1,n) * (y <=0 ) + 0.5) / (rep(1,n) + 1)
  glmFitZi <- suppressWarnings(
    glm(formulaZ, family = "binomial", control = glm.control(maxit = maxit), 
        mustart = mustartzi, data=mfZ))
  pr0 <- glmFitZi$fitted
  
  
  # glm fit for background component
  glmFitBg <- glm(formulaB, weights = 1-priorWeight, family = poisson,
                  control = glm.control(maxit = maxit), mustart = y + y0/6, 
                  data = mfB)
  
  
  # glm fit for enrichment component
  glmFitE <- glm(formulaE, weights = priorWeight, family = poisson, 
                 control = glm.control(maxit = maxit),  mustart = y + y0/6, 
                 data = mfE)
  
  
  # E-step: calculate posterior probabilities
  f1 <- dnbinom(y, size = 1, mu = glmFitBg$fitted)
  f2 <- dnbinom(y, size = 1, mu = glmFitE$fitted)
  
  
  # update posterior probabilities for each component
  t[, 1] <- pr0/(pr0 + (1 - pr0)*(pr1*f1 + pr2*f2))
  t[, 2] <- ((1 - pr0)*pr1*f1)/(pr0*(y <= 0) + (1 - pr0)*(pr1*f1 + pr2*f2))
  t[, 3] <- ((1 - pr0)*pr2*f2)/(pr0*(y <= 0) + (1 - pr0)*(pr1*f1 + pr2*f2))
  t[y > 0, 1] <- 0
  
  
  # handle NA values
  NAs<-which(t[, 2]=='NaN'| t[, 3]=='NaN')      
  if(length(NAs > 0)){
    t[, 2][NAs] <- 0
    t[, 3][NAs] <- 1
  }
  
  
  # update the variable
  df$exp_count_zi <- t[, 1]
  mfZ <- model.frame(formulaZ, data=df)
  pr1 <- sum(t[, 2]) / sum(t[, 2] + t[, 3])
  pr2 <- sum(t[, 3]) / sum(t[, 2] + t[, 3])
  
  
  # M-step: weighted glm fit using posterior probability values as weights
  # zero-inflated component
  mustartzi <- (rep(1,n) * t[,1] + 0.5) / (rep(1,n) + 1)
  glmFitZi <- suppressWarnings(
    glm(formulaZ, family = "binomial", control = glm.control(maxit = maxit), 
        mustart = mustartzi, data = mfZ))
  pr0 <- glmFitZi$fitted
  glmFitZi <- list(fitted = glmFitZi$fitted, coefficients = glmFitZi$coefficients)
  
  
  # background component
  glmFitBg <- suppressWarnings(
    glm.nb(formulaB, weights = t[, 2], control = glm.control(maxit = maxit), 
           init.theta = 1, mustart = glmFitBg$fitted, data = mfB))
  glmFitBg <- list(theta = glmFitBg$theta, fitted = glmFitBg$fitted, 
                   coefficients = glmFitBg$coefficients) 
  
  
  # enrichment component
  glmFitE <- suppressWarnings(
    glm.nb(formulaE, weights = t[, 3], control = glm.control(maxit = maxit), 
           init.theta = 1, mustart = glmFitE$fitted, data = mfE))
  glmFitE <- list(theta = glmFitE$theta, fitted = glmFitE$fitted, 
                  coefficients = glmFitE$coefficients)
  
  
  # Compute complete log-likelihood of a model
  L <- loglikelihood(y, glmFitBg$fitted, glmFitE$fitted, pr0, pr1, pr2, 
                     glmFitBg$theta, glmFitE$theta)
  
  
  # return
  list(t = t, glmFitZi = glmFitZi, glmFitBg = glmFitBg, glmFitE = glmFitE, L = L)
}









# ------------------------------------------------------------------------
# Calculate matrices of conditional probabilities
# ------------------------------------------------------------------------
#
# Args:
#   probi1/2:    matrices of posterior probabilities containing three columns
#                each corresponding to one of the three model components: 
#                zero-inflated, background and enrichment
#
# Returns a list containing:
#   conditProb1/2: matrix of conditional probabilities for track 1/2
#   jointProb:     joined probability matrix
#
# ------------------------------------------------------------------------


getconditionalprob <- function(probi1, probi2){
  
  # initialize variables
  M1 <- matrix(0, 3, 3)
  M2 <- matrix(0, 3, 3)
  
  
  # calculate joint probability matrix
  jointProb<-(t(probi1) %*% probi2) / length(probi1[, 1])
  
  
  # calculate conditional probability matrix for the first track
  M1[,1] <- jointProb[, 1] / sum(jointProb[, 1])
  M1[,2] <- jointProb[, 2] / sum(jointProb[, 2])
  M1[,3] <- jointProb[, 3] / sum(jointProb[, 3])
  M1[M1 == 'NaN'] <- 0
  
  
  # calculate conditional probability matrix for the second track
  M2[1,] <- jointProb[1, ] / sum(jointProb[1, ]) 
  M2[2,] <- jointProb[2, ] / sum(jointProb[2, ]) 
  M2[3,] <- jointProb[3, ] / sum(jointProb[3, ])
  M2[M2 == 'NaN'] <- 0; M2 <- t(M2) 
  
  
  # return list
  list(conditProb1 = M1, conditProb2 = M2, jointProb = jointProb)
}









# ------------------------------------------------------------------------
# Executes E-step of the extended EM-algorithm
# ------------------------------------------------------------------------
#
# Args:
#   y: 
#   conditpr: 
#   glmFitBg: 
#   glmFitE: 
#   priorGiven0: 
#   priorGiven1: 
#   priorGiven2: 
#   probiMatrGiven: 
#
# Returns:
#   A matrix of posterior probability vectors calculated based on the 
#   parameters and posterior probability distributions of the other 
#   ('given') track
#
# ------------------------------------------------------------------------


estep <- function(y, conditpr, glmFitBg, glmFitE, priorGiven0, priorGiven1, 
                  priorGiven2, probiMatrGiven){
  
  n <- length(y)
  tmpMatr <- list(matrix(0, n, 3), matrix(0, n, 3), matrix(0, n, 3))
  
  
  # calculate the likelihood of the negative binomial functions
  f0 <- 1 * (y <= 0)
  f1 <- dnbinom(y, size = glmFitBg$theta, mu = glmFitBg$fitted)
  f2 <- dnbinom(y, size = glmFitE$theta, mu = glmFitE$fitted)
  
  
  # calculate the denominator for expectation values
  T1 <- conditpr[1,1] * f0 * priorGiven0 + 
        conditpr[2,1] * f1 * priorGiven0 + 
        conditpr[3,1] * f2 * priorGiven0
  
  T2 <- conditpr[1,2] * f0 * (1 - priorGiven0) * priorGiven1 + 
        conditpr[2,2] * f1 * (1 - priorGiven0) * priorGiven1 + 
        conditpr[3,2] * f2 * (1 - priorGiven0) * priorGiven1
  
  T3 <- conditpr[1,3] * f0 * (1 - priorGiven0) * priorGiven2 + 
        conditpr[2,3] * f1 * (1 - priorGiven0) * priorGiven2 + 
        conditpr[3,3] * f2 * (1 - priorGiven0) * priorGiven2
  TT <- cbind(T1, T2, T3)
  
  
  # probabilities
  pr <- cbind(priorGiven0, (1 - priorGiven0) * priorGiven1, 
              (1 - priorGiven0) * priorGiven2)
  for (i in 1:3){
    tmpMatr[[1]][,i] <- (conditpr[1,i] * pr[,i] * f0) / TT[,i]
    tmpMatr[[2]][,i] <- (conditpr[2,i] * pr[,i] * f1) / TT[,i]
    tmpMatr[[3]][,i] <- (conditpr[3,i] * pr[,i] * f2) / TT[,i]
  }
  

  # create a new matrix with estimated posterior probability values
  probiMatrNew <- matrix(0, n, 3)
  for (i in 1:3){
    probiMatrNew[, i] <- tmpMatr[[i]][,1] * probiMatrGiven[,1] + 
                         tmpMatr[[i]][,2] * probiMatrGiven[,2] + 
                         tmpMatr[[i]][,3] * probiMatrGiven[,3]
  }
  
  
  # handle NA data to avoid errors
  NAs <- which(probiMatrNew[, 2] == 'NaN'| probiMatrNew[, 3] == 'NaN')
  if(length(NAs > 0)){
    probiMatrNew[, 2][NAs] <- 10 ^ -10
    probiMatrNew[, 3][NAs] <- 1
  }
  probiMatrNew[, 2][probiMatrNew[, 2] < 10 ^ -10] <- 10 ^ -10
  probiMatrNew[, 3][probiMatrNew[, 3] < 10 ^ -10] <- 10 ^ -10
  probiMatrNew[probiMatrNew[, 1] == 'NaN', 1] <- 0
  

  # return value
  probiMatrNew
}









# ------------------------------------------------------------------------
# Executes M-step of the extended EM-algorithm
# ------------------------------------------------------------------------
#
# Args:
#   data: 
#   t: 
#   glmFitBg: 
#   glmFitE: 
#   formulaZ: 
#   formulaB: 
#   formulaE: 
#
# Returns:
#   A matrix of posterior probability vectors calculated based on the 
#   parameters and posterior probability distributions of the other 
#   ('given') track
#
# ------------------------------------------------------------------------


mstep <- function(data, t, glmFitBg, glmFitE, formulaZ, formulaB, formulaE){
  
  maxit <- 26
  formulaZ <- as.formula(formulaZ)
  formulaB <- as.formula(formulaB)
  formulaE <- as.formula(formulaE)
  
  mfE <- model.frame(formulaE, data = data)
  mfB <- model.frame(formulaB, data = data)
  mfZ <- model.frame(formulaZ, data = data)
  
  y <- model.extract(mfE, "response")
  n <- length(y)

  
  # glm fit for zero-infl component using weights estimated on previous iteration
  mustartzi <- as.double((rep(1,n)*t[,1] + 0.5)/(rep(1,n) + 1))
  glmFitZi <- suppressWarnings(
    glm(formulaZ, family = "binomial", control = glm.control(maxit = maxit), 
        mustart = mustartzi, data=mfZ))
  
  
  # glm fit for background component using weights estimated on previous iteration
  glmFitBg <- suppressWarnings(
    glm.nb(formulaB, weights = t[, 2], control = glm.control(maxit = maxit), 
           init.theta = glmFitBg$theta, mustart = glmFitBg$fitted, data=mfB))
  
  
  # glm fit for enrich component using weights estimated on previous iteration
  glmFitE <- suppressWarnings(
    glm.nb(formulaE, weights = t[, 3], control = glm.control(maxit = maxit), 
           init.theta = glmFitE$theta, mustart = glmFitE$fitted, data=mfE))
  
  
  # update the objects
  glmFitZi <- list(fitted = glmFitZi$fitted, 
                   coefficients = glmFitZi$coefficients)
  
  glmFitBg <- list(theta = glmFitBg$theta, 
                   fitted = glmFitBg$fitted, 
                   coefficients = glmFitBg$coefficients) 
  
  glmFitE <- list(theta = glmFitE$theta, 
                  fitted = glmFitE$fitted, 
                  coefficients = glmFitE$coefficients)
  
  
  # return list
  return(list(glmFitZi = glmFitZi, glmFitBg = glmFitBg, glmFitE = glmFitE))
}









# ------------------------------------------------------------------------
# Calculates prior probabilities
# ------------------------------------------------------------------------
#
# Args:
#   t: a matrix of posterior probabilities contatining three columns each 
#      corresponding to one of three model components, i.e. zero-inflated,
#      background and enrichment components
#
# Returns:
#   A list of prior probabilies
#
# ------------------------------------------------------------------------


getpriors <- function (t){
  pr1 <- sum(t[,2]) / sum(t[,2] + t[,3])
  pr2 <- sum(t[,3]) / sum(t[,2] + t[,3])
  list(pr1 = pr1, pr2 = pr2)
}









# ------------------------------------------------------------------------
# Calculates log-likelihood of a model
# ------------------------------------------------------------------------
#
# Args:
#   y:      vector contating read counts quantified across genomic windows
#   mu1:    fitted mean parameter for background component
#   mu2:    fitted mean parameter for enrichmemt component
#   pr0:    vector of prior probabilities for zero-inflated component
#   pr1:    prior probability value for background component
#   pr2:    prior probability value for enrichmemt component
#   theta1: overexpression parameter for background component 
#   theta2: overexpression parameter for background component 
#
# Returns:
#   Log likelihood value calculated for the given parameters
#
# ------------------------------------------------------------------------


loglikelihood<-function(y, mu1, mu2, pr0, pr1, pr2, theta1, theta2){
  y0 <- (y <= 0)
  y1 <- (y > 0)
  loglik0 <- log(pr0 + (1 - pr0) * pr1 * dnbinom(0, size = theta1, mu = mu1) + 
                   (1 - pr0) * pr2 * dnbinom(0, size = theta2, mu = mu2))
  loglik1 <- log((1 - pr0) * pr1 * dnbinom(y, size = theta1, mu = mu1) + 
                   (1 - pr0) * pr2 * dnbinom(y, size = theta2, mu = mu2))
  if(length(which(loglik1 == -Inf) > 0)){
    loglik1[which(loglik1 == -Inf)] <- 
      apply(cbind(log((1 - pr0[which(loglik1 == -Inf)]) * pr1) + 
                    dnbinom(y[which(loglik1==-Inf)], size = theta1,
                            mu = mu1[which(loglik1==-Inf)], log=TRUE),
                  log((1-pr0[which(loglik1 == -Inf)]) * pr2) + 
                    dnbinom(y[which(loglik1 == -Inf)], size = theta2, 
                            mu = mu2[which(loglik1==-Inf)], log=TRUE)), 1, 
            logSumExp)
  }
  loglike <- sum(loglik0 * y0 + loglik1 * y1)
  loglike
}









# ------------------------------------------------------------------------
# Calculates FDR and creates an output of the results
# ------------------------------------------------------------------------
#
# Args:
#   output:       vector contating posterior membership vectors obtained 
#                 from EM algorithm 
#   data1/2:      data frame contating all data for track 1/2
#   threshold:    threshold of posterior probability for selecting peaks
#   outputPath:   path to output folder
#   filename1/2:  output file names 
#   cons:         boolean - specifies whether the consensus track output 
#                 should be created or not, default - FALSE
#
# Returns:
#   Writes the output peaks files to the output folder
#
# ------------------------------------------------------------------------


mergepeaks <- function (output, data1, data2, threshold, outputPath, 
                        filename1, filename2, cons = FALSE, tmpfile){
  
  
  # load the required libraries
  suppressMessages(library(matrixStats))
  suppressMessages(library(MASS))
  suppressMessages(library(BSgenome))
  suppressMessages(library(multicore))
  suppressMessages(library(doMC))
  suppressMessages(library(R.utils))
  

  # create a consensus track if specified
  if (cons == T){
    
    averProb <- .5 * (output$probi1[, 3] + output$probi2[, 3])
    pp <- 1 - averProb
    fdr <- rep(0, length(pp))
    fdr[order(pp)] <- cumsum(pp[order(pp)]) / (1:length(pp))
    fdr[fdr <= 0] <- 0
    data <- data.frame(data1[, 1:3], qVal1 = data1$qVal, 
                  qVal2 = data2$qVal, probi = averProb, fdr = fdr,
                  stringsAsFactors = F)
    
    # consensus track
    consensus <- data[data$probi >= threshold, ]
    consensus <- consensus[which(consensus$qVal1 > 0 | consensus$qVal2 > 0), ]
    cons <- GRanges(consensus[, 1], IRanges(consensus[, 2], consensus[, 3]), 
                    data = consensus[,4:ncol(consensus)])
    cons <- sort(cons)
    consUnion <- union(cons, cons)
    chrname <- as.character(unique(seqnames(cons)))
    winlen <- sort(unique(width(cons)), decreasing = F)[1]
    cat("Creating consensus track for chromosome", chrname, "\n")
    
    
    # create peaks as the output: first track
    d <- as.data.frame(distanceToNearest(cons, consUnion, select = "all"))
    tmp <- cons[d$queryHits]
    meanProbi <- aggregate(tmp$data.probi, list(d$subjectHits), mean)
    minFdr <- aggregate(tmp$data.fdr, list(d$subjectHits), min)
    peaksCons <- data.frame(as.data.frame(consUnion)[, 1:3], prob = meanProbi$x, fdr = minFdr$x)
    
    
    # for console printing in multiple processors mode
    write('0', file = tmpfile, append = T)
    
    
    # write peaks of the consensus track into a separate file
    # write output files
    if(!dir.exists(paste0(outputPath, "/wins"))){
      dir.create(paste0(outputPath, "/wins"))
    }
    # write output files
    if(!dir.exists(paste0(outputPath, "/peaks"))){
      dir.create(paste0(outputPath, "/peaks"))
    }
    write.table(data, sep="\t", row.names = F, col.names = T, quote = F, 
                file=file.path(paste0(outputPath, "/wins/"), paste0("ziifd_consensus_", chrname, "_wins.txt")))
  } else {
    peaksCons = NULL
  }
  
  
  # create new files
  data1 <- data.frame(data1, probi = output$probi1[,3], stringsAsFactors = F)
  data2 <- data.frame(data2, probi = output$probi2[,3], stringsAsFactors = F)
  
  
  # calculating FDR: track 1
  probi1 <- output$probi1[, 3]
  p1 <- 1 - probi1
  fdr1 <- rep(0, length(p1))
  fdr1[order(p1)] <- cumsum(p1[order(p1)]) / (1:length(p1))
  fdr1[fdr1 <= 0] <- 0
  data1 <- data.frame(data1, fdr = fdr1, stringsAsFactors = F)
  
  
  # calculate FDR: track 2
  probi2 <- output$probi2[, 3]
  p2 <- 1 - probi2
  fdr2 <- rep(0,length(p2))
  fdr2[order(p2)] <- cumsum(p2[order(p2)]) / (1:length(p2))
  fdr2[fdr2 <= 0] <- 0
  data2 <- data.frame(data2, fdr = fdr2, stringsAsFactors = F)
  
  
  # write output files
  chrname <- unique(data1[,1])
  time <- paste0("[", strsplit(as.character(Sys.time()), " ", fixed = T)[[1]][2], "]")
  cat(time, "Creating output files for chromosome", chrname, "\n")
  fname1 <- paste("ziifd_", filename1, "_", chrname, "_wins.txt", sep="")
  fname2 <- paste("ziifd_", filename2, "_", chrname, "_wins.txt", sep="")
  # write output files
  if(!dir.exists(paste0(outputPath, "/wins"))){
    dir.create(paste0(outputPath, "/wins"))
  }
  # write output files
  if(!dir.exists(paste0(outputPath, "/peaks"))){
    dir.create(paste0(outputPath, "/peaks"))
  }
  write.table(data1, file = file.path(paste0(outputPath, "/wins/"), fname1), sep = "\t", 
              row.names = F, col.names = T, quote = F)
  write.table(data2, file = file.path(paste0(outputPath, "/wins/"), fname2), sep = "\t", 
              row.names = F, col.names = T, quote = F)
  
  # for console printing in multiple processors mode
  write('0', file = tmpfile, append = T)
  
  # Leave enriched windows with posterior probability higher than a specified threshold
  data1 <- data1[which(data1$probi >= threshold), ]
  data2 <- data2[which(data2$probi >= threshold), ]
  
  
  # filter out all zero count windows
  data1 <- data1[which(data1$qVal > 0), ]
  data2 <- data2[which(data2$qVal > 0), ]
  
  # save the results of run and merge peaks
  time <- paste0("[", strsplit(as.character(Sys.time()), " ", fixed = T)[[1]][2], "]")
  cat(paste0(time, " Merging peaks . . . "))
  pk1 <- GRanges(data1[,1], IRanges(data1[,2], data1[,3]), data = data1[,4:ncol(data1)])
  pk2 <- GRanges(data2[,1], IRanges(data2[,2], data2[,3]), data = data2[,4:ncol(data2)])
  pk1 <- sort(pk1)
  pk2 <- sort(pk2)
  pkUnion1 <- union(pk1, pk1)
  pkUnion2 <- union(pk2, pk2)
  
  
  # create peaks as the output: first track
  d1 <- as.data.frame(distanceToNearest(pk1, pkUnion1, select = "all"))
  tmp1 <- pk1[d1$queryHits]
  meanProbi1 <- aggregate(tmp1$data.probi, list(d1$subjectHits), mean)
  minFdr1 <- aggregate(tmp1$data.fdr, list(d1$subjectHits), min)
  peaksFinal1 <- data.frame(as.data.frame(pkUnion1)[, 1:3],
                            prob = meanProbi1$x, 
                            fdr = minFdr1$x)
  
  # create peaks as the output: second track
  d2 <- as.data.frame(distanceToNearest(pk2, pkUnion2, select = "all"))
  tmp2 <- pk2[d2$queryHits]
  meanProbi2 <- aggregate(tmp2$data.probi, list(d2$subjectHits), mean)
  minFdr2 <- aggregate(tmp2$data.fdr, list(d2$subjectHits), min)
  peaksFinal2 <- data.frame(as.data.frame(pkUnion2)[, 1:3], 
                            prob = meanProbi2$x, 
                            fdr = minFdr2$x)
  
  
  # for console printing in multiple processors mode
  write('0', file = tmpfile, append = T)
  cat("Done!\n")
  
  return(list(pk1 = peaksFinal1, pk2 = peaksFinal2, cons = peaksCons))
}









# ------------------------------------------------------------------------
# Fetching quantified window data using ZINBA 
# buildwindowdata and generateAlignability functions
# ------------------------------------------------------------------------
#
# Args:
#   pathToZinbaLibLoc:  path to zinba library. Default NULL.
#   seq:                path to aligned sample reads in BED file format 
#   inputseq:           path to aligned input reads in BED file format. If 
#                       left blank, then defaults to 'none'. 
#   generateAlign:      boolean value specifying whether the alignability
#                       files should be generated or not
#   inputDirMap:        input directory containing files with mappability 
#                       (if generateAlign = TRUE)
#   outdirAlign:        output directory for writing created output aligna-
#                       bility files (if generateAlign = TRUE)
#   inputDirAlign:      input directory with ready alignability files (if 
#                       generateAlign = FALSE)
#   outdirWins:         output directory for saving window data
#   numProc:            number of processing units. Defaults is 4.
#
# Args for ZINBA functions:
#   twobitfile:         path to build of the genome the reads were mapped to, 
#                       in .2bit format
#   extension:          average length of fragments in fragment library used
#   athresh:        	  Uniqueness threshold, number of occurrences of a given
#                       k-mer imposed during alignment (1 is absolute uniqueness). 
#   winSize:            Window size for build window data, default=250bp
#   offset:             Bp offset, default is no offsets (offset=0). 
#                       If one's winSize is 500 and would like 4 equally spaced 
#                       offsets, then specifying offset=125 would achieve this
#   cnvWinSize:         Size of windows used to calculate CNV activity in 
#                       sample, default size is 100000
#   cnvOffset:          Offset for CNV windows, typically 2500bp is suffcient
#                       although default is no offsets
#
# Returns:
#   Function uses the default setting of zinba parameters for generating 
#   window data.
#
# ------------------------------------------------------------------------


buildwins <- function(pathToZinbaLibLoc = NULL, seq = NULL, inputseq = NULL,
                      generateAlign = NULL, inputDirMap = NULL, outdirAlign = NULL, 
                      inputDirAlign = NULL,  outdirWins =  NULL, twobitfile = NULL,
                      extension = NULL, athresh = 1, winSize = 250, offset = 0,
                      cnvWinSize = 100000, cnvOffset = 2500, numProc = 4){
  
  
  # load the required libraries
  suppressMessages(library(matrixStats))
  suppressMessages(library(MASS))
  suppressMessages(library(BSgenome))
  suppressMessages(library(multicore))
  suppressMessages(library(doMC))
  suppressMessages(library(R.utils))
  
  
  # activate zinba library
  library(zinba, lib.loc = pathToZinbaLibLoc) 
  
  
  # check which files are provided
  if(is.null(inputDirMap) && is.null(inputDirAlign)){
    cat("Either path to raw mappability OR alignability files 
        must be provided! Neither is provided, terminating.\n")
    return(0)
  }
  
  
  # check which files are provided
  if (is.null(inputDirAlign) && !generateAlign) {
    cat("Alignability folder doesn't exist althogh specified 
        generateAlign = FALSE! Terminating.\n")
    return(0)
  }
  
  
  # check input file
  if (is.null(inputseq)) {
    inputseq <- "none"
  }
  
  
  # create new folder for alignability
  outdirWins <-  paste0(outdirWins, "_files")
  if(dir.exists(outdirWins)){
    cat('Overwriting existing directory at', outdirWins, "\n")
    unlink(outdirWins, recursive = TRUE)
  }
  dir.create(outdirWins, recursive = T)
  
  
  # check availability of the twoBitFile
  if(is.null(twobitfile)){
    cat("Two bit file is not provided! Terminating.\n")
    return(0)
  }
  
  
  # generate alignability if specified with ZINBA function
  if(generateAlign){
    generateAlignability(
      mapdir = paste0(inputDirMap, "/"),
      outdir = outdirAlign,
      athresh = atresh,
      extension = extension,
      twoBitFile = twobitfile
    )
    
    inputDirAlign <- outdirAlign
    
  }
  
  # create window data with ZINBA funciton
  tmpFileList <- paste0(outdirWins, "/_tmp.list", collapse = "")
  buildwindowdata(seq = seq,
                  input = inputseq,
                  align = inputDirAlign, 
                  twoBit = twobitfile, 
                  winSize = winSize, 
                  offset = offset, 
                  cnvWinSize = cnvWinSize, 
                  cnvOffset = cnvOffset, 
                  filelist = tmpFileList,
                  filetype = 'bed',
                  extension = extension, 
                  outdir = paste0(outdirWins, "/"), 
                  numProc = numProc)
  
  
  # removing temporary files
  suppressMessages(file.remove(tmpFileList))
  
}



