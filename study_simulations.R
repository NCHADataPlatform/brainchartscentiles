# taken from the Lifespan repo, but refactored to return a dataframe

# Some functions must be available
# source this after calibrate_braincharts_local
#DD <- genSimData()
genSimData <- function(thisseed=222) {
  ## COMMON LIBRARIES AND FUNCTIONS
  #  source("100.common-variables.r")
  #  source("101.common-functions.r")
  
  #  source("200.variables.r")
  #  source("201.functions.r")
  
  ## SCRIPT SPECIFIC LIBRARIES
  
  ## SCRIPT SPECIFIC FUNCTIONS
  
  ##
  ## Set random seed for consistent simulations (over-ride any global setting)
  ##
  oldseed <- set.seed(seed=thisseed)
  
  on.exit(set.seed(oldseed))
  ##
  ## This is the omega simulation scenario
  ##
  DATA.TAG <- "omega"
  
  
  ##
  ## Define simulated dataset structure
  ##
  DATA.spec <- data.frame(Study=LETTERS[1:22],
                          N=c(rep(500,10),rep(1000,5),rep(500,3),rep(1000,2),200,1000),
                          novel=c(rep(FALSE,20),rep(TRUE,2)), ## designate two studies for the novel-script
                          r=c(rep(1,15),rep(5,5),1,5),
                          t0=c(seq(0.0,0.8,length=10),    seq(0,0.8,length=5),    seq(0.2,0.8,length=3),    seq(0.2,0.8,length=2),    0.05,0.32),
                          t1=c(seq(0.0,0.8,length=10)+0.1,seq(0,0.8,length=5)+0.1,seq(0.2,0.8,length=3)+0.1,seq(0.2,0.8,length=2)+0.1,0.05,0.18),
                          row.names=LETTERS[1:22]
  )
  ##
  ## Generate study random-effects
  ##
  RE.SD.STUDY <- 0.075
  if( 1 ) {
    ## Define a finite set of normal steps, in terms of standard normal
    RE.SET <- seq(-2.5,2.5,by=0.25)
    
    DATA.spec$MU.ranef.Study <- sample(x=RE.SET * RE.SD.STUDY, ## MUST RESCALE STD-NORMAL TO RE.SD
                                       size=NROW(DATA.spec),replace=TRUE,prob=dnorm(x=RE.SET))
    
    DATA.spec["U","MU.ranef.Study"] <- -0.25 * RE.SD.STUDY ## MUST RESCALE STD-NORMAL TO RE.SD
    DATA.spec["V","MU.ranef.Study"] <-  3.00 * RE.SD.STUDY ## MUST RESCALE STD-NORMAL TO RE.SD
    
  } else {
    ## OR fully random, but less control on simulation effects
    DATA.spec$MU.ranef.Study <- rnorm(n=NROW(DATA.spec),mean=0,sd=RE.SD.STUDY)
  }
  ##
  ## Build subject-level data structure (DATA0)
  ##
  RE.SD.ID <- 0.1
  DATA.parts <- list()
  for( iTYPE in 1:2 ) {
    lCOUNT <- if( iTYPE==1 ) {DATA.spec$N} else {floor(DATA.spec$N/2)} ## implies a 1:2 ratio for CN:non-CN
    DATA.parts[[iTYPE]]          <- data.frame(Study=factor(rep(DATA.spec$Study,times=lCOUNT)))
    DATA.parts[[iTYPE]][,"Grp"]  <- factor(sample(1:2,NROW(DATA.parts[[iTYPE]]),TRUE),1:2,c("F","M"))
    DATA.parts[[iTYPE]][,"Type"] <- factor( rep( iTYPE, NROW(DATA.parts[[iTYPE]]) ), 1:2, c("CN","notCN"))
    
    DATA.parts[[iTYPE]][,"seq"]  <- Reduce(f=c,sapply( rle(as.numeric(DATA.parts[[iTYPE]]$Study))$lengths, function(X){seq(from=1,to=X)} ))
    DATA.parts[[iTYPE]][,"ID"]   <- factor(sprintf("ID%i%04i",iTYPE,DATA.parts[[iTYPE]][,"seq"]))
    
    DATA.parts[[iTYPE]][,"t0.rand"] <- runif(n=NROW(DATA.parts[[iTYPE]]))
    DATA.parts[[iTYPE]][,"t0"]      <- DATA.spec$t0[(DATA.parts[[iTYPE]]$Study)]
    DATA.parts[[iTYPE]][,"t1"]      <- DATA.spec$t1[(DATA.parts[[iTYPE]]$Study)]
    DATA.parts[[iTYPE]][,"r"]       <- DATA.spec$r[(DATA.parts[[iTYPE]]$Study)]
    DATA.parts[[iTYPE]][,"time0"]      <- with( DATA.parts[[iTYPE]], ((t1 - t0)*t0.rand) + t0 )
    
    DATA.parts[[iTYPE]][,"MU.ranef.ID"] <- rnorm(n=NROW(DATA.parts[[iTYPE]]),mean=0,sd=RE.SD.ID)
    DATA.parts[[iTYPE]][,"MU.ranef.Study"] <- DATA.spec$MU.ranef.Study[DATA.parts[[iTYPE]]$Study]
    
    
    DATA.parts[[iTYPE]][,"SIZE"] <- ifelse((iTYPE==1) & (DATA.parts[[iTYPE]][,"Study"]%in%c("U","V")),
                                           DATA.parts[[iTYPE]][,"seq"],
                                           0)
  }
  DATA.A <- Reduce(rbind,DATA.parts)
  attr(DATA.spec,"re.sd")  <- list(Study=RE.SD.STUDY, ID=RE.SD.ID)
  
  ##
  ## Build dataset
  ##
  DATA <- DATA.A[rep(1:NROW(DATA.A),times=DATA.A$r),]
  
  DATA[,"obs"] <- Reduce( f=c, lapply( rle(as.numeric(DATA$ID))$lengths, function(X){1:X} ) ) - 1
  
  DATA[,"time"] <- DATA[,"time0"] + ( 0.04*DATA[,"obs"])
  
  DATA[,"Time"] <- ( 80 * DATA[,"time"] )
  
  #########
  ## modifying the study random effect of timepoint 5 to simulate a different
  ## scanners
  StudyV_TP5 <- which((DATA$Study == "V") & (DATA$obs==4))
  #DATA[StudyV_TP5, "MU.ranef.Study"] <- DATA[StudyV_TP5, "MU.ranef.Study"] + 0.2

  # introduce a SIGMA random effect for timepoint 0
  DATA[,"SIGMA.ranef.Study"] <- 0
  #StudyV_TP4 <- which((DATA$Study == "V") & (DATA$obs==3))
  DATA[StudyV_TP5, "SIGMA.ranef.Study"] <- DATA[StudyV_TP5, "SIGMA.ranef.Study"] + 1
  #StudyV_TP0 <- which((DATA$Study == "V") & (DATA$obs==0))
  #DATA[StudyV_TP0,"SIGMA.ranef.Study"] <- 1
  #########
  TRANSFORMATIONS <- list()
  TRANSFORMATIONS[[ "X" ]] <- list("OriginalName"="Time",
                                   "TransformedName"="TimeTransformed",
                                   "toTransformed"=function(X) { X/10 }, ## must manually scale X-variable for numerical stability within bfpNA()
                                   "toOriginal"=function(X) { 10 * X }
  )
  DATA[,TRANSFORMATIONS[["X"]][["TransformedName"]]] <- TRANSFORMATIONS[["X"]][["toTransformed"]]( DATA[, TRANSFORMATIONS[["X"]][["OriginalName"]] ] )
  
  ##
  ## Generate outcome (including random-effects)
  ##
  if( 1 ) {
    ## Generate Wand
    ##
    Func.mu.fixed.1 <- function( X ) {
      X <- X/8 ## time -> Time -> TimeTransformed transformation
      OUT <- log(-1*(0.4-X)*(0.5-X)+1.8)
      return(OUT)
    }
    Func.mu.fixed.2 <- function( X ) {
      X <- X/8 ## time -> Time -> TimeTransformed transformation
      OUT <- log(-1*(0.35-X)*(0.3-X)+1.55)
      return(OUT)
    }
    
    Func.mu.fixed.3 <- function( X ) {
      X <- X/8 ## time -> Time -> TimeTransformed transformation
      OUT <- 1 + (0.5*X)
      return(OUT)
    }
    Func.mu.fixed.4 <- function( X ) {
      X <- X/8 ## time -> Time -> TimeTransformed transformation
      OUT <- 0.75 + (0.75*X)
      return(OUT)
    }
    
    
    if( 0 ) {
      SEQ <- seq(0,1,length.out=256)
      plot( x=SEQ, y=exp(Func.mu.fixed.1(SEQ)), type="l", col="black", ylim=c(0,2) )
      lines( x=SEQ, y=exp(Func.mu.fixed.2(SEQ)), col="red" )
      
      lines( x=SEQ, y=Func.mu.fixed.3(SEQ), col="black", lty=2 )
      lines( x=SEQ, y=Func.mu.fixed.4(SEQ), col="red", lty=2 )
    }
    
    
    
    attr(DATA.spec,"truth") <- list("Wand"=list(family="GGalt",
                                                MU=list(TypeBase=Func.mu.fixed.1,TypeOther=Func.mu.fixed.2),
                                                SIGMA=list(TypeBase=function(x){log(0.05)},TypeOther=function(x){log(0.05)}),
                                                NU=list(TypeBase=function(x){2},TypeOther=function(x){2})
    ),
    "Wild"=list(family="NO",
                MU=list(TypeBase=Func.mu.fixed.3,TypeOther=Func.mu.fixed.4),
                SIGMA=list(TypeBase=function(x){log(0.05)},TypeOther=function(x){log(0.05)}),
                NU=list(TypeBase=function(x){2},TypeOther=function(x){2})
    )
    
    )
    
    DATA[,"RAND"] <- runif(NROW(DATA),min=0,max=1) ## common random number
    
    
    
    TRUTH.COLUMNS <- c("RAND")
    
    for( LAB in names(attr(DATA.spec,"truth")) ) {
      
      FAMILY <- get(attr(DATA.spec,"truth")[[LAB]]$family)
      
      ARGS.FULL <- list(p=DATA[,"RAND"])
      
      for( lP in names(FAMILY()$parameters) ) {
        
        if( toupper(lP) %in% names(attr(DATA.spec,"truth")[[LAB]]) ) {
          DATA[,sprintf("%s.%s.fixef",LAB,toupper(lP))] <- ifelse(DATA$Type=="CN",
                                                                  attr(DATA.spec,"truth")[[LAB]][[toupper(lP)]]$TypeBase(DATA[,"TimeTransformed"]),
                                                                  attr(DATA.spec,"truth")[[LAB]][[toupper(lP)]]$TypeOther(DATA[,"TimeTransformed"]) )
          TRUTH.COLUMNS <- append(TRUTH.COLUMNS,sprintf("%s.%s.fixef",LAB,toupper(lP)))
          ASPECTS <- c(sprintf("%s.%s.fixef",LAB,toupper(lP)), sprintf("%s.ranef.ID",toupper(lP)), sprintf("%s.ranef.Study",toupper(lP)) )
          DATA[,sprintf("%s.%s",LAB,toupper(lP))] <- rowSums( DATA[ , ASPECTS[ ASPECTS %in% names(DATA) ], drop=FALSE ] )
          
          TRUTH.COLUMNS <- append( TRUTH.COLUMNS, sprintf("%s.%s",LAB,toupper(lP)) )
          
          ARGS.FULL[[lP]] <- FAMILY()[[sprintf("%s.linkinv",lP)]]( DATA[,sprintf("%s.%s",LAB,toupper(lP))] )
        } else {
          stop("Must specify true-functional forms of all gamlss-components!")
        }
      }
      DATA[,LAB] <- do.call( what=get(sprintf("q%s",attr(DATA.spec,"truth")[[LAB]]$family)), args=ARGS.FULL )
    }
  }
  
  ##
  ## Add special columns (INDEX.ID and INDEX.OB) used by later scripts
  ##
  DATA$INDEX.ID <- factor( paste(DATA$Study,DATA$ID,sep="_") )
  DATA$INDEX.OB <- as.integer(DATA$obs+1)
  DATA$INDEX.TYPE <- DATA$Type
  
  return(DATA)
  ##########
  ## From the orginal
  ##
  ## Select only a few columns
  ##
  COLUMNS <- list(Outcomes=names(attr(DATA.spec,"truth")),
                  Covariates=c("Study","Grp","Type","ID","TimeTransformed"),
                  Additional=c("t0.rand", "t0", "t1", "r", "time0","time","Time",
                               TRUTH.COLUMNS),
                  Index=c("INDEX.ID","INDEX.OB","INDEX.TYPE"), ## Would add 'Type' here, to be consistent with real-dataset, but Type is within Covariates
                  Drop=NULL
  )
  
  
  
  
  ##
  ## Add data spec as an attribute
  ##
  attr(DATA,"spec") <- DATA.spec
  
  attr(DATA,"columns") <- COLUMNS
  
  attr(DATA,"tag") <- DATA.TAG
  
  attr(DATA,"Transformations") <- TRANSFORMATIONS
  
  ##
  ## For simulations, we will generate clone and entirely new novel datasets
  ##
  COMMON.SET <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", 
                  "M", "N", "O", "P", "Q", "R", "S", "T")
  SIZE.SEQ <- c(5,10,15,20,25,30,35,40,45,50,100,150,200)
  attr(DATA,"SCENARIOS") <- list(list(label="__",base=COMMON.SET,added=NULL,sizes=c(0)),
                                 list(label="u_",base=COMMON.SET,added=c("U"),sizes=SIZE.SEQ),
                                 list(label="_v",base=COMMON.SET,added=c("V"),sizes=SIZE.SEQ),
                                 list(label="uv",base=COMMON.SET,added=c("U","V"),sizes=SIZE.SEQ)
  )
  
  
  ##
  ## Save dataset in RDS format for use in later scripts
  ##
  
  DATA.PATH <- file.path( RDS.DIR, DATA.TAG )
  if( !dir.exists( DATA.PATH ) ) {
    dir.create( DATA.PATH, recursive=TRUE)
  }
  TOSAVE <- DATA[ , unlist(attr(DATA,"columns")) ]
  
  attributes(TOSAVE) <- c( attributes(TOSAVE), attributes(DATA)[c("columns","tag","Transformations","spec","SCENARIOS")] )
  
  Check.Attributes(TOSAVE)
  
  return(TOSAVE)  
  
}

# generates data without the site effect offset
genSimDataNoSiteEffect <- function(thisseed=222) {


  ## COMMON LIBRARIES AND FUNCTIONS
  #  source("100.common-variables.r")
  #  source("101.common-functions.r")
  
  #  source("200.variables.r")
  #  source("201.functions.r")
  
  ## SCRIPT SPECIFIC LIBRARIES
  
  ## SCRIPT SPECIFIC FUNCTIONS
  
  ##
  ## Set random seed for consistent simulations (over-ride any global setting)
  ##
  oldseed <- set.seed(seed=thisseed)
  
  on.exit(set.seed(oldseed))
  ##
  ## This is the omega simulation scenario
  ##
  DATA.TAG <- "omega"
  
  
  ##
  ## Define simulated dataset structure
  ##
  DATA.spec <- data.frame(Study=LETTERS[1:22],
                          N=c(rep(500,10),rep(1000,5),rep(500,3),rep(1000,2),200,1000),
                          novel=c(rep(FALSE,20),rep(TRUE,2)), ## designate two studies for the novel-script
                          r=c(rep(1,15),rep(5,5),1,5),
                          t0=c(seq(0.0,0.8,length=10),    seq(0,0.8,length=5),    seq(0.2,0.8,length=3),    seq(0.2,0.8,length=2),    0.05,0.32),
                          t1=c(seq(0.0,0.8,length=10)+0.1,seq(0,0.8,length=5)+0.1,seq(0.2,0.8,length=3)+0.1,seq(0.2,0.8,length=2)+0.1,0.05,0.18),
                          row.names=LETTERS[1:22]
  )
  ##
  ## Generate study random-effects
  ##
  RE.SD.STUDY <- 0.075
  if( 1 ) {
    ## Define a finite set of normal steps, in terms of standard normal
    RE.SET <- seq(-2.5,2.5,by=0.25)
    
    DATA.spec$MU.ranef.Study <- sample(x=RE.SET * RE.SD.STUDY, ## MUST RESCALE STD-NORMAL TO RE.SD
                                       size=NROW(DATA.spec),replace=TRUE,prob=dnorm(x=RE.SET))
    
    DATA.spec["U","MU.ranef.Study"] <- -0.25 * RE.SD.STUDY ## MUST RESCALE STD-NORMAL TO RE.SD
    DATA.spec["V","MU.ranef.Study"] <-  3.00 * RE.SD.STUDY ## MUST RESCALE STD-NORMAL TO RE.SD
    
  } else {
    ## OR fully random, but less control on simulation effects
    DATA.spec$MU.ranef.Study <- rnorm(n=NROW(DATA.spec),mean=0,sd=RE.SD.STUDY)
  }
  ##
  ## Build subject-level data structure (DATA0)
  ##
  RE.SD.ID <- 0.1
  DATA.parts <- list()
  for( iTYPE in 1:2 ) {
    lCOUNT <- if( iTYPE==1 ) {DATA.spec$N} else {floor(DATA.spec$N/2)} ## implies a 1:2 ratio for CN:non-CN
    DATA.parts[[iTYPE]]          <- data.frame(Study=factor(rep(DATA.spec$Study,times=lCOUNT)))
    DATA.parts[[iTYPE]][,"Grp"]  <- factor(sample(1:2,NROW(DATA.parts[[iTYPE]]),TRUE),1:2,c("F","M"))
    DATA.parts[[iTYPE]][,"Type"] <- factor( rep( iTYPE, NROW(DATA.parts[[iTYPE]]) ), 1:2, c("CN","notCN"))
    
    DATA.parts[[iTYPE]][,"seq"]  <- Reduce(f=c,sapply( rle(as.numeric(DATA.parts[[iTYPE]]$Study))$lengths, function(X){seq(from=1,to=X)} ))
    DATA.parts[[iTYPE]][,"ID"]   <- factor(sprintf("ID%i%04i",iTYPE,DATA.parts[[iTYPE]][,"seq"]))
    
    DATA.parts[[iTYPE]][,"t0.rand"] <- runif(n=NROW(DATA.parts[[iTYPE]]))
    DATA.parts[[iTYPE]][,"t0"]      <- DATA.spec$t0[(DATA.parts[[iTYPE]]$Study)]
    DATA.parts[[iTYPE]][,"t1"]      <- DATA.spec$t1[(DATA.parts[[iTYPE]]$Study)]
    DATA.parts[[iTYPE]][,"r"]       <- DATA.spec$r[(DATA.parts[[iTYPE]]$Study)]
    DATA.parts[[iTYPE]][,"time0"]      <- with( DATA.parts[[iTYPE]], ((t1 - t0)*t0.rand) + t0 )
    
    DATA.parts[[iTYPE]][,"MU.ranef.ID"] <- rnorm(n=NROW(DATA.parts[[iTYPE]]),mean=0,sd=RE.SD.ID)
    DATA.parts[[iTYPE]][,"MU.ranef.Study"] <- DATA.spec$MU.ranef.Study[DATA.parts[[iTYPE]]$Study]
    
    
    DATA.parts[[iTYPE]][,"SIZE"] <- ifelse((iTYPE==1) & (DATA.parts[[iTYPE]][,"Study"]%in%c("U","V")),
                                           DATA.parts[[iTYPE]][,"seq"],
                                           0)
  }
  DATA.A <- Reduce(rbind,DATA.parts)
  attr(DATA.spec,"re.sd")  <- list(Study=RE.SD.STUDY, ID=RE.SD.ID)
  
  ##
  ## Build dataset
  ##
  DATA <- DATA.A[rep(1:NROW(DATA.A),times=DATA.A$r),]
  
  DATA[,"obs"] <- Reduce( f=c, lapply( rle(as.numeric(DATA$ID))$lengths, function(X){1:X} ) ) - 1
  
  DATA[,"time"] <- DATA[,"time0"] + ( 0.04*DATA[,"obs"])
  
  DATA[,"Time"] <- ( 80 * DATA[,"time"] )
  
  #########
  ## RB modifying the study random effect of timepoint 5 to simulate a different
  ## scanners
  #StudyV_TP5 <- which((DATA$Study == "V") & (DATA$obs==4))
  #DATA[StudyV_TP5, "MU.ranef.Study"] <- DATA[StudyV_TP5, "MU.ranef.Study"] + 0.2
  #########
  TRANSFORMATIONS <- list()
  TRANSFORMATIONS[[ "X" ]] <- list("OriginalName"="Time",
                                   "TransformedName"="TimeTransformed",
                                   "toTransformed"=function(X) { X/10 }, ## must manually scale X-variable for numerical stability within bfpNA()
                                   "toOriginal"=function(X) { 10 * X }
  )
  DATA[,TRANSFORMATIONS[["X"]][["TransformedName"]]] <- TRANSFORMATIONS[["X"]][["toTransformed"]]( DATA[, TRANSFORMATIONS[["X"]][["OriginalName"]] ] )
  
  ##
  ## Generate outcome (including random-effects)
  ##
  if( 1 ) {
    ## Generate Wand
    ##
    Func.mu.fixed.1 <- function( X ) {
      X <- X/8 ## time -> Time -> TimeTransformed transformation
      OUT <- log(-1*(0.4-X)*(0.5-X)+1.8)
      return(OUT)
    }
    Func.mu.fixed.2 <- function( X ) {
      X <- X/8 ## time -> Time -> TimeTransformed transformation
      OUT <- log(-1*(0.35-X)*(0.3-X)+1.55)
      return(OUT)
    }
    
    Func.mu.fixed.3 <- function( X ) {
      X <- X/8 ## time -> Time -> TimeTransformed transformation
      OUT <- 1 + (0.5*X)
      return(OUT)
    }
    Func.mu.fixed.4 <- function( X ) {
      X <- X/8 ## time -> Time -> TimeTransformed transformation
      OUT <- 0.75 + (0.75*X)
      return(OUT)
    }
    
    
    if( 0 ) {
      SEQ <- seq(0,1,length.out=256)
      plot( x=SEQ, y=exp(Func.mu.fixed.1(SEQ)), type="l", col="black", ylim=c(0,2) )
      lines( x=SEQ, y=exp(Func.mu.fixed.2(SEQ)), col="red" )
      
      lines( x=SEQ, y=Func.mu.fixed.3(SEQ), col="black", lty=2 )
      lines( x=SEQ, y=Func.mu.fixed.4(SEQ), col="red", lty=2 )
    }
    
    
    
    attr(DATA.spec,"truth") <- list("Wand"=list(family="GGalt",
                                                MU=list(TypeBase=Func.mu.fixed.1,TypeOther=Func.mu.fixed.2),
                                                SIGMA=list(TypeBase=function(x){log(0.05)},TypeOther=function(x){log(0.05)}),
                                                NU=list(TypeBase=function(x){2},TypeOther=function(x){2})
    ),
    "Wild"=list(family="NO",
                MU=list(TypeBase=Func.mu.fixed.3,TypeOther=Func.mu.fixed.4),
                SIGMA=list(TypeBase=function(x){log(0.05)},TypeOther=function(x){log(0.05)}),
                NU=list(TypeBase=function(x){2},TypeOther=function(x){2})
    )
    )
    
    DATA[,"RAND"] <- runif(NROW(DATA),min=0,max=1) ## common random number
    
    
    
    TRUTH.COLUMNS <- c("RAND")
    
    for( LAB in names(attr(DATA.spec,"truth")) ) {
      
      FAMILY <- get(attr(DATA.spec,"truth")[[LAB]]$family)
      
      ARGS.FULL <- list(p=DATA[,"RAND"])
      
      for( lP in names(FAMILY()$parameters) ) {
        
        if( toupper(lP) %in% names(attr(DATA.spec,"truth")[[LAB]]) ) {
          DATA[,sprintf("%s.%s.fixef",LAB,toupper(lP))] <- ifelse(DATA$Type=="CN",
                                                                  attr(DATA.spec,"truth")[[LAB]][[toupper(lP)]]$TypeBase(DATA[,"TimeTransformed"]),
                                                                  attr(DATA.spec,"truth")[[LAB]][[toupper(lP)]]$TypeOther(DATA[,"TimeTransformed"]) )
          TRUTH.COLUMNS <- append(TRUTH.COLUMNS,sprintf("%s.%s.fixef",LAB,toupper(lP)))
          ASPECTS <- c(sprintf("%s.%s.fixef",LAB,toupper(lP)), sprintf("%s.ranef.ID",toupper(lP)), sprintf("%s.ranef.Study",toupper(lP)) )
          
          DATA[,sprintf("%s.%s",LAB,toupper(lP))] <- rowSums( DATA[ , ASPECTS[ ASPECTS %in% names(DATA) ], drop=FALSE ] )
          
          TRUTH.COLUMNS <- append( TRUTH.COLUMNS, sprintf("%s.%s",LAB,toupper(lP)) )
          
          ARGS.FULL[[lP]] <- FAMILY()[[sprintf("%s.linkinv",lP)]]( DATA[,sprintf("%s.%s",LAB,toupper(lP))] )
          
        } else {
          stop("Must specify true-functional forms of all gamlss-components!")
        }
      }
      #print(sprintf("q%s",attr(DATA.spec,"truth")[[LAB]]$family))
      # generating random values from the distribution
      DATA[,LAB] <- do.call( what=get(sprintf("q%s",attr(DATA.spec,"truth")[[LAB]]$family)), args=ARGS.FULL )
      browser()
    }
  }
  
  ##
  ## Add special columns (INDEX.ID and INDEX.OB) used by later scripts
  ##
  DATA$INDEX.ID <- factor( paste(DATA$Study,DATA$ID,sep="_") )
  DATA$INDEX.OB <- as.integer(DATA$obs+1)
  DATA$INDEX.TYPE <- DATA$Type
  
  return(DATA)
  ##########
  ## From the orginal
  ##
  ## Select only a few columns
  ##
  COLUMNS <- list(Outcomes=names(attr(DATA.spec,"truth")),
                  Covariates=c("Study","Grp","Type","ID","TimeTransformed"),
                  Additional=c("t0.rand", "t0", "t1", "r", "time0","time","Time",
                               TRUTH.COLUMNS),
                  Index=c("INDEX.ID","INDEX.OB","INDEX.TYPE"), ## Would add 'Type' here, to be consistent with real-dataset, but Type is within Covariates
                  Drop=NULL
  )
  
  
  
  
  ##
  ## Add data spec as an attribute
  ##
  attr(DATA,"spec") <- DATA.spec
  
  attr(DATA,"columns") <- COLUMNS
  
  attr(DATA,"tag") <- DATA.TAG
  
  attr(DATA,"Transformations") <- TRANSFORMATIONS
  
  ##
  ## For simulations, we will generate clone and entirely new novel datasets
  ##
  COMMON.SET <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", 
                  "M", "N", "O", "P", "Q", "R", "S", "T")
  SIZE.SEQ <- c(5,10,15,20,25,30,35,40,45,50,100,150,200)
  attr(DATA,"SCENARIOS") <- list(list(label="__",base=COMMON.SET,added=NULL,sizes=c(0)),
                                 list(label="u_",base=COMMON.SET,added=c("U"),sizes=SIZE.SEQ),
                                 list(label="_v",base=COMMON.SET,added=c("V"),sizes=SIZE.SEQ),
                                 list(label="uv",base=COMMON.SET,added=c("U","V"),sizes=SIZE.SEQ)
  )
  
  
  ##
  ## Save dataset in RDS format for use in later scripts
  ##
  
  DATA.PATH <- file.path( RDS.DIR, DATA.TAG )
  if( !dir.exists( DATA.PATH ) ) {
    dir.create( DATA.PATH, recursive=TRUE)
  }
  TOSAVE <- DATA[ , unlist(attr(DATA,"columns")) ]
  
  attributes(TOSAVE) <- c( attributes(TOSAVE), attributes(DATA)[c("columns","tag","Transformations","spec","SCENARIOS")] )
  
  Check.Attributes(TOSAVE)
  
  return(TOSAVE)  
  
}