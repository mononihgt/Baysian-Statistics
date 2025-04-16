# 贝叶斯统计学基础期末作业
# Author: 毛沛炫 3220102692
# Date: 2025-04-15

# Attention!!
# 
# 调试的过程中发现，
# 因为我的程序中没有保存结果的图像，
# 也没有打开图像窗口，
# 如果RStudio的Plots绘图部分太小的话，
# 会报错，画不出图来，
# 所以请先把Plots的部分调整大一点
# 
# 另外，有4张图，可以点击Plots界面左上角的箭头进行切换

graphics.off() 
rm(list=ls())  

cat("--- pre. 开始定义(从utilities搬运)函数，可跳转至Line 506…… ---\n")

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  if (nCIs <= 0) { return(range(sortedPts))}
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

summarizePost = function( paramSampleVec ,
                          compVal=NULL , ROPE=NULL , credMass=0.95 ) {
  meanParam = mean( paramSampleVec )
  medianParam = median( paramSampleVec )
  dres = density( paramSampleVec )
  modeParam = dres$x[which.max(dres$y)]
  
  mcmcEffSz = tryCatch(round( coda::effectiveSize( paramSampleVec ) , 1 ), error=function(e) NA)
  names(mcmcEffSz) = NULL
  hdiLim = HDIofMCMC( paramSampleVec , credMass=credMass )
  if ( !is.null(compVal) ) {
    pcgtCompVal = ( 100 * sum( paramSampleVec > compVal )
                    / length( paramSampleVec ) )
  } else {
    compVal=NA
    pcgtCompVal=NA
  }
  if ( !is.null(ROPE) ) {
    pcltRope = ( 100 * sum( paramSampleVec < ROPE[1] )
                 / length( paramSampleVec ) )
    pcgtRope = ( 100 * sum( paramSampleVec > ROPE[2] )
                 / length( paramSampleVec ) )
    pcinRope = 100-(pcltRope+pcgtRope)
  } else {
    ROPE = c(NA,NA)
    pcltRope=NA
    pcgtRope=NA
    pcinRope=NA
  }
  return( c( Mean=meanParam , Median=medianParam , Mode=modeParam ,
             ESS=mcmcEffSz ,
             HDImass=credMass , HDIlow=hdiLim[1] , HDIhigh=hdiLim[2] ,
             CompVal=compVal , PcntGtCompVal=pcgtCompVal ,
             ROPElow=ROPE[1] , ROPEhigh=ROPE[2] ,
             PcntLtROPE=pcltRope , PcntInROPE=pcinRope , PcntGtROPE=pcgtRope ) )
}

plotPost = function( paramSampleVec , cenTend=c("mode","median","mean")[1] ,
                     compVal=NULL, ROPE=NULL, credMass=0.95, HDItextPlace=0.7,
                     xlab=NULL , xlim=NULL , yaxt=NULL , ylab=NULL ,
                     main=NULL , cex=NULL , cex.lab=NULL ,
                     col=NULL , border=NULL , showCurve=FALSE , breaks=NULL ,
                     ... ) {
  
  if ( is.null(xlab) ) xlab="Param. Val."
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c( compVal , ROPE , paramSampleVec ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="skyblue"
  if ( is.null(border) ) border="white"
  
  if ( inherits(paramSampleVec, "mcmc.list") ) { # Use inherits for robustness
    paramSampleVec = as.matrix(paramSampleVec)
  }
  
  postSummary = summarizePost( paramSampleVec , compVal=compVal , ROPE=ROPE , credMass=credMass )
  
  cvCol = "darkgreen"
  ropeCol = "darkred"
  if ( is.null(breaks) ) {
    if ( diff(range(paramSampleVec)) > 1e-6 ) { # Check range > 0
      HDIrange = postSummary["HDIhigh"] - postSummary["HDIlow"]
      if(is.na(HDIrange)) HDIrange = diff(range(paramSampleVec))/2 # Fallback if HDI fails
      breaks = c( seq( from=min(paramSampleVec) , to=max(paramSampleVec) ,
                       by = HDIrange/18 ) , max(paramSampleVec) )
      breaks = unique(breaks) # Avoid duplicate breaks
    } else {
      breaks=c(paramSampleVec[1]-1.0E-6, paramSampleVec[1]+1.0E-6) # Handle constant value
      if (is.null(border)) border="skyblue" # Ensure border is visible
    }
  }
  
  histinfo <- NULL
  if ( showCurve ) {
    par(xpd=NA)
    densCurve <- tryCatch(density( paramSampleVec , adjust=2 ), warning = function(w) NULL)
    if (!is.null(densCurve)) {
      plot( densCurve$x , densCurve$y , type="l" , lwd=3 , col=col , bty="n" ,
            xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
            main=main , cex=cex , cex.lab=cex.lab , ... )
      histinfo = list(density = densCurve$y, breaks=densCurve$x) # Use density for heights
    } else { # Fallback if density fails
      showCurve <- FALSE # Revert to histogram
      warning("Density curve could not be computed for: ", main)
    }
  }
  if ( !showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks , plot=(!showCurve), ... ) # Only plot if showCurve is FALSE
  }
  
  if (is.null(histinfo) || is.null(histinfo$density)) {
    cenTendHt <- 1; cvHt <- 0.8; ROPEtextHt <- 0.6
  } else {
    maxDens = max(histinfo$density, na.rm=TRUE)
    if(maxDens == 0) maxDens = 1 # Avoid division by zero if density is flat zero
    cenTendHt = 0.9*maxDens
    cvHt = 0.7*maxDens
    ROPEtextHt = 0.55*maxDens
  }
  
  if ( cenTend=="mode" ){ text( postSummary["Mode"] , cenTendHt , bquote(mode==.(signif(postSummary["Mode"],3))) , adj=c(.5,0) , cex=cex ) }
  if ( cenTend=="median" ){ text( postSummary["Median"] , cenTendHt , bquote(median==.(signif(postSummary["Median"],3))) , adj=c(.5,0) , cex=cex ) }
  if ( cenTend=="mean" ){ text( postSummary["Mean"] , cenTendHt , bquote(mean==.(signif(postSummary["Mean"],3))) , adj=c(.5,0) , cex=cex ) }
  
  if ( !is.null( compVal ) && !is.na(compVal) ) {
    pGtCompVal = postSummary["PcntGtCompVal"] / 100
    pLtCompVal = 1 - pGtCompVal
    lines( c(compVal,compVal) , c(0.96*cvHt,0) , lty="dotted" , lwd=2 , col=cvCol )
    text( compVal , cvHt ,
          bquote( .(round(100*pLtCompVal,1)) * "% < " *
                    .(signif(compVal,3)) * " < " *
                    .(round(100*pGtCompVal,1)) * "%" ) ,
          adj=c(pLtCompVal,0) , cex=0.8*cex , col=cvCol )
  }
  if ( !is.null( ROPE ) && !any(is.na(ROPE)) ) {
    pLtROPE = postSummary["PcntLtROPE"]/100
    pInROPE = postSummary["PcntInROPE"]/100
    pGtROPE = postSummary["PcntGtROPE"]/100
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 , col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 , col=ropeCol)
    text( mean(ROPE) , ROPEtextHt ,
          bquote( .(round(100*pLtROPE,1)) * "% < " * .(ROPE[1]) * " < " *
                    .(round(100*pInROPE,1)) * "% < " * .(ROPE[2]) * " < " *
                    .(round(100*pGtROPE,1)) * "%" ) ,
          adj=c(pLtROPE+.5*pInROPE,0) , cex=1 , col=ropeCol ) # cex=1 from original
  }
  lines( c(postSummary["HDIlow"], postSummary["HDIhigh"]) , c(0,0) , lwd=4 , lend=1 )
  text( mean(c(postSummary["HDIlow"], postSummary["HDIhigh"])) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
        adj=c(.5,-1.7) , cex=cex )
  text( postSummary["HDIlow"] , 0 , bquote(.(signif(postSummary["HDIlow"],3))) ,
        adj=c(HDItextPlace,-0.5) , cex=cex )
  text( postSummary["HDIhigh"] , 0 , bquote(.(signif(postSummary["HDIhigh"],3))) ,
        adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  par(xpd=F)
  
  invisible(postSummary) # Return summary silently as in original
}

DbdaAcfPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) {
    stop("parName must be a column name of coda object")
  }
  nChain = nchain(codaObject) # Use nchain()
  if ( is.null(plColors) ) plColors=1:nChain
  if ( nChain > length(plColors) ) plColors = rep(plColors, length.out=nChain) # Handle enough colors
  
  # Initialize matrices/vectors for storing ACF info
  xMatList = vector("list", nChain)
  yMatList = vector("list", nChain)
  maxLag = 0
  
  # Calculate ACF for each chain
  for ( cIdx in 1:nChain ) {
    # Use tryCatch for robustness if ACF fails for a chain
    acfInfo = tryCatch(acf(codaObject[,c(parName)][[cIdx]], plot=FALSE), error=function(e) NULL)
    if (!is.null(acfInfo)) {
      xMatList[[cIdx]] = acfInfo$lag
      yMatList[[cIdx]] = acfInfo$acf
      maxLag = max(maxLag, length(acfInfo$lag))
    } else {
      xMatList[[cIdx]] = numeric(0) # Empty if failed
      yMatList[[cIdx]] = numeric(0)
    }
  }
  
  # Prepare for plotting: create matrices padded with NA
  if (maxLag == 0) { # Handle case where ACF failed for all chains
    plot.new(); title("ACF Plot Failed\n(All Chains)", col.main="red")
    return()
  }
  xMat = matrix(1:maxLag, ncol=1) # Common lag axis
  yMat = matrix(NA, nrow=maxLag, ncol=nChain)
  
  for ( cIdx in 1:nChain ) {
    if(length(yMatList[[cIdx]]) > 0) {
      len = length(yMatList[[cIdx]])
      yMat[1:len, cIdx] = yMatList[[cIdx]]
    }
  }
  
  # Plotting
  matplot( xMat , yMat , type="o" , pch=20 , col=plColors , ylim=c(min(0, min(yMat, na.rm=TRUE)), 1) , # Adjust ylim dynamically
           main="" , xlab="Lag" , ylab="Autocorrelation" )
  abline(h=0,lty="dashed")
  # Calculate ESS safely
  EffChnLngth = tryCatch(effectiveSize(codaObject[,c(parName)]), error=function(e) NA)
  if (!is.na(EffChnLngth)) {
    text( x=maxLag , y=1 , adj=c(1.0,1.0) , cex=1.25 ,
          labels=paste("ESS =",round(EffChnLngth,1)) )
  }
}

DbdaDensPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) {
    stop("parName must be a column name of coda object")
  }
  nChain = nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  if ( nChain > length(plColors) ) plColors = rep(plColors, length.out=nChain)
  
  # Store density info and HDI limits for each chain
  densListX = vector("list", nChain)
  densListY = vector("list", nChain)
  hdiLims = matrix(NA, nrow=2, ncol=nChain)
  maxPoints = 0 # To find the maximum number of points in density estimates
  
  for ( cIdx in 1:nChain ) {
    chainData = codaObject[,c(parName)][[cIdx]]
    # Use tryCatch for density calculation
    densInfo = tryCatch(density(chainData), error=function(e) NULL)
    if (!is.null(densInfo)) {
      densListX[[cIdx]] = densInfo$x
      densListY[[cIdx]] = densInfo$y
      maxPoints = max(maxPoints, length(densInfo$x))
    } else {
      densListX[[cIdx]] = numeric(0)
      densListY[[cIdx]] = numeric(0)
      warning("Density failed for chain ", cIdx, " of parameter ", parName)
    }
    # Calculate HDI even if density fails
    hdiLims[,cIdx] = HDIofMCMC(chainData)
  }
  
  # Prepare matrices for matplot, interpolating if necessary for consistent X grid
  if (maxPoints == 0) { # Handle case where density failed for all chains
    plot.new(); title("Density Plot Failed\n(All Chains)", col.main="red")
    return()
  }
  # Create a common, fine grid of x values covering the overall range
  overallXrange = range(unlist(densListX), na.rm=TRUE)
  if (!all(is.finite(overallXrange))) overallXrange = range(as.matrix(codaObject[,parName])) # Fallback range
  xMatCommon = seq(overallXrange[1], overallXrange[2], length.out=maxPoints)
  yMat = matrix(NA, nrow=maxPoints, ncol=nChain)
  
  for ( cIdx in 1:nChain ) {
    if (length(densListX[[cIdx]]) > 0) {
      # Interpolate density onto the common x grid
      interpY = tryCatch(approx(densListX[[cIdx]], densListY[[cIdx]], xout=xMatCommon, rule=2)$y, error=function(e) NULL)
      if(!is.null(interpY)) yMat[, cIdx] = interpY
    }
  }
  
  # Plotting
  matplot( xMatCommon , yMat , type="l" , col=plColors ,
           main="" , xlab="Param. Value" , ylab="Density" )
  abline(h=0)
  points( hdiLims[1,] , rep(0,nChain) , col=plColors , pch="|" )
  points( hdiLims[2,] , rep(0,nChain) , col=plColors , pch="|" )
  text( mean(hdiLims, na.rm=TRUE) , 0 , "95% HDI" , adj=c(0.5,-0.2) ) # Use mean HDI for text placement
  
  # Calculate and display MCSE
  EffChnLngth = tryCatch(effectiveSize(codaObject[,c(parName)]), error=function(e) NA)
  if (!is.na(EffChnLngth) && EffChnLngth > 0) {
    MCSE = sd(as.matrix(codaObject[,c(parName)]))/sqrt(EffChnLngth)
    maxY = max(yMat, na.rm=TRUE)
    if (is.finite(maxY)) {
      text( max(xMatCommon, na.rm=TRUE) , maxY , adj=c(1.0,1.0) , cex=1.25 ,
            paste("MCSE =\n",signif(MCSE,3)) )
    }
  }
}

diagMCMC = function( codaObject , parName=varnames(codaObject)[1] ) {
  nChain = nchain(codaObject)
  if (nChain == 0) { warning("No chains found in codaObject."); return() }
  
  defaultColors = c("skyblue","black","royalblue","steelblue")
  if (nChain > length(defaultColors)) {
    DBDAplColors = rep(defaultColors, length.out = nChain)
  } else {
    DBDAplColors = defaultColors[1:nChain]
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); if (!identical(oldpar$mfrow, c(1, 1))) layout(1)}, add = TRUE)
  
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) ,
       cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  
  if (!requireNamespace("coda", quietly = TRUE)) { stop("Package 'coda' is required.") }
  
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                   col=DBDAplColors )
  tryVal = try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE ,
                       col=DBDAplColors, # Use assigned colors
                       ask=FALSE),       # Prevent asking
    silent = TRUE # Suppress console error message
  )
  DbdaAcfPlot(codaObject, parName, plColors=DBDAplColors)
  DbdaDensPlot(codaObject, parName, plColors=DBDAplColors)
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
}

genMCMC = function( data , numSavedSteps=20000 , thinSteps=1 ,
                    parameters = c( "theta") ,
                    adaptSteps = 10000 ,
                    burnInSteps = 10000 ,    
                    nChains = 3 ) {
  y = data$y
  s = as.numeric(data$s) # JAGS needs numeric indices
  Ntotal = length(y)
  Nsubj = max(s)
  
  if ( any( y!=0 & y!=1 ) ) { stop("All y values must be 0 or 1.") }
  if ( Nsubj != length(levels(data$s)) ) { # Check if numeric conversion worked as expected
    warning("Subject levels may not be consecutive integers starting from 1.")
  }
  
  dataList = list( y = y, s = s, Ntotal = Ntotal, Nsubj = Nsubj )
  
  # Define the model string
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dbern( theta[s[i]] ) # Likelihood using subject index
    }
    for ( sIdx in 1:Nsubj ) {
      theta[sIdx] ~ dbeta( 2 , 2 ) # Prior for each subject's theta
    }
  }
  "
  writeLines( modelString , con="TEMPmodel.txt" )
  
  initsList = function() {
    thetaInit = rep(0.5, Nsubj) # Default initialization
    for ( sIdx in 1:Nsubj ) {
      includeRows = ( s == sIdx )
      if(sum(includeRows) > 0) { # Check if subject has data
        yThisSubj = y[includeRows]
        prop = mean(yThisSubj)
        thetaInit[sIdx] = max(0.001, min(0.999, prop + runif(1, -0.05, 0.05)))
      } else {
        thetaInit[sIdx] = 0.5 # Fallback if subject has no data
      }
    }
    return( list( theta=thetaInit ) )
  }
  
  
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Iterations per chain
  
  jagsModel = jags.model( "TEMPmodel.txt", data=dataList, inits=initsList,
                          n.chains=nChains, n.adapt=adaptSteps, quiet=TRUE) # Suppress progress bars
  
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps, progress.bar="none" )
  
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel, variable.names=parameters,
                              n.iter=nIter, thin=thinSteps, progress.bar="none" )
  
  if (file.exists("TEMPmodel.txt")) { file.remove("TEMPmodel.txt") }
  
  return( codaSamples )
}

smryMCMC = function( codaSamples , compVal=NULL, rope=NULL,
                     compValDiff=0.0, ropeDiff=NULL, credMass=0.95 ) {
  mcmcMat = as.matrix(codaSamples, chains=TRUE) # Combine chains for summary
  thetaCols = grep("^theta\\[[0-9]+\\]$", colnames(mcmcMat))
  if (length(thetaCols) == 0) stop("No 'theta' parameters found in MCMC output.")
  
  Ntheta = length(thetaCols)
  parameterNames = colnames(mcmcMat)[thetaCols] # Get actual names
  
  summaryInfo = NULL
  
  for ( tIdx in 1:Ntheta ) {
    parName = parameterNames[tIdx]
    summaryInfo = rbind( summaryInfo ,
                         summarizePost( mcmcMat[,parName], compVal=compVal, ROPE=rope, credMass=credMass ))
    rownames(summaryInfo)[nrow(summaryInfo)] = parName # Assign row name after binding
  }
  
  if ( Ntheta >= 2 ) {
    for ( t1Idx in 1:(Ntheta-1) ) {
      for ( t2Idx in (t1Idx+1):Ntheta ) {
        parName1 = parameterNames[t1Idx]
        parName2 = parameterNames[t2Idx]
        diffSample = mcmcMat[,parName1] - mcmcMat[,parName2]
        summaryInfo = rbind( summaryInfo ,
                             summarizePost( diffSample, compVal=compValDiff, ROPE=ropeDiff, credMass=credMass ))
        # Construct difference name carefully
        diffName = paste0(parName1," - ",parName2)
        rownames(summaryInfo)[nrow(summaryInfo)] = diffName # Assign row name
      }
    }
  }
  return( summaryInfo )
}

plotMCMC = function( codaSamples , data , compVal=NULL, rope=NULL,
                     compValDiff=0.0, ropeDiff=NULL, credMass=0.95 ) {
  y = data$y
  s = as.numeric(data$s) 
  
  mcmcMat = as.matrix(codaSamples, chains=TRUE)
  thetaCols = grep("^theta\\[[0-9]+\\]$", colnames(mcmcMat))
  if (length(thetaCols) == 0) stop("No 'theta' parameters found.")
  
  Ntheta = length(thetaCols)
  parameterNames = colnames(mcmcMat)[thetaCols]
  
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); if (!identical(oldpar$mfrow, c(1, 1))) layout(1)}, add = TRUE)
  
  par( mfrow=c(Ntheta,Ntheta) , mar=c(3.5,3.5,1.5,1.0) , mgp=c(2.0,0.7,0) )
  
  for ( t1Idx in 1:Ntheta ) {
    for ( t2Idx in 1:Ntheta ) {
      parName1 = parameterNames[t1Idx]
      parName2 = parameterNames[t2Idx]
      
      if ( t1Idx > t2Idx ) {
        chainLength = NROW( mcmcMat )
        nToPlot = min( 2000, chainLength) # Limit points (as in original example)
        ptIdx = round(seq(1, chainLength, length.out=nToPlot))
        plot( mcmcMat[ptIdx, parName2] , mcmcMat[ptIdx, parName1] , type="p", pch='.', col="skyblue", # Use pch='.' for smaller points
              cex.lab=1.2, xlab=parName2 , ylab=parName1 ) # Slightly smaller labels
      }
      else if ( t1Idx == t2Idx ) {
        includeRows = ( s == t1Idx )
        if(any(includeRows)) {
          yThisSubj = y[includeRows]
          dataPropor = mean(yThisSubj) # Use mean for proportion
        } else { dataPropor = NA }
        
        plotPost( mcmcMat[,parName1], credMass=credMass,
                  compVal=compVal , ROPE=rope ,
                  xlab=parName1, main=parName1, # Use parameter name as title
                  cex.lab=1.2 ) # Match scatter plot label size
        
        if (!is.na(dataPropor)) {
          points( dataPropor , 0 , pch="+" , col="red" , cex=2 ) # Slightly smaller '+'
        }
      }
      else if ( t1Idx < t2Idx ) {
        includeRows1 = ( s == t1Idx )
        includeRows2 = ( s == t2Idx )
        if (any(includeRows1) && any(includeRows2)) {
          dataPropor1 = mean(y[includeRows1])
          dataPropor2 = mean(y[includeRows2])
          dataDiff = dataPropor1 - dataPropor2
        } else { dataDiff = NA }
        
        
        diffSample = mcmcMat[,parName1] - mcmcMat[,parName2]
        diffName = paste0(parName1," - ",parName2) # Difference name
        plotPost( diffSample, credMass=credMass,
                  compVal=compValDiff , ROPE=ropeDiff ,
                  xlab=diffName, main=diffName, # Use difference name as title
                  cex.lab=1.2 )
        
        if (!is.na(dataDiff)) {
          points( dataDiff , 0 , pch="+" , col="red" , cex=2 ) # Slightly smaller '+'
        }
      }
    } 
  } 
}

cat("--- pre. 函数定义(搬运)完成 ---\n")

if (!require(rjags)) { install.packages("rjags"); library(rjags) }
if (!require(coda)) { install.packages("coda"); library(coda) }

set.seed(2692)

cat("--- 1. 生成模拟数据 ---\n")
n_players <- 3
params <- list(
  player1 = list(N_range = 40:50, theta = 0.1),
  player2 = list(N_range = 50:60, theta = 0.4),
  player3 = list(N_range = 60:70, theta = 0.7)
)

y_list <- list(); s_list <- list(); N_total <- 0
actual_N <- numeric(n_players); actual_theta <- numeric(n_players)

for (i in 1:n_players) {
  player_params <- params[[i]]
  N <- sample(player_params$N_range, 1)
  theta <- player_params$theta
  outcomes <- rbinom(n = N, size = 1, prob = theta)
  y_list[[i]] <- outcomes
  s_list[[i]] <- rep(i, N)
  actual_N[i] <- N
  actual_theta[i] <- theta
  N_total <- N_total + N
}

y <- unlist(y_list)
s <- factor(unlist(s_list))
myData <- data.frame(y = y, s = s)

cat("生成的模拟数据如下:\n")
for(i in 1:n_players) {
  cat(sprintf(" 玩家ID %d: N=%d, 理论值theta=%.2f, 模拟数据得到的投篮命中率=%.3f\n",
              i, actual_N[i], actual_theta[i], mean(y_list[[i]])))
}

cat("--- 2. 进行数据分析 ---\n")

ropeDiff = c(-0.05, 0.05) 
parameters = c( "theta")    
adaptSteps = 1000           
burnInSteps = 10000  
numSavedSteps = 10000
nChains = 3                 

mcmcCoda = genMCMC( data=myData , numSavedSteps=numSavedSteps , thinSteps=2 , 
                    adaptSteps=adaptSteps , burnInSteps=burnInSteps ,
                    nChains=3)
cat("MCMC采样完成.\n")

cat("\n--- 3画不同被试的MCMC的图中...\n")
parameterNames = varnames(mcmcCoda) # get all parameter names
thetaNames = parameterNames[grep("^theta", parameterNames)] # Filter for theta parameters

if (length(thetaNames) > 0) {
  for ( parName in thetaNames ) {
    diagMCMC( codaObject=mcmcCoda , parName=parName )
  }
}

cat("\n计算摘要(Summary)信息...\n")
summaryInfo = smryMCMC( mcmcCoda ,
                        compVal=NULL ,    
                        rope=NULL ,       
                        compValDiff=0.0 , 
                        ropeDiff=ropeDiff ) 

cat("\n--- 后验分布的信息输出 ---\n")
cat("输出来自Gelman-Rubin Diagnostic (PSRF)的R hat:\n")
if (nchain(mcmcCoda) > 1) {
  gelmanResults = try(gelman.diag(mcmcCoda), silent=TRUE)
  print(gelmanResults)
}

cat("\n画图并且输出后验分布的统计表述:\n")
print(round(summaryInfo, 3))
plotMCMC( mcmcCoda , data=myData,
          compVal=NULL, rope=NULL, 
          compValDiff=0.0, ropeDiff=ropeDiff ) 
