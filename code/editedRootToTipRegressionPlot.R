
#' Edited from original RootToTipRegressionPlot from treedater package
#' 
#' Plot evolutionary distance from root to sample times and estimated internal node times and regression lines
#' 
#' If a range of sample times was given, these will be estimated. Red and black respectively indicate sample and internal nodes. 
#' This function will print statistics computed from the linear regression model. 
#' 
#' @param td A fitted treedater object 
#' @param show.tip.labels If TRUE, the names of each sample will be plotted at the their corresponding time and evoutionary distance
#' @param textopts An optional list of parameters for plotted tip labels. Passed to the *text* function. 
#' @param pointopts An optional list of parameters for plotted points if showing tip labels. Passed to the *points* function. 
#' @param ... Additional arguments are passed to plot
#' @return The fitted linear model (class 'lm')
#' @examples
#' ## simulate a random tree and sample times for demonstration
#' # make a random tree:
#' tre <- ape::rtree(50)
#' # sample times based on distance from root to tip:
#' sts <- setNames( ape::node.depth.edgelength( tre )[1:ape::Ntip(tre)], tre$tip.label)
#' # modify edge length to represent evolutionary distance with rate 1e-3:
#' tre$edge.length <- tre$edge.length * 1e-3
#' # treedater: 
#' td <- dater( tre, sts =sts, clock='strict', s = 1000, omega0=.0015 )
#' # root to tip regression: 
#' fit = rootToTipRegressionPlot( td )
#' summary(fit)
#' 
#' @export 
#' 
#' 
#' 
#' 


editedRootToTipRegressionPlot <- function(td, show.tip.labels=FALSE, textopts = NULL, pointopts=NULL, ... ){
  
  # stopifnot( inherits( td, 'treedater')) # remove to avoid error when running the function with inputs from saved RData or RDS objects
  
  dT <- ape::node.depth.edgelength( td  )
  dG <- ape::node.depth.edgelength( td$intree )
  #scatter.smooth( dT, dG )
  sts <- (td$timeOfMRCA+dT[1:ape::Ntip(td)])
  nts <- (td$timeOfMRCA+dT)
  mtip  <- lm( dG[1:ape::Ntip(td)] ~ sts )
  mall  <- lm( dG ~ nts )
  smtip <- summary( mtip ) # summary for tips regression: inyongera
#  smmall <- summary( mall ) # inyongera 2
  if ( !show.tip.labels){
    
    graphics::plot( dT + td$timeOfMRCA, dG
                    , col = c(rep('red', ape::Ntip(td)), rep('black', ape::Nnode(td) ) )
                    , xlab = ''
                    , ylab = 'Evolutionary distance'
                    , ...
    )
    
    # INYONGERA
    # define xlim and ylim and coordonates where to add summary outputs in the figure
    do.call( graphics::text, c( list( x = sort(dT)[2]+5 + td$timeOfMRCA, y = max(dG)-0.001 , paste('Root-to-tip mean rate:', round(coef(mtip)[2], digits = 5), '\n') ), textopts ) )
    do.call( graphics::text, c( list( x = sort(dT)[2]+5 + td$timeOfMRCA, y = max(dG)-0.02 , paste('Root-to-tip R square:', round(smtip$r.squared, digits = 5), '\n') ), textopts ) )
    do.call( graphics::text, c( list( x = sort(dT)[2]+5 + td$timeOfMRCA, y = max(dG)-0.01 , paste('Root-to-tip p value:', round(smtip$coefficients[2, 4 ], digits = 5), '\n') ), textopts ) )
    
  } else {

    i <- 1:ape::Ntip(td)
    j <- (ape::Ntip(td)+1):( ape::Ntip(td) + ape::Nnode(td) ) 
    graphics::plot( x = NULL, y = NULL 
                    , xlab = '' 
                    , ylab = 'Evolutionary distance'
                    , xlim = range( dT + td$timeOfMRCA)
                    , ylim = range( dG ) 
                    , ...
    )
    do.call( graphics::points, c( pointopts, list(x =   dT[j] + td$timeOfMRCA, y = dG[j] )))
    # do.call( graphics::text, c( list( x = dT[i] + td$timeOfMRCA, y = dG[i] , labels=td$tip.label ), textopts ) )
    # inyongera
    # define xlim and ylim and coordonates where to add summary outputs in the figure
    do.call( graphics::text, c( list( x = sort(dT)[2]+5 + td$timeOfMRCA, y = max(dG)-0.001 , paste('Root-to-tip mean rate:', round(coef(mtip)[2], digits = 5), '\n') ), textopts ) )
    do.call( graphics::text, c( list( x = sort(dT)[2]+5 + td$timeOfMRCA, y = max(dG)-0.02 , paste('Root-to-tip R square:', round(smtip$r.squared, digits = 5), '\n') ), textopts ) )
    do.call( graphics::text, c( list( x = sort(dT)[2]+5 + td$timeOfMRCA, y = max(dG)-0.01 , paste('Root-to-tip p value:', round(smtip$coefficients[2, 4 ], digits = 5), '\n') ), textopts ) )
  }
  graphics::abline( a = coef(mtip)[1], b = coef(mtip)[2], col = 'red' ) 
  graphics::abline( a = coef(mall)[1], b = coef(mall)[2], col = 'black' ) 
  smtip <- summary( mtip )
  cat(paste( 'Root-to-tip mean rate:', coef(mtip)[2], '\n'))
  cat(paste( 'Root-to-tip p value:', smtip$coefficients[2, 4 ], '\n'))
  cat(paste( 'Root-to-tip R squared (variance explained):', smtip$r.squared, '\n'))
  cat('Returning fitted linear model.\n')
  invisible(mtip)

  
}

