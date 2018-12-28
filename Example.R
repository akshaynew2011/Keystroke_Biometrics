
library( MASS );
library( ROCR );
library( stats );

euclideanTrain <- function( YTrain ) {
  dmod <- list ( mean  = colMeans( YTrain ) );
  return( dmod );
}

euclideanScore <- function( dmod, YScore ) {
  p <- length( dmod$mean );
  n <- nrow( YScore );
  
  if( ncol(YScore) != p ) stop("Training/test feature length mismatch ");
  
  meanMatrix <- matrix( dmod$mean, byrow=TRUE, nrow=n, ncol=p );
  
  scores <- rowSums( ( YScore - meanMatrix )^2 );
  
  return( scores );
}


manhattanTrain <- function( YTrain ) {
  dmod <- list ( mean  = colMeans( YTrain ) );
  return( dmod );
}

manhattanScore <- function( dmod, YScore ) {
  p <- length( dmod$mean );
  n <- nrow( YScore );
  
  if( ncol(YScore) != p ) stop("Training/test feature length mismatch ");
  
  meanMatrix <- matrix( dmod$mean, byrow=TRUE, nrow=n, ncol=p );
  
  scores <- rowSums( abs( YScore - meanMatrix ) );
  
  return( scores );
}


mahalanobisTrain <- function( YTrain ) {
  dmod <- list( mean  = colMeans( YTrain ),
                covInv = ginv( cov( YTrain ) ) );
  return( dmod );
}

mahalanobisScore <- function( dmod, YScore ) {
  p <- length( dmod$mean );
  n <- nrow( YScore );
  
  if( ncol(YScore) != p ) stop("Training/test feature length mismatch ");
  
  scores <- mahalanobis( YScore, dmod$mean, dmod$covInv, inverted=TRUE );
  
  return( scores );
}

# The detectorSet data structure compiles the train and score
# functions for each detector into a named list, with the name
# corresponding to the name of the detector.

detectorSet = list(
  Euclidean =
    list( train = euclideanTrain,
          score = euclideanScore ),
  Manhattan =
    list( train = manhattanTrain,
          score = manhattanScore ),
  Mahalanobis =
    list( train = mahalanobisTrain,
          score = mahalanobisScore ) );


calculateEqualError <- function( userScores, impostorScores ) {
  
  predictions <- c( userScores, impostorScores );
  labels <- c( rep( 0, length( userScores ) ),
               rep( 1, length( impostorScores ) ) );
  
  pred <- prediction(predictions, labels)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
  plot(perf, col=rainbow(10))
  
  
  pred <- prediction( predictions, labels );
  
  missrates <- pred@fn[[1]] / pred@n.pos[[1]];
  farates <- pred@fp[[1]] / pred@n.neg[[1]];
  
  dists <- missrates - farates;
  idx1 <- which( dists == min( dists[ dists >= 0 ] ) );
  idx2 <- which( dists == max( dists[ dists < 0 ] ) );
  stopifnot( length( idx1 ) == 1 );
  stopifnot( length( idx2 ) == 1 );
  stopifnot( abs( idx1 - idx2 ) == 1 );
  
  
  x <- c( missrates[idx1], farates[idx1] );
  y <- c( missrates[idx2], farates[idx2] );
  a <- ( x[1] - x[2] ) / ( y[2] - x[2] - y[1] + x[1] );
  eer <- x[1] + a * ( y[1] - x[1] );
  
  return( eer );
}




evaluateSubject <- function( X, evalSubject, detectorTrain, detectorScore ) {
  
  
  YTrain <- as.matrix( subset( X,
                               subset = ( subject == evalSubject &
                                            sessionIndex <= 4 ),
                               select = -c( subject, sessionIndex, rep ) ) );
  
  YScore0 <- as.matrix( subset( X,
                                subset = ( subject == evalSubject &
                                             sessionIndex > 4 ),
                                select = -c( subject, sessionIndex, rep ) ) );
  
  YScore1 <- as.matrix( subset( X,
                                subset = ( subject != evalSubject &
                                             sessionIndex == 1 &
                                             rep <= 5 ),
                                select = -c( subject, sessionIndex, rep ) ) );
  
  
  dmod <- detectorTrain( YTrain );
  userScores <- detectorScore( dmod, YScore0 );
  impostorScores <- detectorScore( dmod, YScore1 );
  
  eer <- calculateEqualError( userScores, impostorScores );
  
  return( eer );  
}

cat("Loading the data file\n");

datafile <- 'DSL-StrongPasswordData.txt';
if( ! file.exists(datafile) )
  stop( "Password data file ",datafile," does not exist");

X <- read.table( datafile, header = TRUE );
subjects <- sort( levels( X$subject ) );

eers <- list();
for( detectorName in names( detectorSet ) ) {
  
  
  cat("Evaluating the",detectorName,"detector\n");
  detectorTrain = detectorSet[[ detectorName ]]$train;
  detectorScore = detectorSet[[ detectorName ]]$score;
  
  eers[[ detectorName ]] <- rep( NA, length(subjects) );
  
  n <- length(subjects);
  for( i in 1:n ) {
    
    eer <- evaluateSubject( X, subjects[i],
                            detectorTrain = detectorTrain,
                            detectorScore = detectorScore );
    
    eers[[ detectorName ]][i] <- eer;
    cat("\r  ",i,"/",n,":",eer);
  }
  cat("\r  average equal-error:",mean(eers[[detectorName]]),"\n");
}


cat("Tabulating results:\n");

eers <- data.frame( eers );
rownames( eers ) <- subjects;

res <- data.frame(eer.mean = colMeans(eers),
                  eer.sd   = apply( eers, 2, sd ));



