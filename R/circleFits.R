#' Fit a circle to a point cloud
#'
#' This function fits a circle to a point cloud using the least squares method.
#' The function returns the estimated parameters of the circle, the root mean square error (RMSE) of the distances from the points to the circle, and the x and y coordinates of the center of the circle.
#'
#' @param xp A vector of the x coordinates of the point cloud.
#' @param yp A vector of the y coordinates of the point cloud.
#' @param c0 If TRUE, centers x and y to zero.
#'
#' @return A vector with four elements: The x and y coordinates of the center of the circle, the estimated radius of the circle, and the root mean square error (rmse) of the distances from the points to the circle.
#'
circlefit = function (xp, yp, c0=T){

  if (!is.vector(xp, mode = "numeric") || !is.vector(yp, mode = "numeric"))
    stop("Arguments 'xp' and 'yp' must be numeric vectors.")
  if (length(xp) != length(yp))
    stop("Vectors 'xp' and 'yp' must be of the same length.")
  if (any(is.na(xp)) || any(is.na(yp)))
    stop("Vectors 'xp' and 'yp' must not contain NA values.")
  if (length(xp) < 3)
    stop("At least three points are required to fit a circle.")

  if(c0){
    cen = c(x=mean(xp), y=mean(yp))
    xp = xp-cen[1]
    yp = yp-cen[2]
  }

  n <- length(xp)
  p <- qr.solve(cbind(xp, yp, 1), matrix(xp^2 + yp^2, ncol = 1))
  r <- c(p[1]/2, p[2]/2, sqrt((p[1]^2 + p[2]^2)/4 + p[3]))
  cerr <- function(v) {
    dist <- sqrt((xp - v[1])^2 + (yp - v[2])^2)
    vapply(dist, function(x) (x - v[3])^2, FUN.VALUE = numeric(1))
  }

  q <- optim(p, function(v) sum(cerr(v)))
  r <- q$par
  rmse <- sqrt(q$value/n)

  out = c(unlist(r), rmse)
  if(c0) out[1:2] = out[1:2]+cen

  return(out)
}

#' Fit a circle to a point cloud
#'
#' This function fits a circle to a point cloud using the least squares method.
#' The function returns the estimated parameters of the circle, the root mean square error (RMSE) of the distances from the points to the circle, and the x and y coordinates of the center of the circle.
#'
#' @param xp A vector of the x coordinates of the point cloud.
#' @param yp A vector of the y coordinates of the point cloud.
#' @param c0 If TRUE, centers x and y to zero.
#'
#' @return A vector with four elements: The x and y coordinates of the center of the circle, the estimated radius of the circle, and the root mean square error (rmse) of the distances from the points to the circle.
#'
#' @examples
#'  \dontrun{
#' # Create a set of points on a circle with center (2, 3) and radius 5
#' theta <- seq(0, 2*pi, length.out = 100)
#' xp <- 2 + 5 * cos(theta)
#' yp <- 3 + 5 * sin(theta)
#'
#' # Use circlefit to estimate the center, radius and RMSE
#' result <- circlefit(xp, yp)
#'}
#' @export
least.sq.circle = function (xp, yp, c0=T){

  if (!is.vector(xp, mode = "numeric") || !is.vector(yp, mode = "numeric"))
    stop("Arguments 'xp' and 'yp' must be numeric vectors.")
  if (length(xp) != length(yp))
    stop("Vectors 'xp' and 'yp' must be of the same length.")
  if (any(is.na(xp)) || any(is.na(yp)))
    stop("Vectors 'xp' and 'yp' must not contain NA values.")
  if (length(xp) < 3)
    stop("At least three points are required to fit a circle.")

  if(c0){
    cen = c(x=mean(xp), y=mean(yp))
    xp = xp-cen[1]
    yp = yp-cen[2]
  }

  n <- length(xp)
  p <- qr.solve(cbind(xp, yp, 1), matrix(xp^2 + yp^2, ncol = 1))
  r <- c(p[1]/2, p[2]/2, sqrt((p[1]^2 + p[2]^2)/4 + p[3]))
  cerr <- function(v) {
    dist <- sqrt((xp - v[1])^2 + (yp - v[2])^2)
    vapply(dist, function(x) (x - v[3])^2, FUN.VALUE = numeric(1))
  }

  q <- optim(p, function(v) sum(cerr(v)))
  r <- q$par
  rmse <- sqrt(q$value/n)

  out = c(unlist(r), rmse)
  if(c0) out[1:2] = out[1:2]+cen

  return(out)
}

#' RAndom SAmple Consensus (RANSAC) circle fit
#'
#' This function fits a circle to a point cloud using the RANSAC algorithm.
#' The function returns the estimated parameters of the circle, the root mean square error (RMSE) of the distances from the points to the circle, and the x and y coordinates of the center of the circle.
#'
#' @param slice A matrix of the point cloud coordinates.
#' @param n The number of points to sample on every RANSAC iteration.
#' @param p The estimated proportion of inliers in the dataset.
#' @param P The level of confidence desired.
#'
#' @return A vector with four elements: The x and y coordinates of the center of the circle, the estimated radius of the circle, and the root mean square error (rmse) of the distances from the points to the circle.
#'
#' @examples
#'  \dontrun{
#' # Create a set of points on a circle with center (2, 3) and radius 5
#' theta <- seq(0, 2*pi, length.out = 100)
#' xp <- 2 + 5 * cos(theta)
#' yp <- 3 + 5 * sin(theta)
#' slice <- cbind(xp, yp)
#'
#' # Use RANSAC.circle to estimate the center, radius and RMSE
#' result <- RANSAC.circle(slice, n=15, p=.8, P=.99)
#'}
#' @export
ransac.circle = function(slice, n=15, p=.8, P=.99){

  if (!is.matrix(slice) || ncol(slice) < 2)
    stop("Argument 'slice' must be a matrix with at least two columns.")
  if (!is.numeric(n) || length(n) != 1 || n <= 0)
    stop("Argument 'n' must be a positive number.")
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p > 1)
    stop("Argument 'p' must be a number between 0 and 1.")
  if (!is.numeric(P) || length(P) != 1 || P <= 0 || P > 1)
    stop("Argument 'P' must be a number between 0 and 1.")

  if(nrow(slice) < n) n = nrow(slice)

  N = max(1, round(log(1 - P) / log(1 - p^n)))

  data = matrix(ncol=4, nrow=N)
  for(j in seq_len(N)){
    a = sample(nrow(slice), size = n)

    b = tryCatch(circlefit(slice[a,1], slice[a,2]),
                 error = function(con){ return(NULL) },
                 warning = function(con) return(NULL))

    if(is.null(b)) next

    data[j,] = b
  }

  data = data[!is.na(data[,1]),]

  if(nrow(data) == 0){ dt = NULL }else{
    c = which.min(data[,4])
    dt = data[c,]
  }

  return(dt)
}

#' Fit a circle to a point cloud
#'
#' This function fits a circle to a point cloud using either the RANSAC or least squares method.
#' The function returns the estimated parameters of the circle, the root mean square error (RMSE) of the distances from the points to the circle, and the x, y coordinates of the center of the circle.
#'
#' @param pts_mat A matrix of the point cloud coordinates.
#' @param method The method to use for fitting. Either 'ransac' or 'leastsq'.
#' @param n The number of points sampled in every RANSAC iteration. Only used if method is 'ransac'.
#' @param p The proportion of inliers in P. Only used if method is 'ransac'.
#' @param P The level of confidence desired. Only used if method is 'ransac'.
#' @param c0 If TRUE, centers x and y to zero. Only used if method is 'leastsq'.
#'
#' @return A vector with four elements: The x and y coordinates of the center of the circle, the estimated radius of the circle, and the root mean square error (rmse) of the distances from the points to the circle.
#'
#' @examples
#'  \dontrun{
#' # Create a set of points on a circle with center (2, 3) and radius 5
#' theta <- seq(0, 2*pi, length.out = 100)
#' xp <- 2 + 5 * cos(theta)
#' yp <- 3 + 5 * sin(theta)
#' pts_mat <- cbind(xp, yp)
#'
#' # Use fit_circle to estimate the center, radius and RMSE
#' result <- fit_circle(pts_mat, method = 'ransac', n=15, p=.8, P=.99)
#'}
#' @export
fit_circle = function(pts_mat, method = "ransac", n=20, p=.9, P=.99, c0=T){
  if(method == "ransac"){
    result = ransac.circle(pts_mat, n, p, P)
  } else if(method == "leastsq"){
    result = least.sq.circle(pts_mat[,1], pts_mat[,2], c0)
  } else {
    stop("Invalid method. Choose either 'ransac' or 'leastsq'.")
  }

  return(result)
}
