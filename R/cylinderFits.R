# ------------------------------------------------------------------------------
# 'create_cylinder' function
#
# Creates a cylinder and generates points within it.
#
# @param diameter Numeric. Diameter of the cylinder.
# @param length Numeric. Length (height) of the cylinder.
# @param orientation_angle Numeric. Angle to rotate the cylinder around the z-axis.
# @param n_points Integer. Number of points to be generated within the cylinder.
# @param sd_points Numeric. Standard deviation for generating the points.
#
# @return A data frame containing the coordinates of the points generated within the
#   cylinder. The data frame has three columns: X, Y, and Z, which represent
#   the x, y, and z coordinates of the points, respectively.
# ------------------------------------------------------------------------------
create_cylinder <- function(diameter, length, orientation_angle, n_points,
                            sd_points) {
  # Adjust the orientation angle
  adjusted_orientation_angle <- orientation_angle - 90
  if (adjusted_orientation_angle < 0) {
    adjusted_orientation_angle <- adjusted_orientation_angle + 360
  }
  # Generate cylinder coordinates
  x <- c(0, diameter / 2, diameter / 2, 0, 0,-diameter / 2,-diameter / 2, 0, 0)
  y <- c(0, 0, length, length, 0, 0, length, length, 0)
  z <- rep(0, 9)
  # Rotate cylinder around z-axis
  angle_z <- adjusted_orientation_angle
  x_rot_z <- x * cos(angle_z) - y * sin(angle_z)
  y_rot_z <- x * sin(angle_z) + y * cos(angle_z)
  z_rot_z <- z
  # Rotate cylinder around y-axis
  angle_y <- runif(1, 0, 2 * pi)
  x_rot_y <- x_rot_z * cos(angle_y) + z_rot_z * sin(angle_y)
  y_rot_y <- y_rot_z
  z_rot_y <- -x_rot_z * sin(angle_y) + z_rot_z * cos(angle_y)
  # Rotate cylinder around x-axis
  angle_x <- runif(1, 0, 2 * pi)
  x_rot_x <- x_rot_y
  y_rot_x <- y_rot_y * cos(angle_x) - z_rot_y * sin(angle_x)
  z_rot_x <- y_rot_y * sin(angle_x) + z_rot_y * cos(angle_x)
  # Generate points within cylinder using cylindrical coordinates
  theta <- runif(n_points, 0, 2 * pi)
  r <- rnorm(n_points, diameter / 2, sd_points)
  h <- rnorm(n_points, 0, length)
  x_points <- r * cos(theta)
  y_points <- r * sin(theta)
  z_points <- h
  # Rotate points around x-axis
  x_points_rot_x <- x_points
  y_points_rot_x <- y_points * cos(angle_x) + z_points * sin(angle_x)
  z_points_rot_x <- -y_points * sin(angle_x) + z_points * cos(angle_x)
  # Rotate points around y-axis
  x_points_rot_y <-
    x_points_rot_x * cos(angle_y) + z_points_rot_x * sin(angle_y)
  y_points_rot_y <- y_points_rot_x
  z_points_rot_y <-
    -x_points_rot_x * sin(angle_y) + z_points_rot_x * cos(angle_y)
  # Rotate points around z-axis
  x_points_rot_z <-
    x_points_rot_y * cos(angle_z) - y_points_rot_y * sin(angle_z)
  y_points_rot_z <-
    x_points_rot_y * sin(angle_z) + y_points_rot_y * cos(angle_z)
  z_points_rot_z <- z_points_rot_y
  # Translate points to final position
  x_points_final <- x_points_rot_z + mean(x_rot_x)
  y_points_final <- y_points_rot_z + mean(y_rot_x)
  z_points_final <- z_points_rot_z + mean(z_rot_x)
  # Return cylinder coordinates and points as a list
  data.frame(list(
    X = c(x_rot_x, x_points_final),
    Y = c(y_rot_x, y_points_final),
    Z = c(z_rot_x, z_points_final)
  ))
}

# ------------------------------------------------------------------------------
# 'my.Xprod' function
#
# Computes the cross product of two vectors.
#
# @param a Numeric vector. The first vector.
# @param b Numeric vector. The second vector.
#
# @return A numeric vector that is the cross product of vectors a and b.
# ------------------------------------------------------------------------------
my.Xprod = function(a,b){
  x = c(a[2]*b[3] - a[3]*b[2] ,
        a[3]*b[1] - a[1]*b[3] ,
        a[1]*b[2] - a[2]*b[1])
  return(x)
}

# ------------------------------------------------------------------------------
# -------------------------------ransac-----------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 'my.cyl.fit' function
#
# Calculates the sum of squared distances from a cylinder surface.
#
# @param ang.rad Numeric vector. Contains the 5 cylinder parameters: rho, theta, phi, alpha and radius.
# @param P Matrix. Cylinder point cloud - xyz matrix.
#
# @return Numeric. Sum of squared distances from P to the cylinder surface with ang.rad parameters.
# ------------------------------------------------------------------------------
my.cyl.fit = function(ang.rad, P){

  # Extract cylinder parameters from the input vector
  rho = ang.rad[1]
  theta = ang.rad[2]
  phi = ang.rad[3]
  alpha = ang.rad[4]
  r = ang.rad[5]

  # Calculate the normal vector
  n = c( cos(phi) * sin(theta) , sin(phi) * sin(theta) , cos(theta) )

  # Calculate the derivative of the normal vector with respect to theta and phi
  n.theta = c( cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta) )
  n.phi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0 )

  # Calculate the normalized derivative of the normal vector with respect to phi
  n.phi.bar = n.phi / sin(theta)

  # Calculate the vector a
  a = n.theta * cos(alpha) + n.phi.bar * sin(alpha)

  # Compute the position vector of the cylinder's axis
  Q = crossprod(rho + r,n)

  # Subtract Q from each row of P to get the position vectors of the points with respect to the cylinder's axis
  P_minus_Q = sweep(P, 2, Q, "-")

  # Compute the distances from the points to the cylinder
  dists = apply(P_minus_Q, 1, function(u) sqrt(sum(my.Xprod(u, a)^2)) - r)

  # Calculate the sum of squared distances
  sumsq = sum(dists^2)

  # Return the sum of squared distances
  return(sumsq)

}

# ------------------------------------------------------------------------------
# 'my.cyl.parameters' function
#
# Estimates the cylinder parameters of a point cloud.
#
# @param P Matrix. Cylinder point cloud - xyz matrix.
# @param init Numeric vector (optional). Initial guess of cylinder parameters.
# @param opt.method Character string. Optimization method - passed on to optim.
#
# @return List. Optim output containing the optimal cylinder parameters of P, respectively: rho, theta, phi, alpha and radius.
# ------------------------------------------------------------------------------
my.cyl.parameters = function(P, init = NULL, opt.method = 'Nelder-Mead'){

  # If no initial guess is provided, calculate one
  if(is.null(init)){

    # Calculate the center of the point cloud
    center = colMeans(P)
    center[3] = 0
    rho = sqrt(sum(center^2))

    # Calculate the normal vector
    n = center / rho
    n[3] = 0
    theta = acos(n[3])
    phi = asin(n[2]/sin(theta))

    # Calculate the derivative of the normal vector with respect to theta and phi
    n.theta = c( cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta) )
    n.phi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0 )

    # Calculate the normalized derivative of the normal vector with respect to phi
    n.phi.bar = n.phi / sin(theta)

    # Define a function to calculate the vector a
    a = function(alpha){
      a1 = n.theta * cos(alpha) + n.phi.bar * sin(alpha)
      pr = crossprod(n, a1)
      return(abs(pr))
    }
    # Optimize the function a to find the optimal alpha
    alpha = optimize(f = a, interval = c(0,2*pi))[[1]]

    # Set the initial guess for the cylinder parameters
    init = c(rho, theta, phi, alpha, 0)
    init = c(rho, pi/2, 0, 0, 0)
  }

  # Use the optim function to find the optimal cylinder parameters
  out = optim(fn = my.cyl.fit, par = init, P=P, method = opt.method)

  # Return the output of the optim function
  return(out)
}

#' Fit a cylinder to a point cloud using RANSAC
#'
#' This function fits a cylinder to a point cloud using the RANSAC method. It uses parallel computing to speed up the RANSAC iterations.
#'
#' @param stem.sec A matrix of the point cloud coordinates.
#' @param select_n The number of points sampled in every RANSAC iteration.
#' @param p The proportion of inliers in stem.sec.
#' @param P The level of confidence desired.
#' @param timesN An inflation factor to multiply the number of RANSAC iterations by.
#' @param init An initial guess of cylinder parameters.
#' @param opt.method The optimization method to use. Passed on to optim.
#' @param rmse.threshold The threshold for the root mean square error (RMSE) for the optimization algorithm to stop.
#' @param cores The number of cores to use for parallel computing. Defaults to the number of cores available on the machine.
#'
#' @return A named vector with nine elements: The estimated parameters of the cylinder (rho, theta, phi, alpha, radius), the root mean square error (rmse) of the distances from the points to the cylinder, and the x, y, and z coordinates of the center of the base of the cylinder.
#'
#' @examples
#'  \dontrun{
#' # Define the cylinder attributes
#' npts = 500
#' cyl_length = 0.5
#' radius = 0.2718
#'
#' # Generate the X,Y,Z values
#' Z=runif(n = n, min = 0, max = cyl_length)
#' angs = runif(n, 0, 2*pi)
#' X = sin(angs)*radius
#' Y = cos(angs)*radius
#'
#' # Creation of a LAS object out of external data
#' cloud <- LAS(data.frame(X,Y,Z))
#'
#' # Fit a cylinder and retrun the information
#' cyl_par = ransac.cylinder(cloud, n=5, p=.9, P=.95, timesN=5, opt.method='Nelder-Mead', rmse.threshold=0.001, cores=1)
#'}
#' @export
ransac.cylinder = function(stem.sec, select_n=20, p=.9, P=.99, timesN = 5, init = NULL, opt.method = 'Nelder-Mead', rmse.threshold = 0.001, cores = detectCores(logical = F)){

  # If the number of points is less than the sample size, return NULL
  if(nrow(stem.sec) < select_n){ dt = NULL }else{

    # Calculate the number of iterations needed
    N = log(1 - P) / log(1 - p^select_n)

    # Initialize an empty matrix to store the results
    data = matrix(ncol=6, nrow=timesN*N)

    # Create a cluster with the number of cores available
    cl <- parallel::makeCluster(cores)

    # Export necessary variables to the cluster
    parallel::clusterExport(cl, c("stem.sec", "select_n", "init", "opt.method", "my.cyl.parameters", 'my.cyl.fit', 'my.Xprod'), envir = environment())

    # Use parLapply to run the loop in parallel
    data <- parallel::parLapply(cl, 1:(timesN*N), function(j) {
      # Randomly sample select_n points
      a = sample(1:nrow(stem.sec), size = select_n)
      temp = stem.sec[a,]

      # Estimate the cylinder parameters for the sampled points
      b = my.cyl.parameters(temp, init = init, opt.method = opt.method)

      # Store the parameters and the sum of squared distances
      bp = unlist(b[1:2])
      bp
    })

    # Stop the cluster
    parallel::stopCluster(cl)

    # Convert the list to a matrix
    data = do.call(rbind, data)

    # Select the parameters that give the smallest sum of squared distances
    c = which.min(data[,6])
    dt = if(length(c) > 1) data[sample(c, size = 1),] else data[c,]
  }

  # Compute the RMSE from the sum of squared distances
  rmse = sqrt(dt[6] / select_n)

  # Return the best parameters and get the xyz
  dt[6] = rmse
  dt = c(dt, px = dt[1] * sin(dt[2]) * cos(dt[3]), py = dt[1] * sin(dt[2]) * sin(dt[3]), pz = dt[1] * cos(dt[2]))
  names(dt)<-c("rho", "theta", "phi", "alpha", "radius", "rmse", "px", "py", "pz")
  return(dt)
}

# ------------------------------------------------------------------------------
# -------------------------------leastsq----------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 'my.cyl.dists' function
#
# Computes the distances from points to a cylinder.
#
# @param ang.rad Numeric vector. A vector of five elements representing the parameters of the cylinder.
#                ang.rad[1] is rho, the radial distance from the origin to the axis of the cylinder.
#                ang.rad[2] is theta, the polar angle.
#                ang.rad[3] is phi, the azimuthal angle.
#                ang.rad[4] is alpha, the angle of the cylinder's axis with respect to the z-axis.
#                ang.rad[5] is r, the radius of the cylinder.
# @param P Matrix. A matrix where each row represents a point in 3D space.
#
# @return Numeric vector. A vector of distances from each point in P to the cylinder.
# ------------------------------------------------------------------------------
my.cyl.dists = function(ang.rad, P){

  # Extract the parameters of the cylinder from the input vector
  rho = ang.rad[1]
  theta = ang.rad[2]
  phi = ang.rad[3]
  alpha = ang.rad[4]
  r = ang.rad[5]

  # Compute the direction cosines of the cylinder's axis
  n = c( cos(phi) * sin(theta) , sin(phi) * sin(theta) , cos(theta) )

  # Compute the derivatives of the direction cosines with respect to theta and phi
  n.theta = c( cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta) )
  n.phi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0 )

  # Compute the normalized derivative of the direction cosine with respect to phi
  n.phi.bar = n.phi / sin(theta)

  # Compute the direction vector of the cylinder
  a = n.theta * cos(alpha) + n.phi.bar * sin(alpha)

  # Compute the position vector of the cylinder's axis
  Q = crossprod(rho + r,n)

  # Subtract Q from each row of P to get the position vectors of the points with respect to the cylinder's axis
  P_minus_Q = sweep(P, 2, Q, "-")

  # Compute the distances from the points to the cylinder
  dists = apply(P_minus_Q, 1, function(u) sqrt(sum(my.Xprod(u, a)^2)) - r)

  return(dists)
}

# ------------------------------------------------------------------------------
# 'my.tukey.estimator' function
#
# Computes the Tukey's biweight estimator for robust regression.
#
# @param errors Numeric vector. A vector of residuals (errors) from a regression model.
# @param b Numeric. A tuning constant for the Tukey's biweight function. Default is 5.
#
# @return List. A list with two elements:
#   - Y: The scaled residuals.
#   - weights: The weights computed using Tukey's biweight function.
# ------------------------------------------------------------------------------
my.tukey.estimator = function(errors, b = 5){
  # Compute the median absolute deviation of the errors
  s = mad(errors)

  # Scale the errors by the median absolute deviation
  Y = errors / s

  # Initialize the Tukey's biweight function
  tue = Y

  # Identify the scaled errors that are larger than the tuning constant
  larger = abs(tue) > b

  # For the scaled errors that are larger than the tuning constant, set the Tukey's biweight function to 0
  tue[larger] = 0

  # For the scaled errors that are not larger than the tuning constant, compute the Tukey's biweight function
  tue[!larger] = (1-(tue[!larger]/b)^2)^2

  # Return the scaled errors and the Tukey's biweight function
  return(list(Y = Y, weights = tue))
}

# ------------------------------------------------------------------------------
# 'my.andrews.estimator' function
#
# Computes the Andrew's sine estimator for robust regression.
#
# @param errors Numeric vector. A vector of residuals (errors) from a regression model.
# @param b Numeric. A tuning constant for the Andrew's sine function. Default is 5.
#
# @return List. A list with two elements:
#   - Y: The scaled residuals.
#   - weights: The weights computed using Andrew's sine function.
# ------------------------------------------------------------------------------
my.andrews.estimator = function(errors, b = 5){
  # Compute the median absolute deviation of the errors
  s = mad(errors)

  # Scale the errors by the median absolute deviation
  Y = errors / s

  # Initialize the Andrew's sine function
  andrews = Y

  # Identify the scaled errors that are larger than the tuning constant
  larger = abs(andrews) > b

  # For the scaled errors that are larger than the tuning constant, set the Andrew's sine function to 0
  andrews[larger] = 0

  # For the scaled errors that are not larger than the tuning constant, compute the Andrew's sine function
  andrews[!larger] = sin(andrews[!larger] / b)

  # Return the scaled errors and the Andrew's sine function
  return(list(Y = Y, weights = andrews))
}

# ------------------------------------------------------------------------------
# 'my.hubers.estimator' function
#
# Computes the Huber's M-estimator for robust regression.
#
# @param errors Numeric vector. A vector of residuals (errors) from a regression model.
# @param b Numeric. A tuning constant for the Huber's M-estimation function. Default is 5.
#
# @return List. A list with two elements:
#   - Y: The scaled residuals.
#   - weights: The weights computed using Huber's M-estimation function.
# ------------------------------------------------------------------------------
my.hubers.estimator = function(errors, b = 5){
  # Compute the median absolute deviation of the errors
  s = mad(errors)

  # Scale the errors by the median absolute deviation
  Y = errors / s

  # Initialize the Huber's M-estimation function
  hubers = Y

  # Identify the scaled errors that are larger than the tuning constant
  larger = abs(hubers) > b

  # For the scaled errors that are larger than the tuning constant, set the Huber's M-estimation function to b
  hubers[larger] = b

  # For the scaled errors that are not larger than the tuning constant, compute the Huber's M-estimation function
  hubers[!larger] = abs(hubers[!larger])

  # Return the scaled errors and the Huber's M-estimation function
  return(list(Y = Y, weights = hubers))
}

# ------------------------------------------------------------------------------
# 'my.robust.estimator' function
#
# Computes a robust estimator for regression based on the provided method.
#
# @param errors Numeric vector. A vector of residuals (errors) from a regression model.
# @param b Numeric. A tuning constant for the estimator function. Default is 5.
# @param method Character string. The method to use for the estimator. One of 'tukey', 'andrews', or 'huber'. Default is 'tukey'.
#
# @return List. A list with two elements:
#   - Y: The scaled residuals.
#   - weights: The weights computed using the specified estimator function.
#
# https://online.stat.psu.edu/stat501/lesson/t/t.1/t.1.1-robust-regression-methods
# ------------------------------------------------------------------------------
my.robust.estimator = function(errors, b = 4.7, method = 'tukey'){
  if (method == 'tukey') {
    return(my.tukey.estimator(errors, b = 4.685))
  } else if (method == 'andrews') {
    return(my.andrews.estimator(errors, b = 1.339))
  } else if (method == 'huber') {
    return(my.hubers.estimator(errors, b = 1.345))
  } else {
    stop("Invalid method. Must be one of 'tukey', 'andrews', or 'huber'.")
  }
}

#' Fit a cylinder to a point cloud using least squares
#'
#' This function fits a cylinder to a point cloud using a least squares approach. It also uses a robust estimator to handle outliers in the data.
#'
#' @param P A matrix where each row represents a point in 3D space.
#' @param init An initial guess for the parameters of the cylinder. If NULL, an initial guess is computed.
#' @param select_n The number of points to randomly select from P for the estimation. Default is 20.
#' @param max.iter The maximum number of iterations for the optimization algorithm. Default is 10.
#' @param opt.method The optimization method to use. Default is 'Nelder-Mead'.
#' @param rmse.threshold The threshold for the root mean square error (RMSE) for the optimization algorithm to stop. Default is 0.001.
#' @param m.estimator The type of M-estimator to use for robust estimation. Default is 'tukey'.
#' @return A named vector with six elements: The estimated parameters of the cylinder (rho, theta, phi, alpha, radius) and the root mean square error (rmse) of the distances from the points to the cylinder.
#'
#' @examples
#'  \dontrun{
#' # Define the cylinder attributes
#' npts = 500
#' cyl_length = 0.5
#' radius = 0.2718
#'
#' # Generate the X,Y,Z values
#' Z=runif(n = n, min = 0, max = cyl_length)
#' angs = runif(n, 0, 2*pi)
#' X = sin(angs)*radius
#' Y = cos(angs)*radius
#'
#' # Creation of a LAS object out of external data
#' cloud <- LAS(data.frame(X,Y,Z))
#'
#' # Fit a cylinder and retrun the information
#' cyl_par = least.sq.cylinder(cloud, select_n=20, max.iter=10, opt.method='Nelder-Mead', rmse.threshold=0.001, m.estimator='tukey')
#'}
#' @export
least.sq.cylinder = function(P, init=NULL, select_n=20, max.iter = 10, opt.method = 'Nelder-Mead', rmse.threshold = 0.001, m.estimator = 'tukey'){

  # Initialize the weights
  w = numeric(nrow(P))

  # If no initial guess for the parameters of the cylinder is provided, compute one
  if(is.null(init)){

    # Compute the mean of the points
    center = apply(P, 2, mean)
    center[3] = 0
    rho = sqrt(sum(center^2))

    # Compute the direction cosines of the cylinder's axis
    n = center / rho
    theta = acos(n[3])
    phi = asin(n[2]/sin(theta))

    # Compute the derivatives of the direction cosines with respect to theta and phi
    n.theta = c( cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta) )
    n.phi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0 )
    n.phi.bar = n.phi / sin(theta)

    # Define a function to optimize over alpha
    a = function(alpha){
      a1 = n.theta * cos(alpha) + n.phi.bar * sin(alpha)
      pr = crossprod(n, a1)
      return(abs(pr))
    }
    # Optimize over alpha
    alpha = optimize(f = a, interval = c(0,2*pi))[[1]]

    # Set the initial guess for the parameters of the cylinder
    init = c(rho, theta, phi, alpha, 0)
  }

  # Randomly select a subset of the points
  P = P[sample(1:nrow(P), select_n, replace = F),]

  # Set the weights to 1
  w = rep(1,nrow(P))

  # Initialize the stopping criterion and the iteration counter
  crit = F
  i = 1
  ss = 0
  while(crit == F){
    # Save the previous sum of squared distances
    ss.prev = ss

    # Optimize the sum of squared distances from the points to the cylinder
    opt = optim(par = init, fn = function(u){ a=my.cyl.dists(u, P = P) ; ssq = sum(w*a^2) ; return(ssq) }, method = opt.method)

    # Update the parameters of the cylinder and the sum of squared distances
    init = opt[[1]]
    ss = opt$value

    # Compute the distances from the points to the cylinder
    dst = my.cyl.dists(ang.rad = init, P = P)

    # Compute the Tukey's biweight function for the distances
    wh = my.robust.estimator(dst, b = 4.7, method = m.estimator)

    # Update the scaled residuals and the weights
    Y = wh$Y
    w = wh$weights

    # Increment the iteration counter
    i = i+1

    # Check the stopping criterion
    crit = abs(ss - ss.prev) < rmse.threshold || i == max.iter || ss < rmse.threshold
  }

  # Compute the RMSE from the sum of squared distances
  rmse = sqrt(ss / nrow(P))

  # Return the estimated parameters of the cylinder, the RMSE, and the xyz
  pars = c(init, rmse)
  pars = c(pars, px = pars[1] * sin(pars[2]) * cos(pars[3]), py = pars[1] * sin(pars[2]) * sin(pars[3]), pz = pars[1] * cos(pars[2]))
  names(pars)<-c("rho", "theta", "phi", "alpha", "radius", "rmse", "px", "py", "pz")
  return(pars)
}

#' Fit a cylinder to a point cloud
#'
#' This function fits a cylinder to a point cloud using either the RANSAC or least squares method.
#' The function returns the estimated parameters of the cylinder, the root mean square error (RMSE) of the distances from the points to the cylinder, and the x, y, and z coordinates of the center of the base of the cylinder.
#'
#' @param pts_mat A matrix of the point cloud coordinates.
#' @param method The method to use for fitting. Either 'ransac' or 'leastsq'.
#' @param select_n The number of points sampled in every RANSAC iteration.
#' @param p The proportion of inliers in P. Only used if method is 'ransac'.
#' @param P The level of confidence desired. Only used if method is 'ransac'.
#' @param timesN An inflation factor to multiply the number of RANSAC iterations by. Only used if method is 'ransac'.
#' @param init An initial guess of cylinder parameters.
#' @param opt.method The optimization method to use. Passed on to optim.
#' @param max.iter The maximum number of iterations for the optimization algorithm. Only used if method is 'leastsq'.
#' @param rmse.threshold The threshold for the root mean square error (RMSE) for the optimization algorithm to stop.
#' @param per.allowable.error The percentage of allowable error for the RMSE.
#' @param cores The number of cores to use for parallel computing. Only used if method is 'ransac'.
#' @param m.estimator The type of M-estimator to use for robust estimation. Only used if method is 'leastsq'.
#'
#' @return A vector with nine elements: The estimated parameters of the cylinder (rho, theta, phi, alpha, radius), the root mean square error (rmse) of the distances from the points to the cylinder, and the x, y, and z coordinates of the center of the base of the cylinder.
#'
#' @examples
#'  \dontrun{
#' # Define the cylinder attributes
#' npts = 500
#' cyl_length = 0.5
#' radius = 0.2718
#'
#' # Generate the X,Y,Z values
#' Z=runif(n = n, min = 0, max = cyl_length)
#' angs = runif(n, 0, 2*pi)
#' X = sin(angs)*radius
#' Y = cos(angs)*radius
#'
#' # Creation of a LAS object out of external data
#' cloud <- LAS(data.frame(X,Y,Z))
#'
#' # Fit a cylinder and retrun the information
#' cyl_par = fit_cylinder(cloud, method = 'ransac', n=5, p=.9, P=.95, timesN=5, opt.method='Nelder-Mead', select_n=20, max.iter=10, rmse.threshold=0.001, per.allowable.error=50, cores=1, m.estimator='tukey')
#'}
#' @export
fit_cylinder = function(pts_mat, method = "ransac", select_n=20, p=.9, P=.99, timesN = 5, init = NULL, opt.method = 'Nelder-Mead',
                        max.iter = 10, rmse.threshold = 0.001, per.allowable.error = 50, cores = detectCores(logical = F), m.estimator = 'tukey'){
  if(method == "ransac"){
    result = ransac.cylinder(pts_mat, select_n, p, P, timesN, init, opt.method, rmse.threshold, cores)
  } else if(method == "leastsq"){
    result = least.sq.cylinder(pts_mat, init, select_n, max.iter, opt.method, rmse.threshold, m.estimator)
  } else {
    stop("Invalid method. Choose either 'ransac' or 'leastsq'.")
  }

  if(result[6] <= 0 || ((result[6] / result[5]) * 100) > per.allowable.error){
    error_pars = c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
    names(error_pars)<-c("rho", "theta", "phi", "alpha", "radius", "rmse", "px", "py", "pz")
    return(error_pars)
  }

  return(result)
}
