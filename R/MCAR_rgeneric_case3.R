utils::globalVariables("inv_eps_min")

inla.rgeneric.MCAR.model_case3 <- function(
  cmd = c(
    "graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"
  ), 
  theta = NULL
) {
  ## Implementation of MCAR(B, Sigma) model described in Jin et al. 2007 as Case
  ## n°3. They proved that the joint distribution of phi is given by
  # phi ~ N(0, (A \kr I_n) (I_p \kr D - B \kr W) ^ (-1) (A \kr I_n)')
  # where B is a p x p symmetric matrix, W is the adjacency matrix, D is a
  # diagonal matrix with elements m_i that denote the number of neighbours of
  # area i and AA' = Sigma, i.e A is the upper triangular cholesky square root
  # of Sigma. The distribution is well-defined if 1 / eps_min < xi_j < 1 with
  # indices j = 1, ..., p, where eps_min is the minimum eigenvalue of
  # D^(-0.5)WD^(-0.5) and xi_j are the eigenvalues of B. 
  #
  # This should imply that I can write the joint distribution of phi as
  # phi ~ N(0, [(A ^ (-1) \kr I_n)' (I_p \kr D - B \kr W) (A ^ (-1) \kr I_n)] ^ (-1))
  #
  # It's difficut to design a prior for B that satisfy that condition, so the
  # authors used a different approach, they represented B using it's spectral 
  # decomposition. B = PDeltaP' where P is the orthogonal matrix of eigenvectors
  # and Delta is diagonal matrix of ordered eigenvalues xi_1, ..., xi_p. 
  # The p x p orthogonal matrix P is characterized using the p * (p - 1) / 2
  # Givens angles theta_ij for i = 1, ..., p - 1 and j = i + 1, ..., p. The
  # matrix P is written as the product of p * (p - 1) / 2 matrices, each
  # associated with a Givens angle: P = G12 * G13 * ...  G_(p-1)p where the
  # indices i and j are distinct and G_ij is a p x p identity matrix with the
  # ith and jth diagonal elements replaced by cos(theta_ij), and the (i, j) and
  # (j, i) elements replaced by \pm sin(theta_ij). The Givens angles are unique
  # with a domain [-pi/2, pi/2] and the eigenvalues xi_j are in range
  # (1/eps_min, 1), the authors put a Uniform prior (-pi/2, pi/2) on theta_ij
  # and U(1/eps_min, 1) prior on xi_j.
  # Sigma is a positive definite covariance matrix and the inverse Wishart prior
  # distribution renders itself as a natural choice.
  
  # Start with p = 2
  # The INLA-vector theta is composed by 6 elements:
  # 1, 2: xi_1 and xi_2, the eigenvalues of B. They must lie in the range
  #       (1/eps_min , 1) so we used a "sort-of" inverse logit 
  #       transformation. We assign to them a Uniform prior between 
  #       1 / eps_min and 1. The matrix Delta is a diagonal matrix with 
  #       elements xi_1 and xi_2. 
  # 3   : theta_12, which is the first (and only) Givens angle. It must lie in
  #       the range [-pi/2, pi/2] (or [0, pi/2]) so we assigned a uniform 
  #       prior between -pi/2 and pi/2. The matrix P is created as 
  #       [cos(theta_12) & sin(theta_12) \\ -sin(theta_12) & cos(theta_12)]. 
  #       The matrix B is created as P Delta P'
  # 4, 5: the marginal precisions of Sigma
  # 6   : the correlation parameter. 
  # The elements 4, 5, 6 are combined for creating the matrix Sigma
  interpret.theta <- function() {
    # Function for changing from internal scale to external scale
    
    # theta[1] and theta[2] are the internal representations of xi_1 and xi_2.
    # theta[1] and theta[2] \in R so we must use a "sort of" logist
    # transformation, bounded between 1 / eps_min and 1.
    xi = vapply(
      theta[1:2], 
      function(x) inv_eps_min + (1 - inv_eps_min) * stats::plogis(x), 
      numeric(1)
    )
    
    # Build matrix Delta
    Delta <- Matrix::Diagonal(n = 2, x = xi)
    
    # theta[3] is the internal representation of theta_12. theta[3] \in R while
    # theta_12 \in [-pi/2, pi/2] so we must use a logit transformation. 
    theta_12 = -pi / 2 + pi * stats::plogis(theta[3])
    
    # Build matrix P
    P <- Matrix::Matrix(
      c(cos(theta_12), sin(theta_12), -sin(theta_12), cos(theta_12)), 
      nrow = 2, 
      ncol = 2, 
      byrow = TRUE
      )
    
    # Build matrix B
    B <- P %*% Delta %*% t(P)
    
    # theta[4] and theta[5] are the log-precisions, which are defined in R. I
    # take the exp to obtain the precisions and then I will derive the standard
    # deviations.
    mprec = vapply(
      theta[4:5], 
      exp, 
      numeric(1)
    )
    # theta[6] is the logit transformation of the correlation parameter. I need
    # to transform it back and I will use to determine Sigma.
    corre = 2 * stats::plogis(theta[6]) - 1 
    
    # Now I want to build Sigma
    
    # intial matrix with 1s at the diagonal
    M <- diag(1, 2)
    
    # Adding correlation parameters (lower.tri) and (upper.tri)
    M[lower.tri(M)] <- corre
    M[upper.tri(M)] <- t(M)[upper.tri(M)]
    
    # Preparing the st. dev matrix
    st.dev <- sqrt(1 / mprec)
    
    # Matrix of st. dev.
    st.dev.mat <- matrix(st.dev, ncol = 1) %*% matrix(st.dev, nrow = 1)
    
    # Build Sigma
    Sigma <- M * st.dev.mat
    PREC <- solve(Sigma)
    
    return(list(Delta = Delta, P = P, B = B, Sigma = Sigma, PREC = PREC))    
  }
  
  # Graph of precision function; i.e., a 0/1 representation of precision matrix
  graph <- function() {
    # I use simple symmetric diagonal matrix for Sigma and B
    Sigma <- B <- Diagonal(2, 0.5)
    A <- chol(Sigma)
    D <- Diagonal(nrow(W), rowSums(W))
    
    # I build the left/right part of the product, i.e. A ^ (-1) \kr I_n
    A_kron_I <- kronecker(
      solve(A), 
      Diagonal(nrow(W), 1)
    )
    # I build the central block of the product, i.e. I_p \kr D - B \kr W
    central_block <- kronecker(Diagonal(2, 1), D) - kronecker(B, W)
    G <- t(A_kron_I) %*% central_block %*% A_kron_I
    G
  }
  
  Q <- function() {
    # Parameters in model scale
    param <- interpret.theta()
    
    # I extract Sigma and calculate A
    Sigma <- param$Sigma
    A <- t(chol(Sigma))
    
    # Calculate the left part of the product
    A_kron_I <- kronecker(
      solve(A), 
      Diagonal(nrow(W), 1)
    )
    
    # I calculate the central block of the product
    B <- param$B
    D <- Diagonal(nrow(W), rowSums(W))
    central_block <- kronecker(Diagonal(2, 1), D) - kronecker(B, W)
    
    Q <- t(A_kron_I) %*% central_block %*% A_kron_I
    Q
  }
  
  # Mean of model
  mu <- function() {
    return(numeric(0))
  }
  
  # return the log(normalising constant) for the model
  log.norm.const <- function() {
    return(numeric(0))
  }
  
  log.prior <- function() {
    ## return the log-prior for the hyperparameters
    param <- interpret.theta()
    
    # I set a Uniform prior for xi_1 and xi_2 between -pi/2 and pi/2, so I need
    # to add the log of density of theta[1] and theta[2] which is calculated
    # using the transformation of r.v. theorem.
    # Same stuff for theta_12
    val <- 0
    for (j in 1:3) {
      val <- val - theta[j] - 2 * log(1 + exp(-theta[j]))
    }
    
    # Whishart prior for joint matrix of hyperparameters
    val <- val + log(MCMCpack::dwish(W = param$PREC, v = 2, S = diag(rep(1, 2))))
    # This is for precisions
    val <- val + sum(theta[4:5])
    # This is for correlation terms
    val <- val + sum(
      log(2) + theta[6] - 2 * log(1 + exp(theta[6]))
    )
    
    return(val)
  }
  
  initial <- function() {
    ## return initial values
    # The Initial values form a diagonal matrix
    return(c(rep(0, 2), (1 + inv_eps_min) / 2, rep(0, 3)))
  }
  
  quit <- function() {
    return(invisible())
  }
  
  val <- do.call(match.arg(cmd), args = list())
  return(val)
}

#' B
#'
#' @param ... ABC
#' @return B
#' @export
#'
#' @examples
#' 1 + 1
MCAR.model_case3 <- function(...) {
  INLA::inla.rgeneric.define(model = inla.rgeneric.MCAR.model_case3, ...)
}

