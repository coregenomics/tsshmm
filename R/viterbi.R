viterbi_low_level <- function(model, observations) {
    if (is.list(observations) || is(observations, "List")) {
        lengths <- lengths(observations)
    } else {
        lengths <- length(observations)
    }
    hidden_states <- vector("integer", sum(lengths))
    .Call(C_viterbi, PACKAGE = "tsshmm", hidden_states, model@external_pointer,
          unlist(observations, use.names = FALSE), lengths)
    if (is.list(observations) || is(observations, "List")) {
        as(relist(hidden_states, observations), "IntegerList")
    } else {
        hidden_states
    }
}

## The theory and calculation of this minimum distance between promoters is
## explained in the vignette appendix.
prom_dist <- function(model, tol = 10^-8) {
    params <- parameters(model)
    ## Replace NA values with zeroes for the eigen() function.
    A <- params$trans
    A[is.na(A)] <- 0
    eigen <- eigen(A)
    ## Split A into diagonal matrix L and factorized H and H-inverse.
    l <- eigen$values
    L <- l * diag(length(l))
    H <- eigen$vectors
    H_inv <- solve(H)
    ## Validate factorized matrices multiply into the identity matrix.
    stopifnot(all(zapsmall(H %*% H_inv) == diag(length(l))))
    ## Validate multiplying factorized with diagonal matrices recreate A.
    stopifnot(all(zapsmall(H %*% L %*% H_inv) == A))
    ## Create our HMM emissions "psi" function.
    B <- t(params$emis)
    psi <- function(pi, n = 1) c(B %*% t(H %*% L^n %*% H_inv) %*% unlist(pi))
    ## Calculate psi for all combinations of initial values up to 100 windows.
    pis <- unlist(apply(diag(length(l)), 1, list), recursive = FALSE)
    grid <- expand.grid(n = 1:100, pi = pis)
    dd <- mapply(psi, grid$pi, grid$n)
    ## Set emission names.
    rownames(dd) <- c("bg", "enriched", "depleted")
    ## Create emission probability matrix.
    probs <- cbind(grid, t(dd))
    ## All emission probabilities should sum to 1.
    stopifnot(all(zapsmall(rowSums(probs[, -2:-1])) == 1))
    ## Calculate relative differences between probabilities.
    diffs_df <- probs
    diffs_df[, "pi"] <-
        sapply(diffs_df$pi,
               function(x) which(unlist(x) == 1)) # Parsable for grouping.
    relative_diff <- function(x) diff(c(NA, x)) / x
    for (col in colnames(probs)[-2:-1]) {
        formula <- formula(eval(parse(text = paste(col, "~ pi"))))
        y <- aggregate(formula, diffs_df, relative_diff)
        diffs_df[, col] <- c(y[[2]])
    }
    ## Locate worst case scenario.
    diffs_df[, "max"] <- apply(diffs_df[, -2:-1], 1, max)
    diffs <- aggregate(max ~ n, diffs_df, max)$max
    ## Lowest index (number of windows) from tolerance.
    min(which(diffs < tol))
}
