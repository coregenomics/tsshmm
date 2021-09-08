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
