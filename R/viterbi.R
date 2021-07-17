setGeneric("viterbi",
           function(observations, parallel = TRUE) standardGeneric("viterbi"),
           signature = "observations")

setMethod("viterbi", "integer", function(observations, parallel) {
    hidden_states <- vector("integer", length(observations))
    .Call(C_viterbi, PACKAGE = "tsshmm", hidden_states, observations)
    hidden_states
})

setMethod("viterbi", "IntegerList", function(observations, parallel) {
    lengths <- lengths(observations)
    hidden_states <- vector("integer", sum(lengths))
    if (parallel) {
        .Call(C_viterbi_vectorized_mt, PACKAGE = "tsshmm", hidden_states,
              unlist(observations, use.names = FALSE), lengths)
    } else {
        .Call(C_viterbi_vectorized, PACKAGE = "tsshmm", hidden_states,
              unlist(observations, use.names = FALSE), lengths)
    }
    as(relist(hidden_states, observations), "IntegerList")
})
