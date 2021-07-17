setGeneric("viterbi", function(observations) standardGeneric("viterbi"))

setMethod("viterbi", "integer", function(observations) {
    hidden_states <- vector("integer", length(observations))
    .Call(C_viterbi, PACKAGE = "tsshmm", hidden_states, observations)
    hidden_states
})

setMethod("viterbi", "IntegerList", function(observations) {
    lengths <- lengths(observations)
    hidden_states <- vector("integer", sum(lengths))
    .Call(C_viterbi_vectorized, PACKAGE = "tsshmm", hidden_states,
          unlist(observations, use.names = FALSE), lengths)
    as(relist(hidden_states, observations), "IntegerList")
})
