setGeneric("viterbi", function(observations) standardGeneric("viterbi"))

setMethod("viterbi", "integer", function(observations) {
    hidden_states <- vector("integer", length(observations))
    .Call(C_viterbi, PACKAGE = "tsshmm", hidden_states, observations)
    hidden_states
})

setMethod("viterbi", "IntegerList", function(observations) {
    hidden_states <- as(lapply(lengths(observations), vector,
                               mode = "integer"),
                        "IntegerList")
    .Call(C_viterbi_vectorized, PACKAGE = "tsshmm", hidden_states,
          observations)
    hidden_states
})
