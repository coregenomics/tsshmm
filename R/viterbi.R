viterbi <- function(observations) {
    hidden_states <- vector("integer", length(observations))
    .Call(C_viterbi, PACKAGE = "tsshmm", hidden_states, observations)
    hidden_states
}
