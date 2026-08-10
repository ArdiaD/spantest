# =============================================================================
# Names for the twelve simulation processes.
#
# The processes are a 3 x 4 factorial: innovation law (normal / Student-t /
# skew-t) crossed with dynamics (i.i.d. / GARCH / AR / AR-GARCH). Referring to
# them by position in that list forces every reader to look the number up, and
# the lookup is easy to get wrong -- downstream code and write-ups have ordered
# the AR and AR-GARCH blocks differently from this file. A name carries the cell
# it denotes, so it cannot be transposed the way a number can.
#
# The vocabulary is the one used in the write-ups: N / ST / SKST for the
# innovations, iid / GARCH / AR / AR-GARCH for the dynamics.
# =============================================================================

#' Names of the built-in simulation processes
#'
#' A character vector of length twelve giving a name for each preset of
#' [span_simulate()]. **Position in this vector is the preset number**, so
#' `DGP_NAMES[9]` is the name of `dgp = 9`.
#'
#' @format A character vector of length 12.
#' @seealso [span_dgp_table()], [span_simulate()]
#' @export
DGP_NAMES <- c(
  "iid-N",         "iid-ST",         "iid-SKST",        # presets 1-3
  "GARCH-N",       "GARCH-ST",       "GARCH-SKST",      # presets 4-6
  "AR-GARCH-N",    "AR-GARCH-ST",    "AR-GARCH-SKST",   # presets 7-9
  "AR-N",          "AR-ST",          "AR-SKST"          # presets 10-12
)

#' The simulation processes, by name and preset number
#'
#' Returns the correspondence between the names accepted by [span_simulate()]'s
#' `dgp` argument, the preset numbers, and the process each one denotes. Useful
#' when reading code or output that refers to a process by number.
#'
#' @return A data frame with one row per process and columns `preset` (integer),
#'   `name` (character), `innovation`, `dynamics`, `df` and `standardize`. The
#'   last two record the preset's own settings: the Student-t presets `iid-ST`
#'   and `AR-ST` use *raw*, non-standardised \eqn{t_5} innovations, so they are
#'   not reproduced by passing `innovation = "t"` with the default arguments.
#'
#' @examples
#' span_dgp_table()
#'
#' # the same process, two ways
#' set.seed(1); a <- span_simulate(100, 2, 3, dgp = "AR-GARCH-SKST")
#' set.seed(1); b <- span_simulate(100, 2, 3, dgp = 9)
#' identical(a, b)
#'
#' @seealso [span_simulate()], [DGP_NAMES]
#' @export
span_dgp_table <- function() {
  data.frame(
    preset      = 1:12,
    name        = DGP_NAMES,
    innovation  = c("normal", "t", "skew-t", "normal", "t", "skew-t",
                    "normal", "t", "skew-t", "normal", "t", "skew-t"),
    dynamics    = c("iid", "iid", "iid", "garch", "garch", "garch",
                    "ar-garch", "ar-garch", "ar-garch", "ar", "ar", "ar"),
    df          = c(5, 5, 4, 5, 4, 4, 5, 4, 4, 5, 5, 4),
    standardize = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
                    TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}
