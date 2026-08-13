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

# THE definition of the twelve presets. span_simulate() reads this; nothing
# restates it. It used to be written out twice -- once here and once as a list
# inside span_simulate() -- and the two drifted: the table said one thing about
# the Student presets and the simulator did another.
#
# The grid is a 3 x 4 factorial, one innovation law per family repeated across
# the four dynamics, so that moving down a column changes the DYNAMICS and
# nothing else. Keep it that way: an innovation law that varies with the dynamics
# makes the rows of the manuscript's size tables incomparable.
.DGP_PRESETS <- data.frame(
  preset      = 1:12,
  name        = DGP_NAMES,
  innovation  = rep(c("normal", "t", "skew-t"), times = 4),
  dynamics    = rep(c("iid", "garch", "ar-garch", "ar"), each = 3),
  df          = rep(c(5, 5, 4), times = 4),        # unused for "normal"
  xi          = rep(c(NA, NA, 0.9), times = 4),    # skew-t only
  standardize = TRUE,
  # The twelve are meant to differ in the SHAPE of the dependence and the
  # thickness of the tails, not in how much variance they carry. The GARCH
  # parameters were already picked so that omega/(1 - alpha - beta) = 1; scale =
  # "process" does the same for the AR families, which would otherwise carry
  # 1/(1 - phi^2) = 1.042. It changes no size result -- at ncp = 0 it rescales the
  # whole system and the tests are scale-invariant -- and moves the
  # signal-to-noise ratio under the alternative by 2.1%.
  scale       = "process",
  stringsAsFactors = FALSE
)


#' The simulation processes, by name and preset number
#'
#' Returns the correspondence between the names accepted by [span_simulate()]'s
#' `dgp` argument, the preset numbers, and the process each one denotes. Useful
#' when reading code or output that refers to a process by number.
#'
#' @return A data frame with one row per process and columns `preset` (integer),
#'   `name` (character), `innovation`, `dynamics`, `df`, `standardize` and
#'   `scale`. Together they are enough to rebuild any preset by hand with
#'   [span_simulate()]'s explicit arguments; a table that described a process
#'   only partly is how a discrepancy between the two once went unnoticed.
#'
#'   The twelve are a 3 x 4 factorial of innovation law by dynamics, so a column
#'   varies the dynamics alone: every heavy-tailed preset is standardised, and
#'   `scale = "process"` gives each one unit unconditional variance whatever its
#'   dynamics.
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
# Every column a caller would need to rebuild the preset by hand. A table that
# describes a process only partly is how the last discrepancy survived: the
# reconstruction silently picked up a default the preset does not use.
span_dgp_table <- function() .DGP_PRESETS[, c("preset", "name", "innovation",
                                              "dynamics", "df", "standardize",
                                              "scale")]
