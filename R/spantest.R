#' @name spantest
#' @aliases spantest-package
#' @title spantest: Mean-Variance Spanning Tests for Portfolio Analysis
#' @description
#' \code{spantest} is an R package for conducting rigorous mean-variance spanning tests in financial and econometric analysis.
#' It implements classical and modern methods to assess whether a given set of benchmark assets spans the mean-variance
#' efficient frontier of a larger asset universe.
#' The package includes the well-known tests of Huberman and Kandel (1987), Gibbons, Ross, and Shanken (1989),
#' Kempf and Memmel (2006), Pesaran and Yamagata (2023), Kan and Zhou (2012), and Britten-Jones (1999),
#' along with recent high-dimensional and robust spanning tests by Gungor and Luger (2016) and Ardia and Sessinou (2024).
#'
#' The package supports both asset-level and global null hypothesis testing, and includes extensions for high-dimensional,
#' heteroskedastic, and serially dependent data environments using subseries-based and Cauchy Combination Test (CCT) methodologies.
#'
#' @section Functions:
#' \itemize{
#'   \item Alpha Spanning: \code{\link{span_grs}}, \code{\link{span_py}}, \code{\link{span_f1}}, \code{\link{span_gl_a}}, \code{\link{span_as}};
#'   \item Variance Spanning: \code{\link{span_km}}, \code{\link{span_f2}}, \code{\link{span_bj}}, \code{\link{span_as}};
#'   \item Joint Mean-Variance Spanning: \code{\link{span_hk}}, \code{\link{span_gl_ad}}, \code{\link{span_as}};
#' }
#'
#' @section Update:
#' The latest version of the package is available at \url{https://github.com/ArdiaD/spantest}
#'
#' @author David Ardia and Benjamin Seguin
#'
#' @note By using \code{spantest}, you agree to the following:
#' (1) You must cite Ardia and Sessinou (2025) in any working or published paper using \code{spantest},
#' (2) You must place the URL \url{https://CRAN.R-project.org/package=spantest} in a footnote or appendix,
#' (3) You assume all responsibility for use and interpretation of the results.
#'
#' @references
#' Ardia, D., Sessinou, V. (2024). Robust Inference in Large Panels and Markowitz Portfolios. Working paper.
#'
#' Gungor, S., Luger, R. (2016). Spanning tests with time-varying second moments. \emph{Journal of Financial Econometrics}, \bold{14}(3), 561--589.
#'
#' Huberman, G., Kandel, S. (1987). Mean-variance spanning. \emph{Journal of Finance}, \bold{42}(4), 873--888.
#'
#' Gibbons, M.R., Ross, S.A., Shanken, J. (1989). A test of the efficiency of a given portfolio. \emph{Econometrica}, \bold{57}(5), 1121--1152.
#'
#' Britten-Jones, M. (1999). The sampling error in estimates of mean-variance efficient portfolio weights. \emph{Journal of Finance}, \bold{54}(2), 655--671.
#'
#' Kempf, A., Memmel, C. (2006). Estimating the global minimum variance portfolio. \emph{Schmalenbach Business Review}, \bold{58}(4), 332--348.
#'
#' Bossaerts, P., Jolivet, B. (2008). Learning to coordinate: A study in portfolio efficiency. \emph{Review of Financial Studies}, \bold{21}(4), 1491--1523.
#'
#' Pesaran, M.H., Yamagata, T. (2023). Large cross-section tests of alpha spanning. \emph{Journal of Econometrics}, \bold{236}(2), 623--645.
#'
#' Kan, R., Zhou, G. (2012). Tests of mean-variance spanning. \emph{Annals of Economics and Finance}, \bold{13}(1), 145--193.
#'
#' @import stats
#' @importFrom utils head tail
NULL
