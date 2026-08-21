21st Century Supermoon Estimation in R
================

``` r
knitr::opts_chunk$set(cache = FALSE, results = "hold") # Group output by chunk.

# Suppress inefficiency warnings when converting mpfr values to numerics.
# N.B. "mpfr" is Multiple Precision Floating-Point Reliable Arithmetic, a C
#      library. Rmpfr is an R package that provides an R interface to mpfr.
library(pkgcond)

# Arbitrary-precision arithmetic.
suppressMessages(library(Rmpfr))

# The following packages are also used explicitly...
#
#   extrafont  E.g., for extra fonts in graphs.
#   ggplot2    E.g., for plotting graphs.
#   lubridate  E.g., for date-time manipulation utilities.
#   scales     E.g., for formatting date axes in graphs.

# Load the shared object holding Fortran code for ELP-2000/82.
# Designate the shared library extension in a platform-independent manner.
if (isFALSE(is.loaded("elp82b2", PACKAGE = "elp82b_2", type = "Fortran"))) {
  dyn.load(file.path(paste0("elp82b_2", .Platform$dynlib.ext)))
  stopifnot(is.loaded("elp82b2", PACKAGE = "elp82b_2", type = "Fortran"))
}


# Start global variables. #####################################################

# Set output verbosity.
g_verbosity <- FALSE

# Set lunar apsis.
# Alternatively, g_lunar_apsis <- "apogee"
g_lunar_apsis <- "perigee"

# Set lunar phase.
# Alternatively, g_lunar_phase <- "new moon"
# Alternatively, g_lunar_phase <- "first quarter"
g_lunar_phase <- "full moon"
# Alternatively, g_lunar_phase <- "last quarter"

# Set general purpose numeric global variables.
g_max_indent <- 27
g_2_mpfr <- mpfr("2", 128L)

# End global variables. #######################################################

calculate_with_uncertainty <- function(function_name, value, uncertainty,
                                       receive_uncertainty) {

  # Call another function, propagating uncertainty.
  #   N.B. It assumes the function called also returns a value with an
  #        associated uncertainty.
  #   N.B. If time-dates are being manipulated, pass the uncertainty as a
  #        lubridate duration.
  #
  # Input:
  #   function_name          The function to invoke.
  #   value                  The value to pass to the function.
  #   uncertainty            The value's uncertainty to pass to the function.
  #   receive_uncertainty    A Boolean indicating if the invoked function
  #                          itself returns a value with an uncertainty.

  # Allow for both zero and non-zero uncertainty.

  if ((value + uncertainty) == value) {
    # Zero uncertainty.
    result_mid <- function_name(value)

    if (receive_uncertainty) {
      result_range <- range(
        result_mid$value - result_mid$uncertainty,
        result_mid$value + result_mid$uncertainty
      )
    } else {
      result_range <- range(result_mid)
    }

  } else {
    # Non-zero uncertainty.
    result_low <- function_name(value - uncertainty)
    result_mid <- function_name(value)
    result_high <- function_name(value + uncertainty)

    if (receive_uncertainty) {
      result_range <- range(
        result_low$value - result_low$uncertainty,
        result_low$value + result_low$uncertainty,
        result_mid$value - result_mid$uncertainty,
        result_mid$value + result_mid$uncertainty,
        result_high$value - result_high$uncertainty,
        result_high$value + result_high$uncertainty
      )
    } else {
      result_range <- range(result_low, result_mid, result_high)
    }
  }

  # result_range is a 2-element vector holding the min and max of input args.
  # result_range[1] = min value.
  # result_range[2] = max value.
  list(
    "value" = mean(result_range),
    "uncertainty" = (result_range[2] - result_range[1]) * 0.5
  )
}

call_elp82b_2 <- function(tjj_val, prec_val = 0.0, nul_val = 10) {

  # This function calls a third-party Fortran lunar solution ELP 2000-82B.
  # The geocentric position of the moon at a given Julian date is returned.
  # That position is given in rectangular coordinates (x, y, z), in
  # kilometres, with a JD 2451545.0 reference frame.
  #
  # Input:
  #   tjj_val   A Julian date. See the "tjj" Fortran parameter, below.
  #
  #   prec_val  A truncation level. See the "prec" Fortran parameter, below.
  #               Optional.
  #
  #   nul_val   A file identifier. See the "nul" Fortran parameter, below.
  #               Optional.
  #
  # Output:
  #  <unnamed>  A list holding "r" and "ierr" output from the Fortran function.
  #             See the elp82b2 output Fortran parameters, below.
  #
  # elp82b2 (Fortran function) parameters:
  #
  # Input:
  #   tjj   A Julian date, in barycentric dynamical time.
  #           Type: real double precision.
  #   prec  Truncation level, in radians.
  #           Type: real double precision.
  #   nul   Logical identifier to connect the program to data files.
  #           N.B. Try and keep it >= 10, lower values are often reserved.
  #           Type: integer.
  #
  # Output:
  #   r     Rectangular coordinates, with a JD 2451545.0 reference frame.
  #           Type: 3-element array, real double precision.
  #   ierr  Error code.
  #           Type: integer.

  # N.B. Pass variables, not raw values, as Fortran function parameters.
  # N.B. Unwrap mpfr values before passing them to Fortran.
  tjj <- as.double(as.numeric(tjj_val))
  prec <- as.double(prec_val)
  nul <- as.integer(nul_val)
  r <- double(3)
  ierr <- integer(1)

  # Call the elp82b2 Fortran function.
  # N.B. Using parameter names in ".Fortran" crashes the session.
  #      Potential name mangling issues.
  # The geocentric lunar position vector is returned in element four.
  .Fortran("elp82b2", tjj, prec, nul, r, ierr)[[4]]
}

get_list_estimates <- function(function_name, date_start, date_end, interval) {

  # Get a list of dates (in JDE format) for nominated apsides or phases over a
  # selected date range.
  #
  #   date_start    Start of the date range, in Coordinated Universal Time.
  #   date_end      End of the date range, in Coordinated Universal Time.
  #   interval      Decimal specifying how often, in days, to check for the
  #                 nearest apside or phase, e.g. 29.4.
  #                 N.B. For apsides, ensure it's < the anomalistic month over
  #                      the date range, else apsides will be missed.
  #                 N.B. For phases, ensure it's < the synodic month over the
  #                      date range, else phases will be missed.
  #                 N.B. An mpfr value is recommended.

  # Total list elements needed is:
  #   ceiling((date_end - date_start) / lubridate::ddays(interval)).
  interval_days <- lubridate::ddays(interval)
  list_estimates <- lapply(
    0:(ceiling((date_end - date_start) / interval_days) - 1),
    function(x) {
      target_dynamical <- time_universal_to_dynamical(
        date_start + x * interval_days
      )
      calculate_with_uncertainty(
        function_name, target_dynamical$value,
        target_dynamical$uncertainty, TRUE
      )
    }
  )

  # Remove duplicates.
  if (g_verbosity) {
    cat("\nList length, initial: ", length(list_estimates), sep = "")
    list_estimates <- unique(list_estimates)
    cat("\nList length, final: ", length(list_estimates), sep = "")
  } else {
    list_estimates <- unique(list_estimates)
  }

  list_estimates
}

is_valid_date <- function(date_target) {

  # This program is designed to process 21st century dates.
  #   I.e. 01-Jan-2001 to 31-Dec-2100.
  #
  # Input:
  #   date_target    Date should be a UTC time.
  #
  # Output:
  #   A Boolean indicating if the date is a 21st century date.

  (date_target >= as.POSIXct(
    "01-01-2001 00:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"
  )) & (date_target < as.POSIXct(
    "01-01-2101 00:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"
  ))
}

make_angle_degrees_to_radians <- function() {

  # Closure constants.
  c_360_mpfr <- mpfr("360", 128L)
  c_pi_on_180_mpfr <- Const("pi", prec = 128L) / mpfr("180", 128L)

  # Body of: angle_degrees_to_radians.
  function(angle_degrees) {

    # Convert angles in degrees to radians.
    #
    # Input:
    #   angle_degrees  One or more angles in degrees.
    #
    # Output:
    #   One or more angles in radians.

    # Ensure angles are (0 >= theta > 360) degrees.
    angle_degrees <- angle_degrees %% c_360_mpfr

    # Convert to radians.
    angle_degrees * c_pi_on_180_mpfr
  }
}

make_angle_radians_to_degrees <- function() {

  # Closure constants.
  c_two_pi_mpfr <- Const("pi", prec = 128L) * mpfr("2", 128L)
  c_180_on_pi_mpfr <- mpfr("180", 128L) / Const("pi", prec = 128L)

  # Body of: angle_radians_to_degrees.
  function(angle_radians) {

    # Convert angles in radians to degrees.
    #
    # Input:
    #   angle_radians  One or more angles in radians.
    #
    # Output:
    #   One or more angles in degrees.

    # Ensure angles are (0 >= theta > 2*Pi) radians.
    angle_radians <- angle_radians %% c_two_pi_mpfr

    # Convert to degrees.
    angle_radians * c_180_on_pi_mpfr
  }
}

make_get_del_t <- function() {

  # Closure constants.
  #
  c_0_36525_mpfr <- mpfr("0.36525", 128L)
  c_0_0004_mpfr <- mpfr("0.0004", 128L)
  c_1_mpfr <- mpfr("1", 128L)
  c_2005_mpfr <- mpfr("2005", 128L)
  c_2050_mpfr <- mpfr("2050", 128L)
  #
  c_threshold_get_del_t <- mpfr(
    lubridate::decimal_date(
      as.POSIXct(
        "01-01-2026 00:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"
      ) - lubridate::dmilliseconds(mpfr("900", 128L))
    ), 128L
  )
  #
  c_uncertainty_get_del_t <- mpfr("0.058", 128L) / 3

  # Body of: get_del_t.
  function(date_target) {

    # Calculate the difference, del T, in seconds, between dynamical time and
    # universal time (UT1).
    #   See Five Millennium Canon of Solar Eclipses: -1999 to +3000
    #                                                (2000 BCE to 3000 CE).
    #
    # Input:
    #   date_target    Date should be universal time (UT1).
    #                  It's assumed the target date will be in the 21st century.
    #
    # N.B. An exact decimal year is used here. The reference material uses a
    #      less accurate value rounded to the middle of months.
    # N.B. Where relevant, a lunar secular acceleration of -26 arcsec/cy^2 is
    #      assumed.

    # Convert input date to decimal years.
    #   E.g. 1987.25 is the end of March 1987.
    #
    date_target <- mpfr(lubridate::decimal_date(date_target), 128L)

    # Before 2026...
    if (date_target < c_threshold_get_del_t) {
      # As known del T values exist until 31-Dec-2025 23:59:59, del T shouldn't
      # be estimated on or prior to that.
      # Allow a 900 millisecond buffer for processing the UT1-UTC uncertainty,
      # which is always <= 0.9 s.
      # So, c_threshold_get_del_t is (01-01-2026 00:00:00) - 0.9 s.

      cat("\nError in function get_del_t...")
      cat("\nInvalid date processing; date_target:\n")
      print(date_target, digits = 20)
      stopifnot(FALSE)
    }

    # Before 2050...
    if (date_target < c_2050_mpfr) {
      del_t <- get_del_t_pre_2050(date_target)

    } else { # After 2050...
      del_t <- get_del_t_post_2050(date_target)
    }

    # N.B. Here uncertainty is in seconds; algorithm's intrinsic uncertainty.
    n <- date_target - c_2005_mpfr
    list(
      "value" = del_t,
      "uncertainty" = c_0_36525_mpfr * n *
        sqrt(n * c_uncertainty_get_del_t * (c_1_mpfr + n * c_0_0004_mpfr))
    )
  }
}

make_get_del_t_post_2050 <- function() {

  # Closure constants.
  c_0_0032_mpfr <- mpfr("0.0032", 128L)
  c_11_0852_mpfr <- mpfr("11.0852", 128L)
  c_9369_66_mpfr <- mpfr("9369.66", 128L)

  # Body of: get_del_t_post_2050.
  function(t_mpfr) {

    # Calculate "delta t" for years including and after 2050.
    #
    # Delta t is the difference in the Earth's rotation period and the length
    # of an ideal day measured by International Atomic Time (86,400 seconds).
    #
    # Input
    #   t_mpfr  The value of a target decimal year.
    #           Assumed to be a mpfr value.
    #

    # We have, del t = -20.0 - 0.5628 * (2150.0 - target date) +
    #                  32.0 * (((target date - 1820.0) * 0.01) ^ 2)
    # Or more simply...

    t_mpfr * (c_0_0032_mpfr * t_mpfr - c_11_0852_mpfr) + c_9369_66_mpfr
  }
}

make_get_del_t_pre_2050 <- function() {

  # Closure constants.
  c_0_005589_mpfr <- mpfr("0.005589", 128L)
  c_0_32217_mpfr <- mpfr("0.32217", 128L)
  c_62_92_mpfr <- mpfr("62.92", 128L)
  c_2000_mpfr <- mpfr("2000", 128L)

  # Body of: get_del_t_pre_2050.
  function(t_mpfr) {

    # Calculate "delta t" for years before 2050.
    #
    # Delta t is the difference in the Earth's rotation period and the length
    # of an ideal day measured by International Atomic Time (86,400 seconds).
    #
    # Input
    #   t_mpfr  The value of a target decimal year.
    #           Assumed to be a mpfr value.

    t_mpfr <- t_mpfr - c_2000_mpfr
    t_mpfr * (c_0_005589_mpfr * t_mpfr + c_0_32217_mpfr) + c_62_92_mpfr
  }
}

make_get_eccentric_anomaly <- function() {

  # Closure constants.
  c_0_0000074_mpfr <- mpfr("0.0000074", 128L)
  c_0_002516_mpfr <- mpfr("0.002516", 128L)
  c_1_mpfr <- mpfr("1", 128L)

  # Body of: get_eccentric_anomaly.
  function(time_target) {

    # Get the eccentric anomaly, for a given time
    #   (in Julian centuries from Epoch J2000.0).
    #
    # Input:
    #   time_target    Target time, in Julian centuries from Epoch J2000.0.
    #                  It's assumed it will be a mpfr value.

    c_1_mpfr - time_target * (c_0_002516_mpfr + time_target * c_0_0000074_mpfr)
  }
}

make_get_lunar_apsis_jde <- function() {

  # Closure constants.
  #
  # Lunar apogee calculation coefficients.
  c_vec_apogee_coeff <- mpfr(c(
    "0.4392",         "0.0684",         "0.0456",         "0.0426",
    "0.0212",         "-0.0189",        "0.0144",         "0.0113",
    "0.0047",         "0.0036",         "0.0035",         "0.0034",
    "-0.0034",        "0.0022",         "-0.0017",        "0.0013",
    "0.0011",         "0.0010",         "0.0009",         "0.0007",
    "0.0006",         rep("0.0005", 2), rep("0.0004", 3), rep("-0.0004", 2),
    rep("0.0003", 3), "-0.0003"
  ), 128L)
  #
  c_vec_apogee_subtrahend <- mpfr(c(
    rep("0", 2), rep("-0.00011", 2), rep("0", 28)
  ), 128L)
  #
  # Lunar perigee calculation coefficients.
  c_vec_perigee_coeff <- mpfr(c(
    "-1.6769",         "0.4589",          "-0.1856",
    "0.0883",          "-0.0773",         "0.0502",
    "-0.0460",         "0.0422",          "-0.0256",
    "0.0253",          "0.0237",          "0.0162",
    "-0.0145",         "0.0129",          "-0.0112",
    "-0.0104",         "0.0086",          "0.0069",
    "0.0066",          "-0.0053",         "-0.0052",
    "-0.0046",         "-0.0041",         "0.0040",
    "0.0032",          "-0.0032",         "0.0031",
    "-0.0029",         rep("0.0027", 2),  "-0.0027",
    "0.0024",          rep("-0.0021", 3), "0.0019",
    "-0.0018",         rep("-0.0014", 3), "0.0014",
    "-0.0014",         rep("0.0013", 2),  "0.0011",
    "-0.0011",         "-0.0010",         "-0.0009",
    "-0.0008",         rep("0.0008", 2),  rep("0.0007", 3),
    rep("-0.0006", 2), "0.0006",          rep("0.0005", 2),
    "-0.0004"
  ), 128L)
  c_vec_perigee_subtrahend <- mpfr(c(
    rep("0", 4), "0.00019", "-0.00013", "0", "-0.00011", rep("0", 52)
  ), 128L)
  #
  c_neg_0_0006691_mpfr <- mpfr("-0.0006691", 128L)
  c_0_0000000052_mpfr  <- mpfr("0.0000000052", 128L)
  c_0_000000055_mpfr   <- mpfr("0.000000055", 128L)
  c_0_0000010_mpfr     <- mpfr("0.0000010", 128L)
  c_0_000001098_mpfr   <- mpfr("0.000001098", 128L)
  c_0_00001156_mpfr    <- mpfr("0.00001156", 128L)
  c_0_0000148_mpfr     <- mpfr("0.0000148", 128L)
  c_0_0008130_mpfr     <- mpfr("0.0008130", 128L)
  c_0_0100383_mpfr     <- mpfr("0.0100383", 128L)
  c_0_0125053_mpfr     <- mpfr("0.0125053", 128L)
  c_0_5_mpfr           <- mpfr("0.5", 128L)
  c_13_2555_mpfr       <- mpfr("13.2555", 128L)
  c_27_1577721_mpfr    <- mpfr("27.1577721", 128L)
  c_27_55454989_mpfr   <- mpfr("27.55454989", 128L)
  c_171_9179_mpfr      <- mpfr("171.9179", 128L)
  c_316_6109_mpfr      <- mpfr("316.6109", 128L)
  c_335_9106046_mpfr   <- mpfr("335.9106046", 128L)
  c_347_3477_mpfr      <- mpfr("347.3477", 128L)
  c_364_5287911_mpfr   <- mpfr("364.5287911", 128L)
  c_1325_55_mpfr       <- mpfr("1325.55", 128L)
  c_1999_97_mpfr       <- mpfr("1999.97", 128L)
  c_2451534_6698_mpfr  <- mpfr("2451534.6698", 128L)
  #
  c_apsis_uncertainty_apogee <- mpfr("3", 128L) / mpfr("1440", 128L)
  c_apsis_uncertainty_perigee <- mpfr("31", 128L) / mpfr("1440", 128L)

  # Perigee multiplication coefficients.
  #   Column 1: d. (d = the Moon's mean elongation, in radians.)
  #   Column 2: m. (m = the Sun's mean anomaly, in radians.)
  #   Column 3: f. (f = Moon's argument of latitude, in radians.)
  perigee_matrix <- matrix(c(
    2,   0,  0,  #  1
    4,   0,  0,  #  2
    6,   0,  0,  #  3
    8,   0,  0,  #  4
    2,  -1,  0,  #  5
    0,   1,  0,  #  6
    10,  0,  0,  #  7
    4,  -1,  0,  #  8
    6,  -1,  0,  #  9
    12,  0,  0,  # 10
    1,   0,  0,  # 11
    8,  -1,  0,  # 12
    14,  0,  0,  # 13
    0,   0,  2,  # 14
    3,   0,  0,  # 15
    10, -1,  0,  # 16
    16,  0,  0,  # 17
    12, -1,  0,  # 18
    5,   0,  0,  # 19
    2,   0,  2,  # 20
    18,  0,  0,  # 21
    14, -1,  0,  # 22
    7,   0,  0,  # 23
    2,   1,  0,  # 24
    20,  0,  0,  # 25
    1,   1,  0,  # 26
    16, -1,  0,  # 27
    4,   1,  0,  # 28
    9,   0,  0,  # 29
    4,   0,  2,  # 30
    2,  -2,  0,  # 31
    4,  -2,  0,  # 32
    6,  -2,  0,  # 33
    22,  0,  0,  # 34
    18, -1,  0,  # 35
    6,   1,  0,  # 36
    11,  0,  0,  # 37
    8,   1,  0,  # 38
    4,   0, -2,  # 39
    6,   0,  2,  # 40
    3,   1,  0,  # 41
    5,   1,  0,  # 42
    13,  0,  0,  # 43
    20, -1,  0,  # 44
    3,   2,  0,  # 45
    4,  -2,  2,  # 46
    1,   2,  0,  # 47
    22, -1,  0,  # 48
    0,   0,  4,  # 49
    6,   0, -2,  # 50
    2,   1, -2,  # 51
    0,   2,  0,  # 52
    0,  -1,  2,  # 53
    2,   0,  4,  # 54
    0,  -2,  2,  # 55
    2,   2, -2,  # 56
    24,  0,  0,  # 57
    4,   0, -4,  # 58
    2,   2,  0,  # 59
    1,  -1,  0   # 60
  ), ncol = 3, byrow = TRUE)

  # Apogee multiplication coefficients.
  #   Column 1: d. (d = the Moon's mean elongation, in radians.)
  #   Column 2: m. (m = the Sun's mean anomaly, in radians.)
  #   Column 3: f. (f = Moon's argument of latitude, in radians.)
  apogee_matrix <- matrix(c(
    2,   0,  0,  #  1
    4,   0,  0,  #  2
    0,   1,  0,  #  3
    2,  -1,  0,  #  4
    0,   0,  2,  #  5
    1,   0,  0,  #  6
    6,   0,  0,  #  7
    4,  -1,  0,  #  8
    2,   0,  2,  #  9
    1,   1,  0,  # 10
    8,   0,  0,  # 11
    6,  -1,  0,  # 12
    2,   0, -2,  # 13
    2,  -2,  0,  # 14
    3,   0,  0,  # 15
    4,   0,  2,  # 16
    8,  -1,  0,  # 17
    4,  -2,  0,  # 18
    10,  0,  0,  # 19
    3,   1,  0,  # 20
    0,   2,  0,  # 21
    2,   1,  0,  # 22
    2,   2,  0,  # 23
    6,   0,  2,  # 24
    6,  -2,  0,  # 25
    10, -1,  0,  # 26
    5,   0,  0,  # 27
    4,   0, -2,  # 28
    0,   1,  2,  # 29
    12,  0,  0,  # 30
    2,   1,  2,  # 31
    1,  -1,  0   # 32
  ), ncol = 3, byrow = TRUE)

  # Body of: get_lunar_apsis_jde.
  function(date_target) {

    # Calculate lunar apsis dates.
    #   See Astronomical Algorithms, Meeus, chapter 50.
    #
    # Input:
    #   date_target    Date should be a Dynamical Time, as a date/time value.

    # Determine new moon coefficient "k".
    #   k = 0 corresponds to the perigee of 22nd Dec, 1999.
    #   Negative values of k give lunar apsides before that perigee.
    #   - For perigee, round to the nearest integer.
    #   - For apogee, round to the nearest integer + 0.5.
    #   Assumes target apsis is a global variable.
    #
    #   N.B. Ensure x.5 values are rounded towards zero if negative, but away
    #        from zero if positive.
    #   N.B. Round odd and even numbers consistently (R doesn't, by default).
    #   N.B. Convert target dynamical time to decimal years before use.

    k <- (mpfr(lubridate::decimal_date(date_target), 128L) - c_1999_97_mpfr) *
      c_13_2555_mpfr

    if (g_lunar_apsis == "perigee") {
      k <- floor(c_0_5_mpfr + k)
      is_perigee <- TRUE

    } else if (g_lunar_apsis == "apogee") {
      k <- c_0_5_mpfr + floor(k)
      is_perigee <- FALSE

    } else {
      cat("\nError in function get_lunar_apsis_jde...")
      cat("\nInvalid lunar apsis; g_lunar_apsis:\n")
      print(g_lunar_apsis)
      stopifnot(FALSE)
    }

    # Determine the approximate time "t" in Julian centuries from epoch 2000.
    t <- k / c_1325_55_mpfr
    t_mul_t <- t * t

    # Determine the mean perigee or apogee, in Julian Ephemeris Days "JDE".
    # Specify the calculation using Horner's Method.
    mean_apsis <- t_mul_t * (
      c_neg_0_0006691_mpfr + t * (t * c_0_0000000052_mpfr - c_0_000001098_mpfr)
    ) + c_27_55454989_mpfr * k + c_2451534_6698_mpfr

    # Calculate the Moon's mean elongation, in radians.
    d <- angle_degrees_to_radians(
      t_mul_t * (
        t * (t * c_0_000000055_mpfr - c_0_00001156_mpfr) - c_0_0100383_mpfr
      ) + c_171_9179_mpfr + c_335_9106046_mpfr * k
    )

    # Calculate the Sun's mean anomaly, in radians.
    m <- angle_degrees_to_radians(
      t_mul_t * (-t * c_0_0000010_mpfr - c_0_0008130_mpfr) + c_347_3477_mpfr +
        c_27_1577721_mpfr * k
    )

    # Calculate the Moon's argument of latitude, in radians.
    f <- angle_degrees_to_radians(
      t_mul_t * (-t * c_0_0000148_mpfr - c_0_0125053_mpfr) + c_316_6109_mpfr +
        c_364_5287911_mpfr * k
    )

    # Calculate the first set of periodic corrections, in days.
    if (is_perigee) {
      # Set apsis-specific uncertainty from Astronomical Algorithms, Meeus,
      # chapter 50.
      # N.B. The worst-case uncertainty, not mean uncertainty, is taken here.
      # N.B. Values here are minutes converted into JDE days, using 1 day = 1440
      #      minutes.
      apsis_uncertainty <- c_apsis_uncertainty_perigee

      # Calculate timing periodic corrections for perigee as a matrix-vector
      # product.
      mvp <- d * perigee_matrix[, 1] + m * perigee_matrix[, 2] +
        f * perigee_matrix[, 3]
      corrections_periodic <- sum(
        (c_vec_perigee_coeff + c_vec_perigee_subtrahend * t) * sin(mvp)
      )

    } else {
      # Set apsis-specific uncertainty.
      # N.B. Values here are minutes converted into JDE days, using 1 day = 1440
      #      minutes.
      apsis_uncertainty <- c_apsis_uncertainty_apogee

      # Calculate timing periodic corrections for apogee as a matrix-vector
      # product.
      mvp <- d * apogee_matrix[, 1] + m * apogee_matrix[, 2] +
        f * apogee_matrix[, 3]
      corrections_periodic <- sum(
        (c_vec_apogee_coeff + c_vec_apogee_subtrahend * t) *
          sin(mvp)
      )
    }

    # Diagnostic information.
    if (g_verbosity) {
      cat("   get_lunar_apsis_jde settings...\n", sep = "")
      cat("\n   k:")
      print(k, digits = 20)
      cat("\n   t:")
      print(t, digits = 20)
      cat("\n   d (degrees):")
      print(angle_radians_to_degrees(d), digits = 20)
      cat("\n   m (degrees):")
      print(angle_radians_to_degrees(m), digits = 20)
      cat("\n   f (degrees):")
      print(angle_radians_to_degrees(f), digits = 20)
      cat("\n   mean_apsis: ")
      print(mean_apsis, digits = 20)
      cat("\n   corrections_periodic: ")
      print(corrections_periodic, digits = 20)
    }

    list(
      "value" = mean_apsis + corrections_periodic,
      "uncertainty" = apsis_uncertainty
    )
  }
}

# Mid-chunk pause.
# Add a two-second pause to ameliorate RStudio input/output and RAM spikes.
Sys.sleep(2)

make_get_lunar_distance <- function() {

  # Closure constants.
  #
  # Lunar distance periodic calculation coefficients. 46 elements.
  c_vec_ld_coeff <- mpfr(c(
    "-20905355", "-3699111", "-2955968", "-569925", "48888",   "-3149",
    "246158",    "-152138",  "-170733",  "-204586", "-129620", "108743",
    "104755",    "10321",    "79661",    "-34782",  "-23210",  "-21636",
    "24208",     "30824",    "-8379",    "-16675",  "-12831",  "-10445",
    "-11650",    "14403",    "-7003",    "10056",   "6322",    "-9884",
    "5751",      "-4950",    "4130",     "-3958",   "3258",    "2616",
    "-1897",     "-2117",    "2354",     "-1423",   "-1117",   "-1571",
    "-1739",     "-4421",    "1165",     "8752"
  ), 128L)
  #
  # Lunar distance masks for vectorised operations. 46 elements each.
  c_mask_distance_d <- c(
    0, 2, 2, 0, 0, 0, 2, 2, 2, 2, 0, 1, 0, 2, 0, 4,
    0, 4, 2, 2, 1, 1, 2, 2, 4, 2, 0, 2, 1, 2, 0, 2,
    2, 4, 3, 2, 4, 0, 2, 4, 0, 4, 1, 0, 0, 2
  )
  c_mask_distance_f_mul_2 <- c(
    0,  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, -1, -1, 0,
    0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0,  0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1
  )
  c_mask_distance_m <- c(
    0,  0, 0, 0,  1, 0,  0, -1, 0, -1, 1,  0, 1,  0, 0,  0,
    0,  0, 1, 1,  0, 1, -1,  0, 0,  0, 1, -1, 0, -2, 1, -2,
    0, -1, 0, 1, -1, 2,  2,  0, 0, -1, 0,  0, 2,  0
  )
  c_mask_distance_m_prime <- c(
    1, -1,  0, 2,  0,  0, -2, -1, 1,  0, -1,  0, 1,  0, 1, -1,
    3, -2, -1, 0, -1,  0,  1,  2, 0, -3, -2, -2, 1,  0, 2, -1,
    1, -1, -1, 1, -2, -1, -1,  1, 4,  0, -2,  2, 1, -1
  )

  # Periodic distance term vector and masks. 46 elements.
  #
  # Indices of elements initialised to 1:
  #   1, 2, 3, 4, 6, 7, 9, 12, 14, 15, 16, 17, 18, 21, 24, 25, 26, 29, 33, 35,
  #   40, 41, 43, 44, 46.
  #
  # Indices of elements holding e:
  #   5, 8, 10, 11, 13, 19, 20, 22, 23, 27, 28, 31, 34, 36, 37, 42.
  #
  # Indices of elements holding e_mul_e:
  #   30, 32, 38, 39, 45.
  #
  # Index reference:
  #   1 1         2 1         3 1         4 1         5 e         6 1
  #   7 1         8 e         9 1        10 e        11 e        12 1
  #  13 e        14 1        15 1        16 1        17 1        18 1
  #  19 e        20 e        21 1        22 e        23 e        24 1
  #  25 1        26 1        27 e        28 e        29 1        30 e_mul_e
  #  31 e        32 e_mul_e  33 1        34 e        35 1        36 e
  #  37 e        38 e_mul_e  39 e_mul_e  40 1        41 1        42 e
  #  43 1        44 1        45 e_mul_e  46 1
  #
  c_vec_e_mask_e <- c(
    5, 8, 10, 11, 13, 19, 20, 22, 23, 27, 28, 31, 34, 36, 37, 42
  )
  c_vec_e_mask_e_mul_e <- c(30, 32, 38, 39, 45)

  c_0_mpfr              <- mpfr("0", 128L)
  c_0_0001536_mpfr      <- mpfr("0.0001536", 128L)
  c_0_001_mpfr          <- mpfr("0.001", 128L)
  c_0_0018819_mpfr      <- mpfr("0.0018819", 128L)
  c_0_0036539_mpfr      <- mpfr("0.0036539", 128L)
  c_0_0087414_mpfr      <- mpfr("0.0087414", 128L)
  c_0_5_mpfr            <- mpfr("0.5", 128L)
  c_2_mpfr              <- mpfr("2", 128L)
  c_93_272095_mpfr      <- mpfr("93.272095", 128L)
  c_134_9633964_mpfr    <- mpfr("134.9633964", 128L)
  c_297_8501921_mpfr    <- mpfr("297.8501921", 128L)
  c_357_5291092_mpfr    <- mpfr("357.5291092", 128L)
  c_1000_mpfr           <- mpfr("1000", 128L)
  c_35999_0502909_mpfr  <- mpfr("35999.0502909", 128L)
  c_36525_mpfr          <- mpfr("36525", 128L)
  c_385000_56_mpfr      <- mpfr("385000.56", 128L)
  c_445267_1114034_mpfr <- mpfr("445267.1114034", 128L)
  c_477198_8675055_mpfr <- mpfr("477198.8675055", 128L)
  c_483202_0175233_mpfr <- mpfr("483202.0175233", 128L)
  c_2451545_mpfr        <- mpfr("2451545", 128L)
  c_14712000_mpfr       <- mpfr("14712000", 128L)
  c_24490000_mpfr       <- mpfr("24490000", 128L)
  c_113065000_mpfr      <- mpfr("113065000", 128L)
  c_863310000_mpfr      <- mpfr("863310000", 128L)
  #
  c_lunar_argument_latitude <- mpfr("1", 128L) / mpfr("3526000", 128L)
  c_lunar_mean_anomaly      <- mpfr("1", 128L) / mpfr("69699", 128L)
  c_lunar_mean_elongation   <- mpfr("1", 128L) / mpfr("545868", 128L)

  # Body of: get_lunar_distance.
  function(jde_mpfr) {

    # Get the lunar distance, in kilometres, for a given JDE (Julian Ephemeris
    # Day). See Astronomical Algorithms, Meeus, chapter 47.
    #
    # Input:
    #  jde_mpfr    Target date, in Julian Ephemeris Days, best passed in mpfr
    #              format with, e.g. 128-bit accuracy.
    #              E.g. R stores a passed value of 2443259.9 as
    #                   2443259.8999999999 by default.
    #                   This easily causes errors of +- 1 second.

    if (!is.mpfr(jde_mpfr)) {
      jde_mpfr <- mpfr(jde_mpfr, 128L)
    }

    # Input must be >= 0.
    if (jde_mpfr < c_0_mpfr) {
      cat("\nError in function get_lunar_distance...")
      cat("\nInput dynamical time; jde_mpfr:\n")
      print(jde_mpfr, digits = 20)
      cat("\nNeed dynamical time >= 0.\n")
      stopifnot(FALSE)  # Halts program execution due to invalid input.
    }

    # Express the target date as "t" in Julian centuries from epoch 2000.
    # N.B. See Astronomical Algorithms, Meeus, chapter 22, p. 143.
    t <- (jde_mpfr - c_2451545_mpfr) / c_36525_mpfr

    # Calculate e, to correct for the (changing) eccentricity of the Earth's
    # orbit around the Sun.
    e <- get_eccentric_anomaly(t) # t will be a mpfr value here.

    # Maintain periodic distance data locally, not as a closure variable,
    # for thread-safety (important for vectorisation or parallelisation).
    # N.B. Update the data vector in-place.
    vec_e <- mpfr(rep("1", 46), 128L)
    vec_e[c_vec_e_mask_e] <- e
    vec_e[c_vec_e_mask_e_mul_e] <- e * e

    # Calculate the Moon's mean elongation, in radians.
    d <- angle_degrees_to_radians(
      c_297_8501921_mpfr +
        t * (
          c_445267_1114034_mpfr -
            t * (c_0_0018819_mpfr - t * (
              c_lunar_mean_elongation - t / c_113065000_mpfr
            ))
        )
    )

    # Calculate the Sun's mean anomaly, in radians.
    m <- angle_degrees_to_radians(
      c_357_5291092_mpfr +
        t * (
          c_35999_0502909_mpfr - t * (c_0_0001536_mpfr - t / c_24490000_mpfr)
        )
    )

    # Calculate the Moon's mean anomaly, in radians.
    m_prime <- angle_degrees_to_radians(
      c_134_9633964_mpfr +
        t * (
          c_477198_8675055_mpfr +
            t * (
              c_0_0087414_mpfr +
                t * (c_lunar_mean_anomaly - t / c_14712000_mpfr)
            )
        )
    )

    # Calculate the Moon's argument of latitude, in radians.
    f_mul_2 <- c_2_mpfr * angle_degrees_to_radians(
      c_93_272095_mpfr +
        t * (
          c_483202_0175233_mpfr -
            t * (
              c_0_0036539_mpfr +
                t * (c_lunar_argument_latitude - t / c_863310000_mpfr)
            )
        )
    )

    # Calculate lunar distance.
    scalar_multiples <- d * c_mask_distance_d +
      f_mul_2 * c_mask_distance_f_mul_2 +
      m * c_mask_distance_m +
      m_prime * c_mask_distance_m_prime
    #
    distance_value <- c_385000_56_mpfr +
      sum(c_vec_ld_coeff * vec_e * cos(scalar_multiples)) * c_0_001_mpfr

    # Diagnostic information.
    if (g_verbosity) {
      cat("   get_lunar_distance...\n", sep = "")
      cat("\n   t:")
      print(t, digits = 20)
      cat("\n   e:")
      print(e, digits = 20)
      cat("\n   d (degrees):")
      print(angle_radians_to_degrees(d), digits = 20)
      cat("\n   m (degrees):")
      print(angle_radians_to_degrees(m), digits = 20)
      cat("\n   m_prime (degrees):")
      print(angle_radians_to_degrees(m_prime), digits = 20)
      cat("\n   F (degrees):")
      print(angle_radians_to_degrees(f_mul_2 * c_0_5_mpfr), digits = 20)
      cat("\n   distance_periodic: ")
      print((distance_value - c_385000_56_mpfr) * c_1000_mpfr, digits = 20)
    }

    # Determine uncertainty by comparison to the complete ELP-2000/82 theory.
    #
    # ELP-2000/82 distance is the length derived from its three Cartesian
    # coordinate vectors.
    #
    # N.B. sqrt(sum(vec^2)) = sqrt(vec[1]^2 + vec[2]^2 + vec[3]^2)
    #      but is generally faster.
    #
    distance_uncertainty <- mpfr(sqrt(sum(call_elp82b_2(jde_mpfr) ^ 2)), 128L)

    list(
      "value" = distance_value,
      "uncertainty" = abs(distance_uncertainty - distance_value)
    )
  }
}

make_get_lunar_phase_jde <- function() {

  # Closure constants.

  # Full or new moon: periodic corrections term vector and masks.
  # Set A. 25 elements.
  #
  # Indices of elements initialised to 1:
  #   1, 3, 4, 8, 9, 11, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25.
  #
  # Indices of elements holding e:
  #   2, 5, 6, 10, 12, 13, 14.
  #
  # Indices of elements holding e_mul_e:
  #   7.
  #
  # Index reference:
  #   1 1   2 e   3 1   4 1   5 e   6 e   7 e_mul_e   8 1   9 1  10 e
  #  11 1  12 e  13 e  14 e  15 1  16 1  17 1        18 1  19 1  20 1
  #  21 1  22 1  23 1  24 1  25 1
  #
  c_vec_ph_e_mask_e <- c(2, 5, 6, 10, 12, 13, 14)


  # Quarter moon: periodic corrections term vector and masks.
  # Set A. 25 elements.
  #
  # Indices of elements initialised to 1:
  #   1, 4, 5, 8, 9, 10, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25.
  #
  # Indices of elements holding e:
  #   2, 3, 6, 11, 12, 13, 15.
  #
  # Indices of elements holding e_mul_e:
  #   7, 14.
  #
  # Index reference:
  #   1 1   2 e   3 e   4 1         5 1   6 e   7 e_mul_e   8 1   9 1  10 1
  #  11 e  12 e  13 e  14 e_mul_e  15 e  16 1  17 1        18 1  19 1  20 1
  #  21 1  22 1  23 1  24 1        25 1
  #
  c_vec_ph_qrt_e_mask_e <- c(2, 3, 6, 11, 12, 13, 15)
  c_vec_ph_qrt_e_mask_e_mul_e <- c(7, 14)


  # Supplemental quarter phase corrections A: term vector and mask. 5 elements.
  #
  # Indices of elements initialised to 1:
  #   2, 3, 4, 5.
  #
  # Indices of elements holding e:
  #   1.
  #
  # Index reference:
  #   1 e  2 1  3 1  4 1  5 1


  # Lunar phase (full moon) calculation coefficients.
  c_vec_ph_full_coeff <- mpfr(c(
    "-0.40614",         "0.17302",         "0.01614",  "0.01043",
    "0.00734",          "-0.00515",        "0.00209",  "-0.00111",
    "-0.00057",         "0.00056",         "-0.00042", "0.00042",
    "0.00038",          "-0.00024",        "-0.00017", "-0.00007",
    rep("0.00004", 2),  rep("0.00003", 2), "-0.00003", "0.00003",
    rep("-0.00002", 2), "0.00002"
  ), 128L)

  # Lunar phase (new moon) calculation coefficients.
  c_vec_ph_new_coeff <- mpfr(c(
    "-0.40720",         "0.17241",         "0.01608",  "0.01039",
    "0.00739",          "-0.00514",        "0.00208",  "-0.00111",
    "-0.00057",         "0.00056",         "-0.00042", "0.00042",
    "0.00038",          "-0.00024",        "-0.00017", "-0.00007",
    rep("0.00004", 2),  rep("0.00003", 2), "-0.00003", "0.00003",
    rep("-0.00002", 2), "0.00002"
  ), 128L)

  # Lunar phase (quarter phase) calculation coefficients.
  c_vec_ph_qrt_coeff <- mpfr(c(
    "-0.62801",        "0.17172",  "-0.01183", "0.00862",
    "0.00804",         "0.00454",  "0.00204",  "-0.00180",
    "-0.00070",        "-0.00040", "-0.00034", rep("0.00032", 2),
    "-0.00028",        "0.00027",  "-0.00017", "-0.00005",
    "0.00004",         "-0.00004", "0.00004",  rep("0.00003", 2),
    rep("0.00002", 2), "-0.00002"
  ), 128L)

  # Lunar phase (quarter phase) supplemental corrections.
  c_vec_ph_sup_coeff <- mpfr(c(
    "-0.00038", " 0.00026", "-0.00002", rep("0.00002", 2)
  ), 128L)

  # Lunar phase additional corrections coefficients.
  c_vec_ph_add_coeff_outer <- mpfr(c(
    "0.000325", "0.000165", "0.000164", "0.000126", "0.000110", "0.000062",
    "0.000060", "0.000056", "0.000047", "0.000042", "0.000040", "0.000037",
    "0.000035", "0.000023"
  ), 128L)
  #
  c_vec_ph_add_addend <- mpfr(c(
    "299.77", "251.88", "251.83", "349.42", "84.66",  "141.74", "207.14",
    "154.84", "34.52",  "207.19", "291.34", "161.72", "239.56", "331.55"
  ), 128L)
  #
  c_vec_ph_add_coeff_inner <- mpfr(c(
    "0.107408",  "0.016321",  "26.651886", "36.412478", "18.206239",
    "53.303771", "2.453732",  "7.306860",  "27.261239", "0.121824",
    "1.844379",  "24.198154", "25.513099", "3.592518"
  ), 128L)
  #
  c_vec_ph_add_subtrahend <- mpfr(c("0.009173", rep("0", 13)), 128L)

  c_0_00000000073_mpfr <- mpfr("0.00000000073", 128L)
  c_0_000000011_mpfr   <- mpfr("0.000000011", 128L)
  c_0_000000058_mpfr   <- mpfr("0.000000058", 128L)
  c_0_00000011_mpfr    <- mpfr("0.00000011", 128L)
  c_0_00000015_mpfr    <- mpfr("0.00000015", 128L)
  c_0_0000014_mpfr     <- mpfr("0.0000014", 128L)
  c_0_00000215_mpfr    <- mpfr("0.00000215", 128L)
  c_0_00000227_mpfr    <- mpfr("0.00000227", 128L)
  c_0_00001238_mpfr    <- mpfr("0.00001238", 128L)
  c_0_00015437_mpfr    <- mpfr("0.00015437", 128L)
  c_0_0016118_mpfr     <- mpfr("0.0016118", 128L)
  c_0_0020672_mpfr     <- mpfr("0.0020672", 128L)
  c_0_00306_mpfr       <- mpfr("0.00306", 128L)
  c_0_0107582_mpfr     <- mpfr("0.0107582", 128L)
  c_0_25_mpfr          <- mpfr("0.25", 128L)
  c_0_5_mpfr           <- mpfr("0.5", 128L)
  c_0_75_mpfr          <- mpfr("0.75", 128L)
  c_1_56375588_mpfr    <- mpfr("1.56375588", 128L)
  c_2_mpfr             <- mpfr("2", 128L)
  c_2_5534_mpfr        <- mpfr("2.5534", 128L)
  c_12_3685_mpfr       <- mpfr("12.3685", 128L)
  c_29_1053567_mpfr    <- mpfr("29.1053567", 128L)
  c_29_530588861_mpfr  <- mpfr("29.530588861", 128L)
  c_124_7746_mpfr      <- mpfr("124.7746", 128L)
  c_160_7108_mpfr      <- mpfr("160.7108", 128L)
  c_201_5643_mpfr      <- mpfr("201.5643", 128L)
  c_385_81693528_mpfr  <- mpfr("385.81693528", 128L)
  c_390_67050284_mpfr  <- mpfr("390.67050284", 128L)
  c_1236_85_mpfr       <- mpfr("1236.85", 128L)
  c_2000_mpfr          <- mpfr("2000", 128L)
  c_2451550_09766_mpfr <- mpfr("2451550.09766", 128L)
  #
  c_phase_uncertainty_first_qrt <- mpfr("15.3", 128L) / mpfr("86400", 128L)
  c_phase_uncertainty_full_moon <- mpfr("17.4", 128L) / mpfr("86400", 128L)
  c_phase_uncertainty_last_qrt <- mpfr("13", 128L) / mpfr("86400", 128L)
  c_phase_uncertainty_new_moon <- mpfr("16.4", 128L) / mpfr("86400", 128L)

  # Syzygy multiplication coefficients. (Syzygy is a full or new moon.)
  #   Column 1: m'.    (m' = the Moon's mean anomaly, in radians.)
  #   Column 2: m.     (m = the Sun's mean anomaly, in radians.)
  #   Column 3: f.     (f = the Moon's argument of latitude, in radians.)
  #   Column 4: omega. (omega = longitude of ascending node of the lunar orbit,
  #                      in radians.)
  syzygy_matrix <- matrix(c(
    1,  0,  0, 0,  #  1
    0,  1,  0, 0,  #  2
    2,  0,  0, 0,  #  3
    0,  0,  1, 0,  #  4
    1, -1,  0, 0,  #  5
    1,  1,  0, 0,  #  6
    0,  2,  0, 0,  #  7
    1,  0, -1, 0,  #  8
    1,  0,  1, 0,  #  9
    2,  1,  0, 0,  # 10
    3,  0,  0, 0,  # 11
    0,  1,  1, 0,  # 12
    0,  1, -1, 0,  # 13
    2, -1,  0, 0,  # 14
    0,  0,  0, 1,  # 15
    1,  2,  0, 0,  # 16
    2,  0, -1, 0,  # 17
    0,  3,  0, 0,  # 18
    1,  1, -1, 0,  # 19
    2,  0,  1, 0,  # 20
    1,  1,  1, 0,  # 21
    1, -1,  1, 0,  # 22
    1, -1, -1, 0,  # 23
    3,  1,  0, 0,  # 24
    4,  0,  0, 0   # 25
  ), ncol = 4, byrow = TRUE)

  # Quarter moon multiplication coefficients. (First or last quarter moon.)
  #   Column 1: m'.    (m' = the Moon's mean anomaly, in radians.)
  #   Column 2: m.     (m = the Sun's mean anomaly, in radians.)
  #   Column 3: f.     (f = the Moon's argument of latitude, in radians.)
  #   Column 4: omega. (omega = longitude of ascending node of the lunar orbit,
  #                      in radians.)
  quarter_matrix <- matrix(c(
    1,  0,  0, 0,  #  1
    0,  1,  0, 0,  #  2
    1,  1,  0, 0,  #  3
    2,  0,  0, 0,  #  4
    0,  0,  1, 0,  #  5
    1, -1,  0, 0,  #  6
    0,  2,  0, 0,  #  7
    1,  0, -1, 0,  #  8
    1,  0,  1, 0,  #  9
    3,  0,  0, 0,  # 10
    2, -1,  0, 0,  # 11
    0,  1,  1, 0,  # 12
    0,  1, -1, 0,  # 13
    1,  2,  0, 0,  # 14
    2,  1,  0, 0,  # 15
    0,  0,  0, 1,  # 16
    1, -1, -1, 0,  # 17
    2,  0,  1, 0,  # 18
    1,  1,  1, 0,  # 19
    1, -2,  0, 0,  # 20
    1,  1, -1, 0,  # 21
    0,  3,  0, 0,  # 22
    2,  0, -1, 0,  # 23
    1, -1,  1, 0,  # 24
    3,  1,  0, 0   # 25
  ), ncol = 4, byrow = TRUE)

  # Body of: get_lunar_phase_jde.
  function(date_target) {

    # Calculate lunar phase dates.
    #   See Astronomical Algorithms, Meeus, chapter 49.
    #
    # Input:
    #   date_target    Date should be a Dynamical Time, as a date/time value.

    # Determine new moon coefficient "k".
    #   k = 0 corresponds to the New Moon of 6th Jan, 2000.
    #   Negative values of k give lunar phases before the year 2000.
    #   - If new moon, round to the nearest integer.
    #   - If first quarter moon, round to the nearest integer + 0.25.
    #   - If full moon, round to the nearest integer + 0.5.
    #   - If last quarter moon, round to the nearest integer + 0.75.
    #   Assumes target phase is a global variable.
    #
    #   N.B. Ensure x.5 values are rounded towards zero if negative, but away
    #        from zero if positive.
    #   N.B. Round odd and even numbers consistently (R doesn't, by default).
    #   N.B. Convert target dynamical time to decimal years before use.

    # Set phase-specific uncertainty from Astronomical Algorithms, Meeus,
    # chapter 49, p. 354.
    # N.B. The worst-case uncertainty, not mean uncertainty, is taken here.
    # N.B. Values here are seconds converted into JDE days,
    #      using 1 day = 86,400 seconds.
    #
    k <- (mpfr(lubridate::decimal_date(date_target), 128L) - c_2000_mpfr) *
      c_12_3685_mpfr
    if (g_lunar_phase == "full moon") {
      k <- c_0_5_mpfr + floor(k)
      phase_uncertainty <- c_phase_uncertainty_full_moon
      is_full_moon <- TRUE
      is_new_moon <- FALSE
      # N.B. For a full moon, "is_first_quarter" isn't used.

    } else if (g_lunar_phase == "new moon") {
      k <- floor(c_0_5_mpfr + k)
      phase_uncertainty <- c_phase_uncertainty_new_moon
      is_full_moon <- FALSE
      is_new_moon <- TRUE
      # N.B. For a new moon, "is_first_quarter" isn't used.

    } else if (g_lunar_phase == "first quarter") {
      k <- sign(k) * c_0_25_mpfr + floor(c_0_5_mpfr - sign(k) * c_0_25_mpfr + k)
      phase_uncertainty <- c_phase_uncertainty_first_qrt
      is_full_moon <- FALSE
      is_new_moon <- FALSE
      is_first_quarter <- TRUE

    } else if (g_lunar_phase == "last quarter") {
      k <- sign(k) * c_0_75_mpfr + floor(c_0_5_mpfr - sign(k) * c_0_75_mpfr + k)
      phase_uncertainty <- c_phase_uncertainty_last_qrt
      is_full_moon <- FALSE
      is_new_moon <- FALSE
      is_first_quarter <- FALSE

    } else {
      cat("\nError in function get_lunar_phase_jde...")
      cat("\nInvalid lunar phase; g_lunar_phase:\n")
      print(g_lunar_phase)
      stopifnot(FALSE)
    }

    # Determine the approximate time "t" in Julian centuries from epoch 2000.
    t <- k / c_1236_85_mpfr
    t_mul_t <- t * t

    # Determine the mean lunar phase, in Julian Ephemeris Days "JDE".
    # Specify the calculation using Horner's Method.
    mean_phase <- t_mul_t * (
      c_0_00015437_mpfr + t * (t * c_0_00000000073_mpfr - c_0_00000015_mpfr)
    ) + c_29_530588861_mpfr * k + c_2451550_09766_mpfr

    # Calculate e, to correct for the (changing) eccentricity of the Earth's
    # orbit around the Sun.
    e <- get_eccentric_anomaly(t) # t will be a mpfr value here.

    # Calculate the Sun's mean anomaly, in radians.
    m <- angle_degrees_to_radians(
      t_mul_t * (-t * c_0_00000011_mpfr - c_0_0000014_mpfr) + c_2_5534_mpfr +
        c_29_1053567_mpfr * k
    )

    # Calculate the Moon's mean anomaly, in radians.
    m_prime <- angle_degrees_to_radians(
      t_mul_t * (
        -t_mul_t * c_0_000000058_mpfr + t * c_0_00001238_mpfr + c_0_0107582_mpfr
      ) + c_201_5643_mpfr + c_385_81693528_mpfr * k
    )

    # Calculate the Moon's argument of latitude, in radians.
    f_mul_2 <- c_2_mpfr * angle_degrees_to_radians(
      t_mul_t * (
        t * (t * c_0_000000011_mpfr - c_0_00000227_mpfr) - c_0_0016118_mpfr
      ) + c_160_7108_mpfr + c_390_67050284_mpfr * k
    )

    # Calculate the longitude of ascending node of the lunar orbit, in radians.
    omega <- angle_degrees_to_radians(
      t_mul_t * (t * c_0_00000215_mpfr + c_0_0020672_mpfr) + c_124_7746_mpfr -
        c_1_56375588_mpfr * k
    )

    # Calculate periodic corrections, in days.
    if (isTRUE(is_full_moon) || isTRUE(is_new_moon)) {

      # Maintain phase data locally, not as a closure variable,
      # for thread-safety (important for vectorisation or parallelisation).
      # N.B. Update the data vector in-place.
      vec_ph_e <- mpfr(rep("1", 25), 128L)
      vec_ph_e[c_vec_ph_e_mask_e] <- e
      vec_ph_e[7] <- e * e   # Element 7 stores e * e.

      # Calculate periodic corrections for syzygy as a matrix-vector product.
      mvp <- m_prime * syzygy_matrix[, 1] + m * syzygy_matrix[, 2] +
        f_mul_2 * syzygy_matrix[, 3] + omega * syzygy_matrix[, 4]

      if (isTRUE(is_full_moon)) {
        # A full moon is being processed.
        corrections_periodic <- sum(
          c_vec_ph_full_coeff * vec_ph_e * sin(mvp)
        )
      } else {
        # Assume a new moon is being processed.
        corrections_periodic <- sum(
          c_vec_ph_new_coeff * vec_ph_e * sin(mvp)
        )
      }

    } else {
      # Assume a quarter phase is being processed.

      # Maintain phase data locally, not as a closure variable,
      # for thread-safety (important for vectorisation or parallelisation).
      # N.B. Update the data vector in-place.
      vec_ph_qrt_e <- mpfr(rep("1", 25), 128L)
      vec_ph_qrt_e[c_vec_ph_qrt_e_mask_e] <- e
      vec_ph_qrt_e[c_vec_ph_qrt_e_mask_e_mul_e] <- e * e

      # Calculate periodic corrections for quarter phases as a matrix-vector
      # product.
      mvp <- m_prime * quarter_matrix[, 1] + m * quarter_matrix[, 2] +
        f_mul_2 * quarter_matrix[, 3] + omega * quarter_matrix[, 4]
      corrections_periodic <- sum(c_vec_ph_qrt_coeff * vec_ph_qrt_e * sin(mvp))

      # Supplemental quarter phase corrections.
      # Modify vectors in place.
      #
      vec_ph_sup_e <- mpfr(rep("1", 5), 128L)
      vec_ph_sup_e[1] <- e   # Element 1 stores e.
      #
      scalar_multiples <- m * c(1, 0, -1, 1, 0) + m_prime * c(0, 1, 1, 1, 0) +
        f_mul_2 * c(0, 0, 0, 0, 1)

      w <- c_0_00306_mpfr + sum(
        c_vec_ph_sup_coeff * vec_ph_sup_e * cos(scalar_multiples)
      )

      # Apply supplemental quarter corrections.
      #   - For the first quarter, add "w" to the periodic corrections.
      #   - For the last quarter, subtract "w" from the periodic corrections.
      corrections_periodic <- if (isTRUE(is_first_quarter)) {
        corrections_periodic + w
      } else {
        corrections_periodic - w
      }

      # Diagnostic information.
      if (g_verbosity) {
        cat("   get_lunar_phase_jde settings (quarter-specific)...\n", sep = "")
        cat("\n   w:")
        print(w, digits = 20)
      }
    }

    # Calculate additional corrections, based on "planetary arguments" in days.
    corrections_additional <- sum(
      c_vec_ph_add_coeff_outer * sin(
        angle_degrees_to_radians(
          c_vec_ph_add_addend + c_vec_ph_add_coeff_inner * k -
            c_vec_ph_add_subtrahend * t_mul_t
        )
      )
    )

    # Diagnostic information.
    if (g_verbosity) {
      cat("   get_lunar_phase_jde settings...\n", sep = "")
      cat("\n   k:")
      print(k, digits = 20)
      cat("\n   t:")
      print(t, digits = 20)
      cat("\n   e:")
      print(e, digits = 20)
      cat("\n   mean_phase: ")
      print(mean_phase, digits = 20)
      cat("\n   corrections_periodic: ")
      print(corrections_periodic, digits = 20)
      cat("\n   corrections_additional: ")
      print(corrections_additional, digits = 20)
    }

    list(
      "value" = corrections_additional + corrections_periodic + mean_phase,
      "uncertainty" = phase_uncertainty
    )
  }
}

# Mid-chunk pause.
# Add a two-second pause to ameliorate RStudio input/output and RAM spikes.
Sys.sleep(2)

make_get_utc <- function() {

  # Closure constants.

  # 31536000 seconds = 365 days.
  c_365_day_scaling <- mpfr("1", 128L) / mpfr("31536000", 128L)
  # 31622400 seconds = 366 days.
  c_366_day_scaling <- mpfr("1", 128L) / mpfr("31622400", 128L)

  c_0_5_mpfr <- mpfr("0.5", 128L)
  c_5_minutes <- lubridate::dminutes(mpfr("5", 128L))

  # N.B. Assume 01-Jan-2026 00:00:00.000 UTC = 01-01-2026 00:01:12.150 DT.
  #      ~2026.000002287871666339920000000000000000 DT.
  # N.B. Technically, here, 01-Jan-2026 00:00:00.000 UTC =
  #      01-Jan-2026 00:01:12.150 +- 00:00:03.261 Dynamical Time.
  c_threshold_2026 <- as.POSIXct(
    "01-01-2026 00:01:12.150", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"
  )
  # Assume 01-Jan-2050 00:00:00.000 UTC = 01-Jan-2050 00:01:33.000 DT.
  c_threshold_2050 <- as.POSIXct(
    "01-01-2050 00:01:33.000", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"
  )

  # Body of: get_utc.
  function(date_target) {

    # Estimate UTC from an input Dynamical Time.
    # Effectively reverses the operation in get_del_t, except it returns UTC,
    # not UT1.
    #
    # Input:
    #   date_target    Date should be a Dynamical Time.
    #                  It's assumed the recovered UTC date will be in the 21st
    #                  century.

    # Ensure a valid date has been entered.
    # Known del T values exist until 31-Dec-2025 23:59:59.
    # UTC should be specified exactly for such dates, not calculated here.
    #
    # Before 2026...
    if (date_target < c_threshold_2026) {
      cat("\nError in function get_utc...")
      cat(
        "\nInput date:  ", format(date_target, "%d-%b-%Y %H:%M:%OS"),
        " (assumed Dynamical Time)"
      )
      cat("\nNeed date >= 01-01-2026 00:01:12.150 Dynamical Time", sep = "")
      stopifnot(FALSE)
    }

    dt_decimal <- mpfr(lubridate::decimal_date(date_target), 128L)
    small_diff <- mpfr(lubridate::decimal_date(date_target - c_5_minutes), 128L)

    # Determine how many seconds were in the UTC year, and get a inverse scaling
    # factor.
    #
    # N.B. Remember, del T is a function of the UTC year (not DT year).
    #
    #        Dynamical Time (decimal year),
    #          is equal to: UTC (decimal year) + (ay^2 + by + c) (seconds)
    #
    if (is_utc_leap_year(date_target, dt_decimal)) {
      scaling <- c_366_day_scaling
    } else {
      scaling <- c_365_day_scaling
    }

    # Before 2050...
    if (date_target < c_threshold_2050) {
      f1 <- function(x) {
        # N.B. Input "x" is always already a mpfr value by here.
        x + scaling * get_del_t_pre_2050(x) - dt_decimal
      }

    } else { # After 2050...
      f1 <- function(x) {
        x + scaling * get_del_t_post_2050(x) - dt_decimal
      }
    }

    # Return the UTC value as a date-time.
    #
    #   N.B. By default unirootR warns if the root finding calculation fails to
    #        converge.
    #
    #   N.B. unirootR may take either a vector "interval" argument, or "lower"
    #        and "upper" arguments. Avoid passing mpfr values to "interval"
    #        via the c() operator, as they will be downcast and lose precision.
    #
    #   N.B. Using tol = 1e-20 is generally fine, but the R numerical mantissa
    #        was used here for added rigour.

    utc <- unirootR(
      f1, lower = small_diff, upper = dt_decimal, verbose = FALSE, tol = 1e-36
    )

    utc_decimal <- lubridate::date_decimal(asNumeric(utc$root), tz = "UTC")

    list(
      "value" = utc_decimal,
      "uncertainty" = difftime(
        lubridate::date_decimal(
          asNumeric(utc$root + utc$estim.prec * c_0_5_mpfr), tz = "UTC"
        ),
        utc_decimal, tz = "UTC", units = "auto"
      )
    )
  }
}

make_is_utc_leap_year <- function() {

  # Closure constants.

  # Leap year check transition constants.
  c_vec_leap_year_transition <- mpfr(c(
    "2029.000002440489652144610000000000000000",
    "2033.000002525305490053139999999999999999",
    "2037.000002615792482174580000000000000000",
    "2041.000002711950855882609999999999999999",
    "2045.000002813780383803530000000000000003",
    "2049.000002921281065937360000000000000003",
    "2053.000003143493358948040000000000000003",
    "2057.000003405644292797659999999999999996",
    "2061.000003671042577479970000000000000000",
    "2065.000003939687758247599999999999999994",
    "2069.000004211580289847920000000000000000",
    "2073.000004486719717533559999999999999995",
    "2077.000004765106496051890000000000000002",
    "2081.000005046740170655539999999999999999",
    "2085.000005331621196091869999999999999999",
    "2089.000005619748662866189999999999999998",
    "2093.000005911123935220530000000000000000",
    "2097.000006205745648912850000000000000004"
  ), 128L)

  # Leap year check sequence.
  c_seq_leap_year_check <- seq(2029, 2097, by = 4)

  # Body of: is_utc_leap_year.
  function(date_target, date_target_decimal_mpfr) {

    # Test if a Dynamical Time falls in a UTC leap year.
    # This information is used to convert DT to UTC date-times > 31-Dec-2025
    # 23:59:59.999 (UTC).
    #
    # Input:
    #   date_target          Date should be a UTC time.
    #   date_target_decimal  Decimal equivalent of date_target, as a mpfr value.
    #
    #     N.B. An input date equivalent to a 21st century UTC date is assumed.
    #     N.B. date_target should be >= 01-Jan-2026 UTC.
    #          I.e. >= 01-Jan-2026 00:01:12.150 Dynamical Time.
    #          I.e. >= 2026.000002287871666339920000000000000000 DT.
    #
    # Output:
    #   A Boolean indicating if the Dynamical Time falls in a UTC leap year.

    # Get the year of the input Dynamical Time.
    year_target <- lubridate::year(date_target)

    # Allow for Dynamical Times falling shortly after New Year's Eve of a UTC
    # leap year.
    #
    #   E.g. 2024 is a leap year.
    #   Assume 01-Jan-2025 00:00:00.000 UTC -->
    #            01-Jan-2025 00:01:14.467 +- 4.561s Dynamical Time.
    #            ~2025.000002361344968448980000000000000003.
    #
    #   So 2025 Dynamical Time dates < 2025.000002361344968448980000000000000003
    #   actually started in 2024 UTC.
    #

    # Get the index of the input year in the vector c_seq_leap_year_check, if
    # it's there. N.B. R indices start from 1, not 0.
    x <- match(year_target, c_seq_leap_year_check, nomatch = 0)
    if (x > 0) {
      if (date_target_decimal_mpfr < c_vec_leap_year_transition[x]) {
        year_target <- year_target - 1
      }
    }

    lubridate::leap_year(year_target)
  }
}

make_time_dynam_to_jde <- function() {

  # Closure constants.
  c_0_01_mpfr    <- mpfr("0.01", 128L)
  c_0_25_mpfr    <- mpfr("0.25", 128L)
  c_1_mpfr       <- mpfr("1", 128L)
  c_2_mpfr       <- mpfr("2", 128L)
  c_3_mpfr       <- mpfr("3", 128L)
  c_12_mpfr      <- mpfr("12", 128L)
  c_24_mpfr      <- mpfr("24", 128L)
  c_30_6001_mpfr <- mpfr("30.6001", 128L)
  c_365_25_mpfr  <- mpfr("365.25", 128L)
  c_1440_mpfr    <- mpfr("1440", 128L)
  c_1524_5_mpfr  <- mpfr("1524.5", 128L)
  c_4716_mpfr    <- mpfr("4716", 128L)
  c_86400_mpfr   <- mpfr("86400", 128L)

  # Body of: time_dynamical_to_jde.
  function(date_target) {

    # Convert Dynamical Time to JDE (Julian Ephemeris Day).
    #   See Astronomical Algorithms, Meeus, chapter 7.
    #   N.B. Assumes a Gregorian calendar.
    #
    # Input:
    #   date_target    A date-time, in Dynamical Time format.

    yr <- mpfr(lubridate::year(date_target), 128L)
    mn <- mpfr(lubridate::month(date_target), 128L)
    if (mn < c_3_mpfr) {
      yr <- yr - c_1_mpfr
      mn <- mn + c_12_mpfr
    }

    # Get dy as a decimal day.
    #
    #   Minutes per day = 1440.
    #   Seconds per day = 86400.
    #   Sum elements from smallest to largest to reduce roundoff.
    #
    dy <- mpfr(lubridate::second(date_target), 128L) / c_86400_mpfr + # Seconds.
      mpfr(lubridate::minute(date_target), 128L) / c_1440_mpfr +      # Minutes.
      mpfr(lubridate::hour(date_target), 128L) / c_24_mpfr +          # Hours.
      mpfr(lubridate::day(date_target), 128L)                         # Days.

    a <- floor(yr * c_0_01_mpfr)
    b <- c_2_mpfr + floor(a * c_0_25_mpfr) - a

    b + dy + floor(c_30_6001_mpfr * (mn + c_1_mpfr)) - c_1524_5_mpfr +
      floor(c_365_25_mpfr * (yr + c_4716_mpfr))
  }
}

make_time_dynam_to_universal <- function() {

  # Closure constants.

  c_zero_seconds  <- lubridate::dseconds(mpfr("0", 128L))

  c_offsets <- c(
    lubridate::dseconds(mpfr("43.18", 128L)),  #  1 (42.18 + 1 sec)
    lubridate::dseconds(mpfr("44.18", 128L)),  #  2 (42.18 + 2 sec)
    lubridate::dseconds(mpfr("45.18", 128L)),  #  3 (42.18 + 3 sec)
    lubridate::dseconds(mpfr("46.18", 128L)),  #  4 (42.18 + 4 sec)
    lubridate::dseconds(mpfr("47.18", 128L)),  #  5 (42.18 + 5 sec)
    lubridate::dseconds(mpfr("48.18", 128L)),  #  6 (42.18 + 6 sec)
    lubridate::dseconds(mpfr("49.18", 128L)),  #  7 (42.18 + 7 sec)
    lubridate::dseconds(mpfr("50.18", 128L)),  #  8 (42.18 + 8 sec)
    lubridate::dseconds(mpfr("51.18", 128L)),  #  9 (42.18 + 9 sec)
    lubridate::dseconds(mpfr("52.18", 128L)),  # 10 (42.18 + 10 sec)
    lubridate::dseconds(mpfr("53.18", 128L)),  # 11 (42.18 + 11 sec)
    lubridate::dseconds(mpfr("54.18", 128L)),  # 12 (42.18 + 12 sec)
    lubridate::dseconds(mpfr("55.18", 128L)),  # 13 (42.18 + 13 sec)
    lubridate::dseconds(mpfr("56.18", 128L)),  # 14 (42.18 + 14 sec)
    lubridate::dseconds(mpfr("57.18", 128L)),  # 15 (42.18 + 15 sec)
    lubridate::dseconds(mpfr("58.18", 128L)),  # 16 (42.18 + 16 sec)
    lubridate::dseconds(mpfr("59.18", 128L)),  # 17 (42.18 + 17 sec)
    lubridate::dseconds(mpfr("60.18", 128L)),  # 18 (42.18 + 18 sec)
    lubridate::dseconds(mpfr("61.18", 128L)),  # 19 (42.18 + 19 sec)
    lubridate::dseconds(mpfr("62.18", 128L)),  # 20 (42.18 + 20 sec)
    lubridate::dseconds(mpfr("63.18", 128L)),  # 21 (42.18 + 21 sec)
    lubridate::dseconds(mpfr("64.18", 128L)),  # 22 (42.18 + 22 sec)
    lubridate::dseconds(mpfr("65.18", 128L)),  # 23 (42.18 + 23 sec)
    lubridate::dseconds(mpfr("66.18", 128L)),  # 24 (42.18 + 24 sec)
    lubridate::dseconds(mpfr("67.18", 128L)),  # 25 (42.18 + 25 sec)
    lubridate::dseconds(mpfr("68.18", 128L)),  # 26 (42.18 + 26 sec)
    lubridate::dseconds(mpfr("69.18", 128L)),  # 27 (42.18 + 27 sec)
    lubridate::dseconds(mpfr("42.18", 128L))   # 28 Put the default value last.
  )

  # Thresholds correspond to the (irregular) addition of leap seconds, as
  # managed by the International Earth Rotation and Reference Systems Service.
  c_thresholds <- as.POSIXct(c(
    "01-07-1972 00:00:42.18", "01-01-1973 00:00:43.18",  #  1,  2
    "01-01-1974 00:00:44.18", "01-01-1975 00:00:45.18",  #  3,  4
    "01-01-1976 00:00:46.18", "01-01-1977 00:00:47.18",  #  5,  6
    "01-01-1978 00:00:48.18", "01-01-1979 00:00:49.18",  #  7,  8
    "01-01-1980 00:00:50.18", "01-07-1981 00:00:51.18",  #  9, 10
    "01-07-1982 00:00:52.18", "01-07-1983 00:00:53.18",  # 11, 12
    "01-07-1985 00:00:54.18", "01-01-1988 00:00:55.18",  # 13, 14
    "01-01-1990 00:00:56.18", "01-01-1991 00:00:57.18",  # 15, 16
    "01-07-1992 00:00:58.18", "01-07-1993 00:00:59.18",  # 17, 18
    "01-07-1994 00:01:00.18", "01-01-1996 00:01:01.18",  # 19, 20
    "01-07-1997 00:01:02.18", "01-01-1999 00:01:03.18",  # 21, 22
    "01-01-2006 00:01:04.18", "01-01-2009 00:01:05.18",  # 23, 24
    "01-07-2012 00:01:06.18", "01-07-2015 00:01:07.18",  # 25, 26
    "01-01-2017 00:01:08.18", "01-01-2026 00:01:09.18"   # 27, 28
  ), format = "%d-%m-%Y %H:%M:%OS",  tz = "UTC")

  # Body of: time_dynamical_to_universal.
  function(date_target)   {

    # Convert Dynamical Time to Coordinated Universal Time.
    #
    #   N.B. This function reverses the process in time_universal_to_dynamical.
    #   N.B. It's assumed a 21st century date will be processed.
    #        I.e. 01-Jan-2001 to 31-Dec-2100.
    #        Limited support for other dates was provided, for testing purposes.
    #
    # Input:
    #   date_target    A date-time, in Dynamical Time format.
    # Output:
    #   UTC, where UTC = date_target - offset.

    # These values are given in UTC, to avoid accidental assignment to other
    # time zones. If the UTC label was omitted, they would be valid dynamical
    # time values.

    if (date_target >= c_thresholds[28]) {
      # Get UTC via algorithmic estimation.
      #   N.B. get_utc returns a date-time.
      #        get_utc returns uncertainty as a time difference.
      return(get_utc(date_target))
    }

    # In R, findInterval(test_val, vector_intervals), returns the index of the
    #   interval in vector_interval, for which test_val is greater than.
    #
    # If test_val is < than every interval, it returns 0.
    # If test_val is >= interval n, but < interval (n + 1), it returns n.
    # If test_val is > every interval, it returns length(vector_intervals).
    #
    #  N.B. It assumes data in vector_intervals are sorted small to large.

    date_idx <- findInterval(date_target, c_thresholds)
    offset_match <- if (date_idx == 0) c_offsets[28] else c_offsets[date_idx]
    list("value" = date_target - offset_match, "uncertainty" = c_zero_seconds)
  }
}

make_time_jde_to_dynam <- function() {

  # Closure constants.
  c_0_mpfr          <- mpfr("0", 128L)
  c_0_25_mpfr       <- mpfr("0.25", 128L)
  c_0_5_mpfr        <- mpfr("0.5", 128L)
  c_1_mpfr          <- mpfr("1", 128L)
  c_2_mpfr          <- mpfr("2", 128L)
  c_13_mpfr         <- mpfr("13", 128L)
  c_14_mpfr         <- mpfr("14", 128L)
  c_30_6001_mpfr    <- mpfr("30.6001", 128L)
  c_122_1_mpfr      <- mpfr("122.1", 128L)
  c_365_25_mpfr     <- mpfr("365.25", 128L)
  c_1524_mpfr       <- mpfr("1524", 128L)
  c_4715_mpfr       <- mpfr("4715", 128L)
  c_4716_mpfr       <- mpfr("4716", 128L)
  c_36524_25_mpfr   <- mpfr("36524.25", 128L)
  c_1867216_25_mpfr <- mpfr("1867216.25", 128L)
  c_2299161_mpfr    <- mpfr("2299161", 128L)

  # Body of: time_jde_to_dynamical.
  function(jde_mpfr) {

    # Convert JDE (Julian Ephemeris Day) to Dynamical Time.
    #   See Astronomical Algorithms, Meeus, chapter 7.
    #   N.B. Error on p. 63: "2291 161" should be "2299 161".
    #
    # Input:
    #   jde_mpfr    Target date, in Julian Ephemeris Days, best passed in mpfr
    #               format with, e.g. 128-bit accuracy.
    #               E.g. R stores a passed value of 2443259.9 as
    #                    2443259.8999999999 by default.
    #                    This easily causes errors of +- 1 second.

    if (!is.mpfr(jde_mpfr)) {
      jde_mpfr <- mpfr(jde_mpfr, 128L)
    }

    # Input must be >= 0.
    if (jde_mpfr < c_0_mpfr) {
      cat("\nError in function time_jde_to_dynamical.\nInput dynamical time:\n")
      print(jde_mpfr, digits = 20)
      cat("\nNeed dynamical time >= 0.\n")
      stopifnot(FALSE)
    }

    jde <- c_0_5_mpfr + jde_mpfr
    z <- trunc(jde)   # Integer part of the input JDE.
    f <- jde - z      # Fractional part of the input JDE.

    if (z < c_2299161_mpfr) {
      a <- z
    } else {
      alpha <- floor((z - c_1867216_25_mpfr) / c_36524_25_mpfr)
      a <- c_1_mpfr + alpha - floor(alpha * c_0_25_mpfr) + z
    }

    b <- c_1524_mpfr + a
    c <- floor((b - c_122_1_mpfr) / c_365_25_mpfr)
    d <- floor(c_365_25_mpfr * c)
    e <- floor((b - d) / c_30_6001_mpfr)

    # Day value, with time, as a decimal.
    day <- b - d - floor(c_30_6001_mpfr * e) + f
    # Day value, without time.
    day_without_time <- trunc(day)

    # Month, as an integer.
    if (e < c_14_mpfr) {
      month <- e - c_1_mpfr
    } else {
      month <- e - c_13_mpfr
    }

    # Year, as an integer.
    if (month > c_2_mpfr) {
      year <- c - c_4716_mpfr
    } else {
      year <- c - c_4715_mpfr
    }

    # The Dynamical Time returned here is stored in UTC date-time format.
    # If the time zone label is neglected however, the date and time are the
    # correct Dynamical Time.
    #
    # N.B. Any fractional part of the day gets added to the output.
    #
    # N.B. If no time zone is specified in a lubridate::make_datetime call, it
    #      defaults to UTC.
    #
    as.POSIXct(
      lubridate::make_datetime(
        year = as.integer(year),
        month = as.integer(month),
        day = as.integer(day_without_time)
      ) + lubridate::ddays(day - day_without_time),
      format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"
    )
  }
}

make_time_universal_to_dynam <- function() {

  # Closure constants.

  c_900_milliseconds <- lubridate::dmilliseconds(mpfr("900", 128L))

  c_zero_seconds <- lubridate::dseconds(mpfr("0", 128L))

  c_offsets <- c(
    lubridate::dseconds(mpfr("43.18", 128L)),  #  1 (42.18 + 1 sec)
    lubridate::dseconds(mpfr("44.18", 128L)),  #  2 (42.18 + 2 sec)
    lubridate::dseconds(mpfr("45.18", 128L)),  #  3 (42.18 + 3 sec)
    lubridate::dseconds(mpfr("46.18", 128L)),  #  4 (42.18 + 4 sec)
    lubridate::dseconds(mpfr("47.18", 128L)),  #  5 (42.18 + 5 sec)
    lubridate::dseconds(mpfr("48.18", 128L)),  #  6 (42.18 + 6 sec)
    lubridate::dseconds(mpfr("49.18", 128L)),  #  7 (42.18 + 7 sec)
    lubridate::dseconds(mpfr("50.18", 128L)),  #  8 (42.18 + 8 sec)
    lubridate::dseconds(mpfr("51.18", 128L)),  #  9 (42.18 + 9 sec)
    lubridate::dseconds(mpfr("52.18", 128L)),  # 10 (42.18 + 10 sec)
    lubridate::dseconds(mpfr("53.18", 128L)),  # 11 (42.18 + 11 sec)
    lubridate::dseconds(mpfr("54.18", 128L)),  # 12 (42.18 + 12 sec)
    lubridate::dseconds(mpfr("55.18", 128L)),  # 13 (42.18 + 13 sec)
    lubridate::dseconds(mpfr("56.18", 128L)),  # 14 (42.18 + 14 sec)
    lubridate::dseconds(mpfr("57.18", 128L)),  # 15 (42.18 + 15 sec)
    lubridate::dseconds(mpfr("58.18", 128L)),  # 16 (42.18 + 16 sec)
    lubridate::dseconds(mpfr("59.18", 128L)),  # 17 (42.18 + 17 sec)
    lubridate::dseconds(mpfr("60.18", 128L)),  # 18 (42.18 + 18 sec)
    lubridate::dseconds(mpfr("61.18", 128L)),  # 19 (42.18 + 19 sec)
    lubridate::dseconds(mpfr("62.18", 128L)),  # 20 (42.18 + 20 sec)
    lubridate::dseconds(mpfr("63.18", 128L)),  # 21 (42.18 + 21 sec)
    lubridate::dseconds(mpfr("64.18", 128L)),  # 22 (42.18 + 22 sec)
    lubridate::dseconds(mpfr("65.18", 128L)),  # 23 (42.18 + 23 sec)
    lubridate::dseconds(mpfr("66.18", 128L)),  # 24 (42.18 + 24 sec)
    lubridate::dseconds(mpfr("67.18", 128L)),  # 25 (42.18 + 25 sec)
    lubridate::dseconds(mpfr("68.18", 128L)),  # 26 (42.18 + 26 sec)
    lubridate::dseconds(mpfr("69.18", 128L)),  # 27 (42.18 + 27 sec)
    lubridate::dseconds(mpfr("42.18", 128L))   # 28 Put the default value last.
  )

  # Thresholds correspond to the (irregular) addition of leap seconds, as
  # managed by the International Earth Rotation and Reference Systems Service.
  c_thresholds_timeless <- as.POSIXct(c(
    "01-07-1972 00:00:00", "01-01-1973 00:00:00",  #  1,  2
    "01-01-1974 00:00:00", "01-01-1975 00:00:00",  #  3,  4
    "01-01-1976 00:00:00", "01-01-1977 00:00:00",  #  5,  6
    "01-01-1978 00:00:00", "01-01-1979 00:00:00",  #  7,  8
    "01-01-1980 00:00:00", "01-07-1981 00:00:00",  #  9, 10
    "01-07-1982 00:00:00", "01-07-1983 00:00:00",  # 11, 12
    "01-07-1985 00:00:00", "01-01-1988 00:00:00",  # 13, 14
    "01-01-1990 00:00:00", "01-01-1991 00:00:00",  # 15, 16
    "01-07-1992 00:00:00", "01-07-1993 00:00:00",  # 17, 18
    "01-07-1994 00:00:00", "01-01-1996 00:00:00",  # 19, 20
    "01-07-1997 00:00:00", "01-01-1999 00:00:00",  # 21, 22
    "01-01-2006 00:00:00", "01-01-2009 00:00:00",  # 23, 24
    "01-07-2012 00:00:00", "01-07-2015 00:00:00",  # 25, 26
    "01-01-2017 00:00:00", "01-01-2026 00:00:00"   # 27, 28
  ), format = "%d-%m-%Y %H:%M:%OS",  tz = "UTC")

  # Body of: time_universal_to_dynamical.
  function(utc) {

    # Convert Coordinated Universal Time to Dynamical Time.
    #   N.B. It's assumed a 21st century date will be processed.
    #        I.e. 01-Jan-2001 to 31-Dec-2100.
    #        Limited support for other dates was provided for testing purposes.
    #
    # Input:
    #   utc    Target date, in Coordinated Universal Time.
    #
    # N.B. The output value is stored as a UTC date-time format.
    #      If the time zone label is neglected however, the date and time are
    #      the correct Dynamical Time.
    #
    #      DT = UTC + 42.18 s + # leap seconds
    #
    # Ref: Polynomial approximations to Delta T, 1620-2000 AD (Meeus & Simons).
    #
    # N.B. Dynamical time can be determined exactly for dates from 1972 to the
    #      present, estimated after.
    # N.B. Number of leap seconds added from 1972-2020 = 27.
    # N.B. Timings are available from Wikipedia, the Astronomical Almanac, etc.

    if (utc >= c_thresholds_timeless[28]) {
      # Get delta T via algorithmic estimation.
      # N.B. The algorithm finds delta T for a UT (or UT1).
      # N.B. Thus, allow for uncertainty in UTC versus UT1 (always <= 0.9 s; as
      #      maintained by leap seconds).
      del_t <- calculate_with_uncertainty(
        get_del_t, utc, c_900_milliseconds, TRUE
      )
      return(list(
        "value" = utc + lubridate::dseconds(del_t$value),
        "uncertainty" = lubridate::dseconds(del_t$uncertainty)
      ))
    }

    # In R, findInterval(test_val, vector_intervals), returns the index of the
    #   interval in vector_interval, for which test_val is greater than.
    #
    # If test_val is < than every interval, it returns 0.
    # If test_val is >= interval n, but < interval (n + 1), it returns n.
    # If test_val is > every interval, it returns length(vector_intervals).
    #
    #  N.B. It assumes data in vector_intervals are sorted small to large.

    date_idx <- findInterval(utc, c_thresholds_timeless)
    offset_match <- if (date_idx == 0) c_offsets[28] else c_offsets[date_idx]
    list("value" = utc - offset_match, "uncertainty" = c_zero_seconds)
  }
}

###############################################################################

# Pre closures instantiation pause.
# Add a two-second pause to ameliorate RStudio input/output and RAM spikes.
Sys.sleep(2)

# Implement closures.
angle_degrees_to_radians <- make_angle_degrees_to_radians()
angle_radians_to_degrees <- make_angle_radians_to_degrees()
get_del_t <- make_get_del_t()
get_del_t_post_2050 <- make_get_del_t_post_2050()
get_del_t_pre_2050 <- make_get_del_t_pre_2050()
get_eccentric_anomaly <- make_get_eccentric_anomaly()
get_lunar_apsis_jde <- make_get_lunar_apsis_jde()
get_lunar_distance <- make_get_lunar_distance()
get_lunar_phase_jde <- make_get_lunar_phase_jde()
get_utc <- make_get_utc()
is_utc_leap_year <- make_is_utc_leap_year()
time_dynamical_to_jde <- make_time_dynam_to_jde()
time_dynamical_to_universal <- make_time_dynam_to_universal()
time_jde_to_dynamical <- make_time_jde_to_dynam()
time_universal_to_dynamical <- make_time_universal_to_dynam()

# Post closures instantiation pause.
# Add a two-second pause to ameliorate RStudio input/output and RAM spikes.
Sys.sleep(2)
```

``` r
# Precision Testing
# -----------------
#
# Quantify the effect of using R's default double (53-bit) precision versus
# mpfr 128-bit precision. Here, converting a date-time in JDE (Julian Ephemeris
# Day) to dynamical time is the test case.

g_lunar_apsis <- "perigee"
g_lunar_phase <- "full moon"

jde_a <- 2443259.9
jde_b <- mpfr("2443259.9", 128L)

result_a <- time_jde_to_dynamical(jde_a)
result_b <- time_jde_to_dynamical(jde_b)

cat("Precision Testing\n-----------------\n\n", sep = "")

cat("JDE, target value...\n2443259.9\n\n", sep = "")
cat("JDE, stored as a double variable (low precision)...\n", sep = "")
print(jde_a, digits = 17)
cat("\nJDE, stored as a MPFR variable (high precision)...\n", sep = "")
print(jde_b, digits = 17)

cat(
  paste0(
    "\n\nConverting JDE to dynamical time, ",
    "low precision (top) versus high precision (bottom)...\n"
  ),
  sep = ""
)
print(result_a, digits = 17)
print(result_b, digits = 17)

cat("\nWhat's the corresponding difference in lunar distance?\n")
lunar_distance_a <- get_lunar_distance(jde_a)
lunar_distance_b <- get_lunar_distance(jde_b)

difference <- abs(lunar_distance_a$value - lunar_distance_b$value) *
  mpfr("1000000", 128L)
cat("\U2248", formatMpfr(difference, 1), " mm\n", sep = "")
```

    ## Precision Testing
    ## -----------------
    ## 
    ## JDE, target value...
    ## 2443259.9
    ## 
    ## JDE, stored as a double variable (low precision)...
    ## [1] 2443259.8999999999
    ## 
    ## JDE, stored as a MPFR variable (high precision)...
    ## 1 'mpfr' number of precision  128   bits 
    ## [1] 2443259.9
    ## 
    ## 
    ## Converting JDE to dynamical time, low precision (top) versus high precision (bottom)...
    ## [1] "1977-04-26 09:35:59.999991 UTC"
    ## [1] "1977-04-26 09:36:00 UTC"
    ## 
    ## What's the corresponding difference in lunar distance?
    ## ≈0.5 mm

``` r
# Testing
# -------

cat("Testing\n-------\n\n", sep = "")
cat(
  paste0(
    "Test the system's functions against examples in ",
    "'Astronomical Algorithms'.\n"
  )
)
cat(
  paste0(
    "Ref: Meeus J., Astronomical Algorithms, 2nd Ed, Willmann-Bell, ",
    "Richmond, VA, 1998.\n"
  )
)

######
cat("\n-------\nTesting: New Moon Feb 1977 (page 353)\n\n")
g_lunar_phase <- "new moon"
date_target_utc <- lubridate::date_decimal(1977.13)
date_target_dynamical <- time_universal_to_dynamical(date_target_utc)
phase_date <- calculate_with_uncertainty(
  get_lunar_phase_jde, date_target_dynamical$value,
  date_target_dynamical$uncertainty, TRUE
)
label <- paste0("Next ", g_lunar_phase, ":")
cat(
  label, replicate(g_max_indent - nchar(label), " "),
  formatMpfr(phase_date$value, 11),
  " +- ", formatMpfr(phase_date$uncertainty, 1), " JDE\n",
  sep = ""
)
phase_date_dynamical <- calculate_with_uncertainty(
  time_jde_to_dynamical, phase_date$value, phase_date$uncertainty, FALSE
)
cat(
  replicate(g_max_indent - 2, " "), "= ",
  format(
    lubridate::round_date(phase_date_dynamical$value, unit = "seconds"),
    "%d-%b-%Y %H:%M:%S"
  ),
  " +- ",
  format(
    hms::as_hms(
      signif(as.numeric(phase_date_dynamical$uncertainty, units = "secs"), 2)
    ),
    "%H:%M:%S"
  ),
  " Dynamical Time\n",
  sep = ""
)
phase_date_utc <- calculate_with_uncertainty(
  time_dynamical_to_universal, phase_date_dynamical$value,
  phase_date_dynamical$uncertainty, TRUE
)
cat(
  replicate(g_max_indent - 2, " "), "= ",
  format(
    lubridate::round_date(phase_date_utc$value, unit = "seconds"),
    "%d-%b-%Y %H:%M:%S"
  ),
  " +- ",
  format(
    hms::as_hms(
      signif(as.numeric(phase_date_utc$uncertainty, units = "secs"), 2)
    ),
    "%H:%M:%S"
  ),
  " UTC\n\n",
  sep = ""
)

label <- paste("Expected:")
cat(
  label, replicate(g_max_indent - nchar(label), " "), "18-Feb-1977 03:37 UT.\n",
  sep = ""
)
cat(replicate(g_max_indent, " "), "N.B. UT = UTC +- 0.9 s.\n", sep = "")
######


######
cat("\n-------\nTesting: Last Quarter Moon Jan 2044 (page 353)\n\n")
g_lunar_phase <- "last quarter"
date_target_utc <- lubridate::date_decimal(2044.026)
date_target_dynamical <- time_universal_to_dynamical(date_target_utc)
phase_date <- calculate_with_uncertainty(
  get_lunar_phase_jde, date_target_dynamical$value,
  date_target_dynamical$uncertainty, TRUE
)
label <- paste0("Next ", g_lunar_phase, ":")
cat(
  label, replicate(g_max_indent - nchar(label), " "),
  formatMpfr(phase_date$value, 11), " +- ",
  formatMpfr(phase_date$uncertainty, 1), " JDE\n",
  sep = ""
)
phase_date_dynamical <- calculate_with_uncertainty(
  time_jde_to_dynamical, phase_date$value, phase_date$uncertainty, FALSE
)
cat(
  replicate(g_max_indent - 2, " "), "= ",
  format(
    lubridate::round_date(phase_date_dynamical$value, unit = "seconds"),
    "%d-%b-%Y %H:%M:%S"
  ),
  " +- ",
  format(
    hms::as_hms(
      signif(as.numeric(phase_date_dynamical$uncertainty, units = "secs"), 2)
    ),
    "%H:%M:%S"
  ),
  " Dynamical Time\n\n",
  sep = ""
)

label <- paste("Expected:")
cat(
  label, replicate(g_max_indent - nchar(label), " "),
  "21-Jan-2044 23:48:17 Dynamical Time.\n",
  sep = ""
)
######


######
cat("\n-------\nTesting: Lunar Apogee Oct 1988 (page 357)\n\n")
g_lunar_apsis <- "apogee"
date_target_utc <- lubridate::date_decimal(1988.75)
date_target_dynamical <- time_universal_to_dynamical(date_target_utc)
apsis_date <- calculate_with_uncertainty(
  get_lunar_apsis_jde, date_target_dynamical$value,
  date_target_dynamical$uncertainty, TRUE
)
label <- paste0("Next ", g_lunar_apsis, ":")
cat(
  label, replicate(g_max_indent - nchar(label), " "),
  formatMpfr(apsis_date$value, 10), " +- ",
  formatMpfr(apsis_date$uncertainty, 1), " JDE\n",
  sep = ""
)
apsis_date_dynamical <- calculate_with_uncertainty(
  time_jde_to_dynamical, apsis_date$value, apsis_date$uncertainty, FALSE
)
cat(
  replicate(g_max_indent - 2, " "), "= ",
  format(
    lubridate::round_date(apsis_date_dynamical$value, unit = "seconds"),
    "%d-%b-%Y %H:%M:%S"
  ),
  " +- ",
  format(
    hms::as_hms(
      signif(as.numeric(apsis_date_dynamical$uncertainty, units = "secs"), 3)
    ),
    "%H:%M:%S"
  ),
  " Dynamical Time\n\n",
  sep = ""
)

label <- paste("Expected:")
cat(
  label, replicate(g_max_indent - nchar(label), " "),
  "07-Oct-1988 20:29 Dynamical Time.\n",
  sep = ""
)
######


######
cat(
  paste0(
    "\n-------\nTesting: ",
    "Lunar Distance for 12 April 1992, at 00:00:00 Dynamical Time ",
    "(pages 342-343).\n"
  )
)
lunar_distance <- calculate_with_uncertainty(
  get_lunar_distance, mpfr("2448724.5", 128L), mpfr("0", 128L), TRUE
)
label <- paste0("Estimated lunar distance:")
cat(
  "\n", label, replicate(g_max_indent - nchar(label), " "),
  signif(asNumeric(lunar_distance$value), 6), " +- ",
  signif(asNumeric(lunar_distance$uncertainty), 1), " km\n\n",
  sep = ""
)

label <- paste("Expected:")
cat(
  label, replicate(g_max_indent - nchar(label), " "), "368409.7 km.\n",
  sep = ""
)
######
```

    ## Testing
    ## -------
    ## 
    ## Test the system's functions against examples in 'Astronomical Algorithms'.
    ## Ref: Meeus J., Astronomical Algorithms, 2nd Ed, Willmann-Bell, Richmond, VA, 1998.
    ## 
    ## -------
    ## Testing: New Moon Feb 1977 (page 353)
    ## 
    ## Next new moon:             2443192.6512 +- 0.0002 JDE
    ##                          = 18-Feb-1977 03:37:42 +- 00:00:16 Dynamical Time
    ##                          = 18-Feb-1977 03:36:54 +- 00:00:16 UTC
    ## 
    ## Expected:                  18-Feb-1977 03:37 UT.
    ##                            N.B. UT = UTC +- 0.9 s.
    ## 
    ## -------
    ## Testing: Last Quarter Moon Jan 2044 (page 353)
    ## 
    ## Next last quarter:         2467636.4919 +- 0.0002 JDE
    ##                          = 21-Jan-2044 23:48:17 +- 00:00:13 Dynamical Time
    ## 
    ## Expected:                  21-Jan-2044 23:48:17 Dynamical Time.
    ## 
    ## -------
    ## Testing: Lunar Apogee Oct 1988 (page 357)
    ## 
    ## Next apogee:               2447442.354 +- 0.002 JDE
    ##                          = 07-Oct-1988 20:29:41 +- 00:03:00 Dynamical Time
    ## 
    ## Expected:                  07-Oct-1988 20:29 Dynamical Time.
    ## 
    ## -------
    ## Testing: Lunar Distance for 12 April 1992, at 00:00:00 Dynamical Time (pages 342-343).
    ## 
    ## Estimated lunar distance:  368410 +- 4 km
    ## 
    ## Expected:                  368409.7 km.

``` r
# Testing against grey literature examples.
# -----------------------------------------
#
# Compare estimated lunar distance estimations.
# Compare this system versus astropixels.com and fourmilab.com.

g_lunar_apsis <- "perigee"
g_lunar_phase <- "full moon"

label_date <- "Target date:"
chars_label_date <- nchar(label_date)
label_distance <- "Lunar distance"

vec_dates_astropixels <- c(
  as.POSIXct("09-01-2001 20:24:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"),
  as.POSIXct("04-07-2050 18:51:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"),
  as.POSIXct("23-05-2100 17:25:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC")
)
vec_dates_fourmilab <- c(
  as.POSIXct("10-01-2001 09:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"),
  as.POSIXct("07-07-2050 02:26:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"),
  as.POSIXct("22-05-2100 11:08:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC")
)

vec_distances_astropixels <- c(357406, 367058, 360904)
vec_distances_fourmilab   <- c(357131, 363255, 359497)

# Get lunar distance at target dates, for astropixels dates.
# N.B. Convert dates from UTC --> Dynamical Time--> JDE.
vec_calculated_distances <- lapply(
  lapply(
    lapply(vec_dates_astropixels, time_universal_to_dynamical),
    function(x) {
      calculate_with_uncertainty(
        time_dynamical_to_jde, x$value, x$uncertainty, FALSE
      )
    }
  ),
  function(x) {
    calculate_with_uncertainty(
      get_lunar_distance, x$value, x$uncertainty, TRUE
    )
  }
)

# Display astropixels.com data.
cat("Lunar distance benchmarking: astropixels.com\n\n")
for (i in seq_along(vec_dates_astropixels)) {

  # Display target date.
  cat(
    "------\n", label_date, replicate(g_max_indent - chars_label_date, " "),
    format(vec_dates_astropixels[i], "%d-%b-%Y %H:%M:%S %Z"), "\n",
    sep = ""
  )

  # Display astropixels.com's estimate.
  cat("\n", label_distance, sep = "")
  cat(
    "\nEstimate, astropixels.com:",
    replicate(g_max_indent - nchar("Estimate, astropixels.com:"), " "),
    vec_distances_astropixels[i], " +- <unspecified> km",
    sep = ""
  )

  # Display lunar distance at target date.
  cat(
    "\nEstimate, R:", replicate(g_max_indent - nchar("Estimate, R:"), " "),
    signif(asNumeric(vec_calculated_distances[[i]]$value), 6), " +- ",
    signif(asNumeric(vec_calculated_distances[[i]]$uncertainty), 1), " km\n\n",
    sep = ""
  )
}


# Get lunar distance at target dates, for fourmilab dates.
# N.B. That source states that its times are only accurate to +- 2 minutes.
#      That error is incorporated, below.
vec_calculated_distances <- lapply(
  lapply(
    lapply(
      vec_dates_fourmilab,
      function(x) {
        calculate_with_uncertainty(
          time_universal_to_dynamical, x, lubridate::dminutes(g_2_mpfr), TRUE
        )
      }
    ),
    function(x) {
      calculate_with_uncertainty(
        time_dynamical_to_jde, x$value, x$uncertainty, FALSE
      )
    }
  ),
  function(x) {
    calculate_with_uncertainty(get_lunar_distance, x$value, x$uncertainty, TRUE)
  }
)

# Display fourmilab.com data.
# N.B. The first test distance needs a different number of decimal places than
#      the rest.
cat(replicate(g_max_indent, "-"), sep = "")
cat("\n\nLunar distance benchmarking: www.fourmilab.ch\n\n")
#
# Display target date.
cat(
  "------\n", label_date, replicate(g_max_indent - chars_label_date, " "),
  format(vec_dates_fourmilab[1], "%d-%b-%Y %H:%M:%S"), " +- 00:00:02 UTC\n",
  sep = ""
)

# Display fourmilab.com's estimate.
cat("\n", label_distance, sep = "")
cat(
  "\nEstimate, fourmilab.com:",
  replicate(g_max_indent - nchar("Estimate, fourmilab.com:"), " "),
  format(vec_distances_fourmilab[1], nsmall = 1), " +- ",
  signif(asNumeric(vec_calculated_distances[[1]]$uncertainty), 1),
  " km (assumed uncertainty, from ELP 2000-82)",
  sep = ""
)

# Display lunar distance at target date.
cat(
  "\nEstimate, R:", replicate(g_max_indent - nchar("Estimate, R:"), " "),
  format(signif(asNumeric(vec_calculated_distances[[1]]$value), 7), nsmall = 1),
  " +- ", signif(asNumeric(vec_calculated_distances[[1]]$uncertainty), 1),
  " km\n\n",
  sep = ""
)
#
for (i in 2:length(vec_dates_fourmilab)) {

  # Display target date.
  cat(
    "------\n", label_date, replicate(g_max_indent - chars_label_date, " "),
    format(vec_dates_fourmilab[i], "%d-%b-%Y %H:%M:%S"), " +- 00:00:02 UTC\n",
    sep = ""
  )

  # Display fourmilab.com's estimate.
  cat("\n", label_distance, sep = "")
  cat(
    "\nEstimate, fourmilab.com:",
    replicate(g_max_indent - nchar("Estimate, fourmilab.com:"), " "),
    vec_distances_fourmilab[i], " +- ",
    signif(asNumeric(vec_calculated_distances[[i]]$uncertainty), 1),
    " km (assumed uncertainty, from ELP 2000-82)",
    sep = ""
  )

  # Display lunar distance at target date.
  cat(
    "\nEstimate, R:", replicate(g_max_indent - nchar("Estimate, R:"), " "),
    signif(asNumeric(vec_calculated_distances[[i]]$value), 6), " +- ",
    signif(asNumeric(vec_calculated_distances[[i]]$uncertainty), 1), " km\n\n",
    sep = ""
  )
}
```

    ## Lunar distance benchmarking: astropixels.com
    ## 
    ## ------
    ## Target date:               09-Jan-2001 20:24:00 UTC
    ## 
    ## Lunar distance
    ## Estimate, astropixels.com: 357406 +- <unspecified> km
    ## Estimate, R:               357408 +- 1 km
    ## 
    ## ------
    ## Target date:               04-Jul-2050 18:51:00 UTC
    ## 
    ## Lunar distance
    ## Estimate, astropixels.com: 367058 +- <unspecified> km
    ## Estimate, R:               367057 +- 4 km
    ## 
    ## ------
    ## Target date:               23-May-2100 17:25:00 UTC
    ## 
    ## Lunar distance
    ## Estimate, astropixels.com: 360904 +- <unspecified> km
    ## Estimate, R:               360904 +- 1 km
    ## 
    ## ---------------------------
    ## 
    ## Lunar distance benchmarking: www.fourmilab.ch
    ## 
    ## ------
    ## Target date:               10-Jan-2001 09:00:00 +- 00:00:02 UTC
    ## 
    ## Lunar distance
    ## Estimate, fourmilab.com:   357131.0 +- 0.8 km (assumed uncertainty, from ELP 2000-82)
    ## Estimate, R:               357131.0 +- 0.8 km
    ## 
    ## ------
    ## Target date:               07-Jul-2050 02:26:00 +- 00:00:02 UTC
    ## 
    ## Lunar distance
    ## Estimate, fourmilab.com:   363255 +- 2 km (assumed uncertainty, from ELP 2000-82)
    ## Estimate, R:               363255 +- 2 km
    ## 
    ## ------
    ## Target date:               22-May-2100 11:08:00 +- 00:00:02 UTC
    ## 
    ## Lunar distance
    ## Estimate, fourmilab.com:   359497 +- 3 km (assumed uncertainty, from ELP 2000-82)
    ## Estimate, R:               359493 +- 3 km

``` r
# Example 1
# ---------
#
# For a given date, determine the time of the next...
# a) full moon,
# b) perigee and associated lunar distance.

cat("Example 1\n---------\n\n", sep = "")
cat("For a target date, determine the time of the next....\n")
cat("a) full moon,\n")
cat("b) perigee and associated lunar distance.\n\n")

g_lunar_apsis <- "perigee"
g_lunar_phase <- "full moon"

# Set target date.
# Assume format is dd-mm-yyyy hh:mm:ss, with UTC time zone.
date_target_utc <- as.POSIXct(
  "18-02-2029 00:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"
)
label <- "Target date:"
cat(
  label, replicate(g_max_indent - nchar(label), " "),
  format(date_target_utc, "%d-%b-%Y %H:%M:%S %Z"), "\n",
  sep = ""
)

# Ensure a 21st century date is being processed.
if (!is_valid_date(date_target_utc)) {
  cat("\nInvalid date; date_target_utc:\n")
  print(date_target_utc, digits = 20)
  cat("\nNeed 21st century dates.")
  stopifnot(FALSE)
}

# Get the dynamical date equivalent of the target date.
date_target_dynamical <- time_universal_to_dynamical(date_target_utc)
cat(
  replicate(g_max_indent - 2, " "), "= ",
  format(
    lubridate::round_date(date_target_dynamical$value, unit = "seconds"),
    "%d-%b-%Y %H:%M:%S"
  ),
  " +- ", format(signif(date_target_dynamical$uncertainty, 1), "%H:%M:%S"),
  " Dynamical Time\n\n",
  sep = ""
)

# Get the date of the next nominated lunar phase (set to full moon by default).
phase_date <- calculate_with_uncertainty(
  get_lunar_phase_jde, date_target_dynamical$value,
  date_target_dynamical$uncertainty, TRUE
)
label <- paste0("Next ", g_lunar_phase, ":")
cat(
  label, replicate(g_max_indent - nchar(label), " "),
  format(signif(asNumeric(phase_date$value), 11), nsmall = 4), " +- ",
  signif(asNumeric(phase_date$uncertainty), 1), " JDE\n",
  sep = ""
)

# Convert JDE to dynamical time.
# N.B. Handle uncertainty differently than for target date; it used a decimal
#      uncertainty, this uses a difftime value.
phase_date_dynamical <- calculate_with_uncertainty(
  time_jde_to_dynamical, phase_date$value, phase_date$uncertainty, FALSE
)
cat(
  replicate(g_max_indent - 2, " "), "= ",
  format(
    lubridate::round_date(phase_date_dynamical$value, unit = "seconds"),
    "%d-%b-%Y %H:%M:%S"
  ),
  " +- ",
  format(
    hms::as_hms(
      signif(as.numeric(phase_date_dynamical$uncertainty, units = "secs"), 2)
    ),
    "%H:%M:%S"
  ),
  " Dynamical Time\n",
  sep = ""
)

# Convert dynamical time to UTC.
# N.B. As for the JDE to dynamical time conversion, the uncertainty here uses a
#      difftime value.
phase_date_utc <- calculate_with_uncertainty(
  time_dynamical_to_universal, phase_date_dynamical$value,
  phase_date_dynamical$uncertainty, TRUE
)
cat(
  replicate(g_max_indent - 2, " "), "= ",
  format(
    lubridate::round_date(phase_date_utc$value, unit = "seconds"),
    "%d-%b-%Y %H:%M:%S"
  ),
  " +- ",
  format(
    hms::as_hms(
      signif(as.numeric(phase_date_utc$uncertainty, units = "secs"), 2)
    ),
    "%H:%M:%S"
  ),
  " UTC\n",
  sep = ""
)


# Get the perigee estimation.
apsis_date <- calculate_with_uncertainty(
  get_lunar_apsis_jde, date_target_dynamical$value,
  date_target_dynamical$uncertainty, TRUE
)
label <- paste0("Next ", g_lunar_apsis, ":")
cat(
  "\n", label, replicate(g_max_indent - nchar(label), " "),
  format(signif(asNumeric(apsis_date$value), 9), nsmall = 2), " +- ",
  signif(asNumeric(apsis_date$uncertainty), 1), " JDE\n",
  sep = ""
)

# Convert JDE to dynamical time.
# N.B. Handle uncertainty differently than for target date; it used a decimal
#      uncertainty, this uses a difftime value.
apsis_date_dynamical <- calculate_with_uncertainty(
  time_jde_to_dynamical, apsis_date$value, apsis_date$uncertainty, FALSE
)
cat(
  replicate(g_max_indent - 2, " "), "= ",
  format(
    lubridate::round_date(apsis_date_dynamical$value, unit = "seconds"),
    "%d-%b-%Y %H:%M:%S"
  ),
  " +- ",
  format(
    hms::as_hms(
      signif(as.numeric(apsis_date_dynamical$uncertainty, units = "secs"), 4)
    ),
    "%H:%M:%S"
  ),
  " Dynamical Time\n",
  sep = ""
)

# Convert dynamical time to UTC.
# N.B. As for the JDE to dynamical time conversion, the uncertainty here uses a
#      difftime value.
apsis_date_utc <- calculate_with_uncertainty(
  time_dynamical_to_universal, apsis_date_dynamical$value,
  apsis_date_dynamical$uncertainty, TRUE
)
cat(
  replicate(g_max_indent - 2, " "), "= ",
  format(
    lubridate::round_date(apsis_date_utc$value, unit = "seconds"),
    "%d-%b-%Y %H:%M:%S"
  ),
  " +- ",
  format(
    hms::as_hms(
      signif(as.numeric(apsis_date_utc$uncertainty, units = "secs"), 4)
    ),
    "%H:%M:%S"
  ),
  " UTC\n",
  sep = ""
)


# Get lunar distance at target apsis.
lunar_distance <- calculate_with_uncertainty(
  get_lunar_distance, apsis_date$value, apsis_date$uncertainty, TRUE
)
label <- paste0("Lunar distance at ", g_lunar_apsis, ":")
cat(
  "\n", label, replicate(g_max_indent - nchar(label), " "),
  signif(asNumeric(lunar_distance$value), 6), " +- ",
  signif(asNumeric(lunar_distance$uncertainty), 1), " km\n",
  sep = ""
)
```

    ## Example 1
    ## ---------
    ## 
    ## For a target date, determine the time of the next....
    ## a) full moon,
    ## b) perigee and associated lunar distance.
    ## 
    ## Target date:               18-Feb-2029 00:00:00 UTC
    ##                          = 18-Feb-2029 00:01:17 +- 6s Dynamical Time
    ## 
    ## Next full moon:            2462196.2163 +- 2e-04 JDE
    ##                          = 28-Feb-2029 17:11:27 +- 00:00:17 Dynamical Time
    ##                          = 28-Feb-2029 17:10:10 +- 00:00:56 UTC
    ## 
    ## Next perigee:              2462197.27 +- 0.02 JDE
    ##                          = 01-Mar-2029 18:30:52 +- 00:31:00 Dynamical Time
    ##                          = 01-Mar-2029 18:29:35 +- 00:31:00 UTC
    ## 
    ## Lunar distance at perigee: 358630 +- 3 km

``` r
# Example 2a
# ----------
#
# Graph variation in lunar distance for 2019-2020.
# Overlay full and new moons.
# First step: generate required data.

cat("Example 2\n---------\n\n", sep = "")
cat("Graph variation in lunar distance for 2019-2020.\n")
cat("Overlay full and new moons.\n")

# Prepare poster data.
date_start <- as.POSIXct(
  "01-01-2019 00:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"
)
date_end <- as.POSIXct(
  "31-12-2020 00:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"
)

# Series 1: Get all full moons and their distances.
g_lunar_phase <- "full moon"
list_full_moon <- get_list_estimates(
  get_lunar_phase_jde, date_start, date_end, mpfr("29.4", 128L)
)
#
list_full_moon_distances <- lapply(
  list_full_moon,
  function(x) {
    calculate_with_uncertainty(get_lunar_distance, x$value, x$uncertainty, TRUE)
  }
)
#
# Convert full moon JDE dates to UTC.
list_full_moon <- lapply(
  list_full_moon,
  function(x) {
    target_dynamical <- calculate_with_uncertainty(
      time_jde_to_dynamical, x$value, x$uncertainty, FALSE
    )
    calculate_with_uncertainty(
      time_dynamical_to_universal, target_dynamical$value,
      target_dynamical$uncertainty, TRUE
    )
  }
)


# Series 2: Get all new moons and their distances.
g_lunar_phase <- "new moon"
list_new_moon <- get_list_estimates(
  get_lunar_phase_jde, date_start, date_end, mpfr("29.4", 128L)
)
#
list_new_moon_distances <- lapply(
  list_new_moon,
  function(x) {
    calculate_with_uncertainty(get_lunar_distance, x$value, x$uncertainty, TRUE)
  }
)
#
# Convert new moon JDE dates to UTC.
list_new_moon <- lapply(
  list_new_moon,
  function(x) {
    target_dynamical <- calculate_with_uncertainty(
      time_jde_to_dynamical, x$value, x$uncertainty, FALSE
    )
    calculate_with_uncertainty(
      time_dynamical_to_universal, target_dynamical$value,
      target_dynamical$uncertainty, TRUE
    )
  }
)


# Series 3: Get enough other lunar distances to graph.
# Try a point every 0.5 days.
interval_days <- lubridate::ddays(mpfr("0.5", 128L))
list_filler_dates <- 0:(ceiling((date_end - date_start) / interval_days) - 1)
list_filler_dates <- date_start + list_filler_dates * interval_days
list_filler_distances <- lapply(
  list_filler_dates,
  function(x) {
    target_dynamical <- time_universal_to_dynamical(x)
    target_jde <- calculate_with_uncertainty(
      time_dynamical_to_jde, target_dynamical$value,
      target_dynamical$uncertainty, FALSE
    )
    calculate_with_uncertainty(
      get_lunar_distance, target_jde$value, target_jde$uncertainty, TRUE
    )
  }
)


# Convert lists to dataframes.
#
lunar_data <- cbind(
  data.frame(
    # Filler distances (values and uncertainties).
    day = lubridate::as_datetime(unlist(list_filler_dates)),
    distance = unlist(lapply(
      list_filler_distances,
      function(x) {
        suppress_warnings(asNumeric(x$value), "*")
      }
    )),
    distance_uncertainty = unlist(lapply(
      list_filler_distances,
      function(x) {
        suppress_warnings(asNumeric(x$uncertainty), "*")
      }
    ))
  )
)

full_moon_data <- cbind(
  data.frame(
    # Full moon UTC dates (values and uncertainties).
    day_full_moon = lubridate::as_datetime(unlist(lapply(
      list_full_moon,
      function(x) {
        x$value
      }
    ))),
    day_full_moon_uncertainty = lubridate::as_datetime(unlist(lapply(
      list_full_moon,
      function(x) {
        x$uncertainty
      }
    ))),
    # Full moon distances (values and uncertainties).
    distance_full_moon = unlist(lapply(
      list_full_moon_distances,
      function(x) {
        suppress_warnings(asNumeric(x$value), "*")
      }
    )),
    distance_full_moon_uncertainty = unlist(lapply(
      list_full_moon_distances,
      function(x) {
        suppress_warnings(asNumeric(x$uncertainty), "*")
      }
    ))
  )
)

new_moon_data <- cbind(
  data.frame(
    # New moon UTC dates (values and uncertainties).
    day_new_moon = lubridate::as_datetime(unlist(lapply(
      list_new_moon,
      function(x) {
        x$value
      }
    ))),
    day_new_moon_uncertainty = lubridate::as_datetime(unlist(lapply(
      list_new_moon,
      function(x) {
        x$uncertainty
      }
    ))),
    # New moon distances (values and uncertainties).
    distance_new_moon = unlist(lapply(
      list_new_moon_distances,
      function(x) {
        suppress_warnings(asNumeric(x$value), "*")
      }
    )),
    distance_new_moon_uncertainty = unlist(lapply(
      list_new_moon_distances,
      function(x) {
        suppress_warnings(asNumeric(x$uncertainty), "*")
      }
    ))
  )
)

cat("\nData prepared.")
```

    ## Example 2
    ## ---------
    ## 
    ## Graph variation in lunar distance for 2019-2020.
    ## Overlay full and new moons.
    ## 
    ## Data prepared.

``` r
# Example 2b
# ----------
#
# Graph variation in lunar distance for 2019-2020.
# Overlay full and new moons.
# Second step: display poster-style graph.

# Ensure the extrafont database has been instantiated.
# N.B. Instantiation may take a few minutes if not previously run.
if (is.null(extrafont::fonts())) {
  extrafont::font_import(prompt = FALSE)
}

# Import Candara font, if not already available.
if (!("Candara" %in% extrafont::fonts())) {
  extrafont::font_import(pattern = "Candara.ttf")
}

# Register non-standard fonts with the R session.
extrafont::loadfonts(quiet = TRUE)

# Compose graph.
p <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data = lunar_data, ggplot2::aes(x = day, y = distance, color = "distance")
  ) +
  ggplot2::geom_point(
    data = full_moon_data,
    ggplot2::aes(
      x = day_full_moon, y = distance_full_moon, color = "full moon"
    ),
    shape = 16, size = 4
  ) +
  ggplot2::geom_point(
    data = new_moon_data,
    ggplot2::aes(x = day_new_moon, y = distance_new_moon, color = "new moon"),
    shape = 4, size = 4, stroke = 2
  ) +
  ggplot2::ggtitle(
    "Lunar Distance versus Date",
    subtitle = "(with full and new moons overlaid)"
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(shape = c(NA, 16, 4)))
  ) +
  ggplot2::scale_color_manual(
    values = rep("#18677E", 3), name = "Legend"
  ) +
  ggplot2::scale_x_datetime(labels = scales::date_format("%d-%b-%Y")) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::theme_bw() +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(color = "#000000"),
    axis.text.y = ggplot2::element_text(color = "#000000"),
    legend.background = ggplot2::element_rect(color = NA, fill = "transparent"),
    legend.box.background = ggplot2::element_rect(
      color = NA, fill = "transparent"
    ),
    panel.background = ggplot2::element_rect(fill = "transparent"),
    plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    text = ggplot2::element_text(
      color = "#717171", family = "Candara", size = 23.76097441
    )
  ) +
  ggplot2::xlab("\ndate (UTC)") +
  ggplot2::ylab("distance (km)\n")

# Save graph, large format for poster use.
ggplot2::ggsave(
  "supermoon_graph.png", bg = "transparent", width = 11.9, height = 8.499,
  units = "in"
)
```

``` r
# Example 2c
# ----------
#
# Graph variation in lunar distance for 2019-2020.
# Overlay full and new moons.
# Third step: display research paper-style graph.

# Compose graph.
p <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data = lunar_data, ggplot2::aes(x = day, y = distance, color = "distance")
  ) +
  ggplot2::geom_point(
    data = full_moon_data,
    ggplot2::aes(
      x = day_full_moon, y = distance_full_moon, color = "full moon"
    ),
    shape = 16
  ) +
  ggplot2::geom_point(
    data = new_moon_data,
    ggplot2::aes(x = day_new_moon, y = distance_new_moon, color = "new moon"),
    shape = 4
  ) +
  ggplot2::ggtitle(
    "Lunar Distance versus Date",
    subtitle = "(with full and new moons overlaid)"
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(shape = c(NA, 16, 4)))
  ) +
  ggplot2::scale_color_manual(
    values = rep("#000000", 3), name = "Legend"
  ) +
  ggplot2::scale_x_datetime(labels = scales::date_format("%d-%b-%Y")) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::theme_bw() +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(color = "#000000"),
    axis.text.y = ggplot2::element_text(color = "#000000"),
    legend.background = ggplot2::element_rect(color = NA, fill = "transparent"),
    legend.box.background = ggplot2::element_rect(
      color = NA, fill = "transparent"
    ),
    panel.background = ggplot2::element_rect(fill = "transparent"),
    plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    text = ggplot2::element_text(color = "#000000", family = "Times New Roman")
  ) +
  ggplot2::xlab("\ndate (UTC)") +
  ggplot2::ylab("distance (km)\n")

# Display graph.
p
```

![](supermoon_files/figure-gfm/example-graph-paper-style-1.png)<!-- -->

``` r
# Save graph, small format for research paper use.
ggplot2::ggsave(
  file = "supermoon_graph.eps", device = cairo_ps, family = "Times", scale = 0.7
)
```

    ## Saving 4.9 x 3.5 in image

``` r
# Example 3
# ---------
#
# Get all 21st century supermoons.
# Define "supermoon" as a full moon and lunar perigee separated by <= 12 h.

cat(
  "Example 3\n",
  "---------\n\n",
  "Get all 21st century supermoons.\n",
  "Define 'supermoon' as a full moon and lunar perigee separated by <= 12 h.\n",
  sep = ""
)

g_lunar_apsis <- "perigee"
g_lunar_phase <- "full moon"

# Get the 21st century full moons.
list_full_moons <- get_list_estimates(
  get_lunar_phase_jde,
  as.POSIXct("01-01-2001 00:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"),
  as.POSIXct("31-12-2100 00:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"),
  mpfr("29.4", 128L)
)

# Get the 21st century lunar perigees.
list_lunar_perigees <- get_list_estimates(
  get_lunar_apsis_jde,
  as.POSIXct("01-01-2001 00:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"),
  as.POSIXct("31-12-2100 00:00:00", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC"),
  mpfr("27.4", 128L)
)

# Build the supermoons list.
# Each list element will hold:
#   (date-time of full moon, date-time of perigee, distance).
# N.B. A 12 h UTC or DT time difference = 0.5 JDE.
list_supermoons <- list()
num_full_moons <- length(list_full_moons)
num_lunar_perigees <- length(list_lunar_perigees)
match_found <- FALSE
#
for (x in 1:num_full_moons) {
  for (y in 1:num_lunar_perigees) {

    # Only assess, in detail, relatively close data points.
    if (abs(list_full_moons[[x]]$value - list_lunar_perigees[[y]]$value) <= 2) {

      if (list_full_moons[[x]]$value < list_lunar_perigees[[y]]$value) {
        if (
          (
            list_lunar_perigees[[y]]$value -
              list_lunar_perigees[[y]]$uncertainty
          ) -
            (
              list_full_moons[[x]]$value +
                list_full_moons[[x]]$uncertainty
            ) <= 0.5
        ) {
          match_found <- TRUE
          break
        }
      } else {
        if (
          (
            list_full_moons[[x]]$value -
              list_full_moons[[x]]$uncertainty
          ) -
            (
              list_lunar_perigees[[y]]$value +
                list_lunar_perigees[[y]]$uncertainty
            ) <= 0.5
        ) {
          match_found <- TRUE
          break
        }
      }

    }
  }

  if (match_found) {
    # a) Convert time of full moon from JDE to UTC.
    # b) Convert time of perigee from JDE to UTC.
    # c) Get lunar distance at perigee.

    phase_date_dynamical <- calculate_with_uncertainty(
      time_jde_to_dynamical, list_full_moons[[x]]$value,
      list_full_moons[[x]]$uncertainty, FALSE
    )

    apsis_date_dynamical <- calculate_with_uncertainty(
      time_jde_to_dynamical, list_lunar_perigees[[y]]$value,
      list_lunar_perigees[[y]]$uncertainty, FALSE
    )

    list_supermoons[[length(list_supermoons) + 1]] <- list(
      "full.moon" = calculate_with_uncertainty(
        time_dynamical_to_universal, phase_date_dynamical$value,
        phase_date_dynamical$uncertainty, TRUE
      ),
      "perigee" = calculate_with_uncertainty(
        time_dynamical_to_universal, apsis_date_dynamical$value,
        apsis_date_dynamical$uncertainty, TRUE
      ),
      "distance" = calculate_with_uncertainty(
        get_lunar_distance, list_lunar_perigees[[y]]$value,
        list_lunar_perigees[[y]]$uncertainty, TRUE
      )
    )

    # Delete the perigee from the perigee list.
    list_lunar_perigees[- y]
    num_lunar_perigees <- length(list_lunar_perigees)

    match_found <- FALSE
  }
}


# Display results.
cat("\n---------------------------------------------\n")
#
for (x in seq_along(list_supermoons)) {

  cat("supermoon #", x, "\n", sep = "")

  label <- paste0(" ", g_lunar_phase, ":")
  if (list_supermoons[[x]]$full.moon$uncertainty < 60) {
    significant_figures <- 2
  } else {
    significant_figures <- 3
  }
  cat(
    label, replicate(g_max_indent - nchar(label), " "),
    format(
      lubridate::round_date(
        list_supermoons[[x]]$full.moon$value, unit = "seconds"
      ),
      "%d-%b-%Y %H:%M:%S"
    ),
    " +- ",
    format(
      hms::as_hms(
        signif(
          as.numeric(
            list_supermoons[[x]]$full.moon$uncertainty, units = "secs"
          ),
          significant_figures
        )
      ),
      "%H:%M:%S"
    ),
    " UTC\n",
    sep = ""
  )

  label <- paste0(" ", g_lunar_apsis, ":")
  cat(
    label, replicate(g_max_indent - nchar(label), " "),
    format(
      lubridate::round_date(
        list_supermoons[[x]]$perigee$value, unit = "seconds"
      ),
      "%d-%b-%Y %H:%M:%S"
    ),
    " +- ",
    format(
      hms::as_hms(
        signif(
          as.numeric(list_supermoons[[x]]$perigee$uncertainty, units = "secs"),
          4
        )
      ),
      "%H:%M:%S"
    ),
    " UTC\n",
    sep = ""
  )

  label <- paste0(" distance at ", g_lunar_apsis, ":")
  if (list_supermoons[[x]]$distance$uncertainty < 1) {
    significant_figures <- 7
  } else {
    significant_figures <- 6
  }
  cat(
    label, replicate(g_max_indent - nchar(label), " "),
    signif(asNumeric(list_supermoons[[x]]$distance$value), significant_figures),
    " +- ",
    signif(asNumeric(list_supermoons[[x]]$distance$uncertainty), 1), " km\n\n",
    sep = ""
  )
}

cat("---------------------------------------------")
```

    ## Example 3
    ## ---------
    ## 
    ## Get all 21st century supermoons.
    ## Define 'supermoon' as a full moon and lunar perigee separated by <= 12 h.
    ## 
    ## ---------------------------------------------
    ## supermoon #1
    ##  full moon:                08-Feb-2001 07:11:33 +- 00:00:17 UTC
    ##  perigee:                  07-Feb-2001 22:18:46 +- 00:31:00 UTC
    ##  distance at perigee:      356854 +- 2 km
    ## 
    ## supermoon #2
    ##  full moon:                27-Feb-2002 09:16:38 +- 00:00:17 UTC
    ##  perigee:                  27-Feb-2002 19:47:12 +- 00:31:00 UTC
    ##  distance at perigee:      356903 +- 5 km
    ## 
    ## supermoon #3
    ##  full moon:                28-Mar-2002 18:24:50 +- 00:00:17 UTC
    ##  perigee:                  28-Mar-2002 07:41:41 +- 00:31:00 UTC
    ##  distance at perigee:      357014 +- 4 km
    ## 
    ## supermoon #4
    ##  full moon:                16-Apr-2003 19:35:34 +- 00:00:17 UTC
    ##  perigee:                  17-Apr-2003 04:57:37 +- 00:31:00 UTC
    ##  distance at perigee:      357156.6 +- 0.5 km
    ## 
    ## supermoon #5
    ##  full moon:                16-May-2003 03:35:56 +- 00:00:17 UTC
    ##  perigee:                  15-May-2003 15:39:16 +- 00:31:00 UTC
    ##  distance at perigee:      357451 +- 1 km
    ## 
    ## supermoon #6
    ##  full moon:                03-Jun-2004 04:19:33 +- 00:00:17 UTC
    ##  perigee:                  03-Jun-2004 13:09:47 +- 00:31:00 UTC
    ##  distance at perigee:      357250 +- 3 km
    ## 
    ## supermoon #7
    ##  full moon:                02-Jul-2004 11:08:51 +- 00:00:17 UTC
    ##  perigee:                  01-Jul-2004 22:59:39 +- 00:31:00 UTC
    ##  distance at perigee:      357448.8 +- 0.8 km
    ## 
    ## supermoon #8
    ##  full moon:                21-Jul-2005 11:00:08 +- 00:00:17 UTC
    ##  perigee:                  21-Jul-2005 19:44:29 +- 00:31:00 UTC
    ##  distance at perigee:      357161 +- 2 km
    ## 
    ## supermoon #9
    ##  full moon:                19-Aug-2005 17:52:54 +- 00:00:17 UTC
    ##  perigee:                  19-Aug-2005 05:32:04 +- 00:31:00 UTC
    ##  distance at perigee:      357399 +- 6 km
    ## 
    ## supermoon #10
    ##  full moon:                07-Sep-2006 18:41:56 +- 00:00:17 UTC
    ##  perigee:                  08-Sep-2006 03:07:14 +- 00:31:00 UTC
    ##  distance at perigee:      357178 +- 2 km
    ## 
    ## supermoon #11
    ##  full moon:                26-Oct-2007 04:51:34 +- 00:00:17 UTC
    ##  perigee:                  26-Oct-2007 11:51:14 +- 00:31:00 UTC
    ##  distance at perigee:      356755 +- 3 km
    ## 
    ## supermoon #12
    ##  full moon:                12-Dec-2008 16:37:12 +- 00:00:17 UTC
    ##  perigee:                  12-Dec-2008 21:36:59 +- 00:31:00 UTC
    ##  distance at perigee:      356569 +- 4 km
    ## 
    ## supermoon #13
    ##  full moon:                30-Jan-2010 06:17:31 +- 00:00:17 UTC
    ##  perigee:                  30-Jan-2010 09:02:39 +- 00:31:00 UTC
    ##  distance at perigee:      356594 +- 1 km
    ## 
    ## supermoon #14
    ##  full moon:                19-Mar-2011 18:10:02 +- 00:00:17 UTC
    ##  perigee:                  19-Mar-2011 19:09:04 +- 00:31:00 UTC
    ##  distance at perigee:      356579 +- 4 km
    ## 
    ## supermoon #15
    ##  full moon:                06-May-2012 03:35:07 +- 00:00:17 UTC
    ##  perigee:                  06-May-2012 03:32:45 +- 00:31:00 UTC
    ##  distance at perigee:      356954 +- 1 km
    ## 
    ## supermoon #16
    ##  full moon:                23-Jun-2013 11:32:13 +- 00:00:17 UTC
    ##  perigee:                  23-Jun-2013 11:09:24 +- 00:31:00 UTC
    ##  distance at perigee:      356991.5 +- 0.3 km
    ## 
    ## supermoon #17
    ##  full moon:                10-Aug-2014 18:09:17 +- 00:00:17 UTC
    ##  perigee:                  10-Aug-2014 17:43:01 +- 00:31:00 UTC
    ##  distance at perigee:      356898 +- 3 km
    ## 
    ## supermoon #18
    ##  full moon:                28-Sep-2015 02:50:29 +- 00:00:17 UTC
    ##  perigee:                  28-Sep-2015 01:45:55 +- 00:31:00 UTC
    ##  distance at perigee:      356876.7 +- 0.6 km
    ## 
    ## supermoon #19
    ##  full moon:                14-Nov-2016 13:52:08 +- 00:00:17 UTC
    ##  perigee:                  14-Nov-2016 11:22:36 +- 00:31:00 UTC
    ##  distance at perigee:      356512 +- 3 km
    ## 
    ## supermoon #20
    ##  full moon:                02-Jan-2018 02:24:04 +- 00:00:17 UTC
    ##  perigee:                  01-Jan-2018 21:54:30 +- 00:31:00 UTC
    ##  distance at perigee:      356568 +- 3 km
    ## 
    ## supermoon #21
    ##  full moon:                19-Feb-2019 15:53:29 +- 00:00:17 UTC
    ##  perigee:                  19-Feb-2019 09:05:45 +- 00:31:00 UTC
    ##  distance at perigee:      356763 +- 2 km
    ## 
    ## supermoon #22
    ##  full moon:                08-Apr-2020 02:34:57 +- 00:00:17 UTC
    ##  perigee:                  07-Apr-2020 18:08:21 +- 00:31:00 UTC
    ##  distance at perigee:      356911 +- 4 km
    ## 
    ## supermoon #23
    ##  full moon:                27-Apr-2021 03:31:30 +- 00:00:17 UTC
    ##  perigee:                  27-Apr-2021 15:24:14 +- 00:31:00 UTC
    ##  distance at perigee:      357377.5 +- 1 km
    ## 
    ## supermoon #24
    ##  full moon:                26-May-2021 11:13:56 +- 00:00:17 UTC
    ##  perigee:                  26-May-2021 01:51:56 +- 00:31:00 UTC
    ##  distance at perigee:      357312 +- 1 km
    ## 
    ## supermoon #25
    ##  full moon:                14-Jun-2022 11:51:39 +- 00:00:17 UTC
    ##  perigee:                  14-Jun-2022 23:21:25 +- 00:31:00 UTC
    ##  distance at perigee:      357435 +- 2 km
    ## 
    ## supermoon #26
    ##  full moon:                13-Jul-2022 18:37:31 +- 00:00:17 UTC
    ##  perigee:                  13-Jul-2022 09:08:03 +- 00:31:00 UTC
    ##  distance at perigee:      357263.4 +- 0.8 km
    ## 
    ## supermoon #27
    ##  full moon:                01-Aug-2023 18:31:32 +- 00:00:17 UTC
    ##  perigee:                  02-Aug-2023 05:51:50 +- 00:31:00 UTC
    ##  distance at perigee:      357312 +- 1 km
    ## 
    ## supermoon #28
    ##  full moon:                31-Aug-2023 01:35:32 +- 00:00:17 UTC
    ##  perigee:                  30-Aug-2023 15:51:07 +- 00:31:00 UTC
    ##  distance at perigee:      357185 +- 4 km
    ## 
    ## supermoon #29
    ##  full moon:                18-Sep-2024 02:34:24 +- 00:00:17 UTC
    ##  perigee:                  18-Sep-2024 13:26:30 +- 00:31:00 UTC
    ##  distance at perigee:      357286 +- 0.5 km
    ## 
    ## supermoon #30
    ##  full moon:                17-Oct-2024 11:26:29 +- 00:00:17 UTC
    ##  perigee:                  17-Oct-2024 00:45:36 +- 00:31:00 UTC
    ##  distance at perigee:      357173 +- 2 km
    ## 
    ## supermoon #31
    ##  full moon:                05-Nov-2025 13:19:20 +- 00:00:17 UTC
    ##  perigee:                  05-Nov-2025 22:28:58 +- 00:31:00 UTC
    ##  distance at perigee:      356833.2 +- 0.4 km
    ## 
    ## supermoon #32
    ##  full moon:                04-Dec-2025 23:14:07 +- 00:00:17 UTC
    ##  perigee:                  04-Dec-2025 11:06:01 +- 00:31:00 UTC
    ##  distance at perigee:      356965 +- 3 km
    ## 
    ## supermoon #33
    ##  full moon:                24-Dec-2026 01:28:17 +- 00:00:47 UTC
    ##  perigee:                  24-Dec-2026 08:29:32 +- 00:31:38 UTC
    ##  distance at perigee:      356652 +- 1 km
    ## 
    ## supermoon #34
    ##  full moon:                10-Feb-2028 15:03:35 +- 00:00:56 UTC
    ##  perigee:                  10-Feb-2028 19:53:16 +- 00:31:00 UTC
    ##  distance at perigee:      356679.7 +- 0.6 km
    ## 
    ## supermoon #35
    ##  full moon:                30-Mar-2029 02:26:28 +- 00:00:47 UTC
    ##  perigee:                  30-Mar-2029 05:39:48 +- 00:31:39 UTC
    ##  distance at perigee:      356666 +- 2 km
    ## 
    ## supermoon #36
    ##  full moon:                17-May-2030 11:19:03 +- 00:00:56 UTC
    ##  perigee:                  17-May-2030 13:45:20 +- 00:31:39 UTC
    ##  distance at perigee:      357017 +- 2 km
    ## 
    ## supermoon #37
    ##  full moon:                04-Jul-2031 19:01:11 +- 00:00:57 UTC
    ##  perigee:                  04-Jul-2031 21:13:33 +- 00:31:39 UTC
    ##  distance at perigee:      357009.2 +- 0.9 km
    ## 
    ## supermoon #38
    ##  full moon:                21-Aug-2032 01:46:37 +- 00:00:57 UTC
    ##  perigee:                  21-Aug-2032 03:51:38 +- 00:31:40 UTC
    ##  distance at perigee:      356881 +- 2 km
    ## 
    ## supermoon #39
    ##  full moon:                08-Oct-2033 10:58:01 +- 00:00:57 UTC
    ##  perigee:                  08-Oct-2033 12:11:20 +- 00:31:40 UTC
    ##  distance at perigee:      356824 +- 2 km
    ## 
    ## supermoon #40
    ##  full moon:                25-Nov-2034 22:32:00 +- 00:00:58 UTC
    ##  perigee:                  25-Nov-2034 22:06:15 +- 00:31:40 UTC
    ##  distance at perigee:      356448 +- 3 km
    ## 
    ## supermoon #41
    ##  full moon:                13-Jan-2036 11:15:59 +- 00:00:58 UTC
    ##  perigee:                  13-Jan-2036 08:47:17 +- 00:31:41 UTC
    ##  distance at perigee:      356519 +- 2 km
    ## 
    ## supermoon #42
    ##  full moon:                02-Mar-2037 00:27:56 +- 00:00:59 UTC
    ##  perigee:                  01-Mar-2037 19:47:57 +- 00:31:00 UTC
    ##  distance at perigee:      356711 +- 3 km
    ## 
    ## supermoon #43
    ##  full moon:                19-Apr-2038 10:35:40 +- 00:00:50 UTC
    ##  perigee:                  19-Apr-2038 04:30:20 +- 00:31:42 UTC
    ##  distance at perigee:      356843 +- 6 km
    ## 
    ## supermoon #44
    ##  full moon:                06-Jun-2039 18:47:35 +- 00:01:00 UTC
    ##  perigee:                  06-Jun-2039 12:01:15 +- 00:31:00 UTC
    ##  distance at perigee:      357207 +- 2 km
    ## 
    ## supermoon #45
    ##  full moon:                24-Jul-2040 02:05:23 +- 00:01:00 UTC
    ##  perigee:                  23-Jul-2040 19:14:58 +- 00:31:21 UTC
    ##  distance at perigee:      357112 +- 2 km
    ## 
    ## supermoon #46
    ##  full moon:                10-Sep-2041 09:23:36 +- 00:00:17 UTC
    ##  perigee:                  10-Sep-2041 02:12:15 +- 00:31:22 UTC
    ##  distance at perigee:      357006 +- 6 km
    ## 
    ## supermoon #47
    ##  full moon:                28-Oct-2042 19:48:32 +- 00:00:43 UTC
    ##  perigee:                  28-Oct-2042 11:27:14 +- 00:31:43 UTC
    ##  distance at perigee:      356973 +- 2 km
    ## 
    ## supermoon #48
    ##  full moon:                16-Nov-2043 21:52:28 +- 00:00:53 UTC
    ##  perigee:                  17-Nov-2043 09:11:04 +- 00:31:44 UTC
    ##  distance at perigee:      356948 +- 1 km
    ## 
    ## supermoon #49
    ##  full moon:                16-Dec-2043 08:01:48 +- 00:01:01 UTC
    ##  perigee:                  15-Dec-2043 22:00:39 +- 00:31:44 UTC
    ##  distance at perigee:      356772 +- 4 km
    ## 
    ## supermoon #50
    ##  full moon:                03-Jan-2045 10:20:13 +- 00:01:02 UTC
    ##  perigee:                  03-Jan-2045 19:24:57 +- 00:31:22 UTC
    ##  distance at perigee:      356774 +- 1 km
    ## 
    ## supermoon #51
    ##  full moon:                01-Feb-2045 21:05:30 +- 00:00:44 UTC
    ##  perigee:                  01-Feb-2045 08:42:54 +- 00:31:44 UTC
    ##  distance at perigee:      357105 +- 1 km
    ## 
    ## supermoon #52
    ##  full moon:                20-Feb-2046 23:43:59 +- 00:01:02 UTC
    ##  perigee:                  21-Feb-2046 06:42:54 +- 00:31:22 UTC
    ##  distance at perigee:      356806 +- 2 km
    ## 
    ## supermoon #53
    ##  full moon:                10-Apr-2047 10:35:12 +- 00:00:54 UTC
    ##  perigee:                  10-Apr-2047 16:08:00 +- 00:31:00 UTC
    ##  distance at perigee:      356790 +- 4 km
    ## 
    ## supermoon #54
    ##  full moon:                27-May-2048 18:57:03 +- 00:01:03 UTC
    ##  perigee:                  27-May-2048 23:56:00 +- 00:31:46 UTC
    ##  distance at perigee:      357115 +- 1 km
    ## 
    ## supermoon #55
    ##  full moon:                15-Jul-2049 02:29:15 +- 00:01:04 UTC
    ##  perigee:                  15-Jul-2049 07:18:19 +- 00:31:23 UTC
    ##  distance at perigee:      357062 +- 1 km
    ## 
    ## supermoon #56
    ##  full moon:                01-Sep-2050 09:30:34 +- 00:00:56 UTC
    ##  perigee:                  01-Sep-2050 14:02:31 +- 00:31:24 UTC
    ##  distance at perigee:      356899 +- 4 km
    ## 
    ## supermoon #57
    ##  full moon:                19-Oct-2051 19:12:50 +- 00:01:06 UTC
    ##  perigee:                  19-Oct-2051 22:40:51 +- 00:31:48 UTC
    ##  distance at perigee:      356809 +- 2 km
    ## 
    ## supermoon #58
    ##  full moon:                06-Dec-2052 07:17:40 +- 00:01:07 UTC
    ##  perigee:                  06-Dec-2052 08:52:07 +- 00:31:49 UTC
    ##  distance at perigee:      356425 +- 4 km
    ## 
    ## supermoon #59
    ##  full moon:                23-Jan-2054 20:07:38 +- 00:01:08 UTC
    ##  perigee:                  23-Jan-2054 19:37:41 +- 00:31:51 UTC
    ##  distance at perigee:      356512 +- 1 km
    ## 
    ## supermoon #60
    ##  full moon:                13-Mar-2055 08:56:47 +- 00:01:09 UTC
    ##  perigee:                  13-Mar-2055 06:24:41 +- 00:31:26 UTC
    ##  distance at perigee:      356698 +- 4 km
    ## 
    ## supermoon #61
    ##  full moon:                29-Apr-2056 18:30:49 +- 00:01:10 UTC
    ##  perigee:                  29-Apr-2056 14:48:16 +- 00:31:27 UTC
    ##  distance at perigee:      356811 +- 5 km
    ## 
    ## supermoon #62
    ##  full moon:                17-Jun-2057 02:18:10 +- 00:01:03 UTC
    ##  perigee:                  16-Jun-2057 22:07:57 +- 00:31:54 UTC
    ##  distance at perigee:      357136 +- 2 km
    ## 
    ## supermoon #63
    ##  full moon:                04-Aug-2058 09:37:28 +- 00:01:13 UTC
    ##  perigee:                  04-Aug-2058 05:22:10 +- 00:31:55 UTC
    ##  distance at perigee:      356996 +- 2 km
    ## 
    ## supermoon #64
    ##  full moon:                21-Sep-2059 17:18:17 +- 00:01:14 UTC
    ##  perigee:                  21-Sep-2059 12:34:20 +- 00:31:57 UTC
    ##  distance at perigee:      356863 +- 5 km
    ## 
    ## supermoon #65
    ##  full moon:                08-Nov-2060 04:17:14 +- 00:01:15 UTC
    ##  perigee:                  07-Nov-2060 22:11:18 +- 00:31:58 UTC
    ##  distance at perigee:      356812 +- 3 km
    ## 
    ## supermoon #66
    ##  full moon:                26-Dec-2061 16:52:34 +- 00:01:16 UTC
    ##  perigee:                  26-Dec-2061 08:55:32 +- 00:31:59 UTC
    ##  distance at perigee:      356619 +- 3 km
    ## 
    ## supermoon #67
    ##  full moon:                14-Jan-2063 19:11:30 +- 00:01:17 UTC
    ##  perigee:                  15-Jan-2063 06:21:14 +- 00:31:30 UTC
    ##  distance at perigee:      356937 +- 2 km
    ## 
    ## supermoon #68
    ##  full moon:                13-Feb-2063 05:48:26 +- 00:01:18 UTC
    ##  perigee:                  12-Feb-2063 19:32:29 +- 00:32:00 UTC
    ##  distance at perigee:      356964 +- 3 km
    ## 
    ## supermoon #69
    ##  full moon:                03-Mar-2064 08:18:35 +- 00:01:19 UTC
    ##  perigee:                  03-Mar-2064 17:31:58 +- 00:31:31 UTC
    ##  distance at perigee:      356971.3 +- 0.4 km
    ## 
    ## supermoon #70
    ##  full moon:                01-Apr-2064 17:40:07 +- 00:01:10 UTC
    ##  perigee:                  01-Apr-2064 05:28:34 +- 00:32:01 UTC
    ##  distance at perigee:      357235.9 +- 0.7 km
    ## 
    ## supermoon #71
    ##  full moon:                20-Apr-2065 18:35:47 +- 00:01:20 UTC
    ##  perigee:                  21-Apr-2065 02:33:38 +- 00:32:02 UTC
    ##  distance at perigee:      356952 +- 2 km
    ## 
    ## supermoon #72
    ##  full moon:                08-Jun-2066 02:30:39 +- 00:01:21 UTC
    ##  perigee:                  08-Jun-2066 10:05:32 +- 00:32:04 UTC
    ##  distance at perigee:      357248 +- 2 km
    ## 
    ## supermoon #73
    ##  full moon:                26-Jul-2067 09:58:15 +- 00:01:22 UTC
    ##  perigee:                  26-Jul-2067 17:23:09 +- 00:32:05 UTC
    ##  distance at perigee:      357148.8 +- 0.4 km
    ## 
    ## supermoon #74
    ##  full moon:                11-Sep-2068 17:18:56 +- 00:01:23 UTC
    ##  perigee:                  12-Sep-2068 00:16:36 +- 00:31:33 UTC
    ##  distance at perigee:      356953 +- 4 km
    ## 
    ## supermoon #75
    ##  full moon:                30-Oct-2069 03:35:07 +- 00:01:25 UTC
    ##  perigee:                  30-Oct-2069 09:15:28 +- 00:31:34 UTC
    ##  distance at perigee:      356831 +- 1 km
    ## 
    ## supermoon #76
    ##  full moon:                17-Dec-2070 16:05:22 +- 00:01:26 UTC
    ##  perigee:                  17-Dec-2070 19:40:29 +- 00:32:09 UTC
    ##  distance at perigee:      356443 +- 4 km
    ## 
    ## supermoon #77
    ##  full moon:                04-Feb-2072 04:55:16 +- 00:01:27 UTC
    ##  perigee:                  04-Feb-2072 06:25:22 +- 00:31:35 UTC
    ##  distance at perigee:      356545.9 +- 0.6 km
    ## 
    ## supermoon #78
    ##  full moon:                23-Mar-2073 17:16:51 +- 00:01:28 UTC
    ##  perigee:                  23-Mar-2073 16:56:50 +- 00:31:35 UTC
    ##  distance at perigee:      356722 +- 3 km
    ## 
    ## supermoon #79
    ##  full moon:                11-May-2074 02:17:28 +- 00:01:30 UTC
    ##  perigee:                  11-May-2074 01:01:44 +- 00:31:36 UTC
    ##  distance at perigee:      356815 +- 4 km
    ## 
    ## supermoon #80
    ##  full moon:                28-Jun-2075 09:46:14 +- 00:01:31 UTC
    ##  perigee:                  28-Jun-2075 08:12:49 +- 00:32:13 UTC
    ##  distance at perigee:      357100 +- 1 km
    ## 
    ## supermoon #81
    ##  full moon:                14-Aug-2076 17:11:27 +- 00:01:32 UTC
    ##  perigee:                  14-Aug-2076 15:30:05 +- 00:31:37 UTC
    ##  distance at perigee:      356914 +- 3 km
    ## 
    ## supermoon #82
    ##  full moon:                02-Oct-2077 01:20:21 +- 00:01:33 UTC
    ##  perigee:                  01-Oct-2077 22:58:31 +- 00:32:16 UTC
    ##  distance at perigee:      356757 +- 7 km
    ## 
    ## supermoon #83
    ##  full moon:                19-Nov-2078 12:52:08 +- 00:01:26 UTC
    ##  perigee:                  19-Nov-2078 08:56:51 +- 00:32:17 UTC
    ##  distance at perigee:      356690 +- 1 km
    ## 
    ## supermoon #84
    ##  full moon:                07-Jan-2080 01:44:33 +- 00:01:36 UTC
    ##  perigee:                  06-Jan-2080 19:49:31 +- 00:32:18 UTC
    ##  distance at perigee:      356508 +- 4 km
    ## 
    ## supermoon #85
    ##  full moon:                23-Feb-2081 14:26:48 +- 00:01:28 UTC
    ##  perigee:                  23-Feb-2081 06:18:34 +- 00:31:40 UTC
    ##  distance at perigee:      356861 +- 2 km
    ## 
    ## supermoon #86
    ##  full moon:                14-Mar-2082 16:45:03 +- 00:01:30 UTC
    ##  perigee:                  15-Mar-2082 04:16:53 +- 00:32:21 UTC
    ##  distance at perigee:      357176 +- 1 km
    ## 
    ## supermoon #87
    ##  full moon:                13-Apr-2082 01:45:08 +- 00:01:38 UTC
    ##  perigee:                  12-Apr-2082 15:53:20 +- 00:32:21 UTC
    ##  distance at perigee:      357106 +- 2 km
    ## 
    ## supermoon #88
    ##  full moon:                02-May-2083 02:29:20 +- 00:01:31 UTC
    ##  perigee:                  02-May-2083 12:56:04 +- 00:31:41 UTC
    ##  distance at perigee:      357150 +- 3 km
    ## 
    ## supermoon #89
    ##  full moon:                31-May-2083 09:41:53 +- 00:01:40 UTC
    ##  perigee:                  30-May-2083 23:06:39 +- 00:31:41 UTC
    ##  distance at perigee:      357247 +- 1 km
    ## 
    ## supermoon #90
    ##  full moon:                18-Jun-2084 10:00:15 +- 00:01:32 UTC
    ##  perigee:                  18-Jun-2084 20:15:01 +- 00:32:23 UTC
    ##  distance at perigee:      357415 +- 2 km
    ## 
    ## supermoon #91
    ##  full moon:                17-Jul-2084 17:01:15 +- 00:01:40 UTC
    ##  perigee:                  17-Jul-2084 06:12:54 +- 00:31:42 UTC
    ##  distance at perigee:      357471.7 +- 0.7 km
    ## 
    ## supermoon #92
    ##  full moon:                05-Aug-2085 17:28:51 +- 00:01:33 UTC
    ##  perigee:                  06-Aug-2085 03:30:21 +- 00:32:25 UTC
    ##  distance at perigee:      357270 +- 2 km
    ## 
    ## supermoon #93
    ##  full moon:                04-Sep-2085 00:40:57 +- 00:01:40 UTC
    ##  perigee:                  03-Sep-2085 13:42:42 +- 00:32:25 UTC
    ##  distance at perigee:      357232 +- 2 km
    ## 
    ## supermoon #94
    ##  full moon:                23-Sep-2086 01:14:37 +- 00:01:40 UTC
    ##  perigee:                  23-Sep-2086 10:35:14 +- 00:32:26 UTC
    ##  distance at perigee:      357041 +- 5 km
    ## 
    ## supermoon #95
    ##  full moon:                22-Oct-2086 09:55:31 +- 00:01:35 UTC
    ##  perigee:                  21-Oct-2086 21:59:24 +- 00:31:43 UTC
    ##  distance at perigee:      357178 +- 5 km
    ## 
    ## supermoon #96
    ##  full moon:                10-Nov-2087 12:04:45 +- 00:01:40 UTC
    ##  perigee:                  10-Nov-2087 19:53:59 +- 00:32:27 UTC
    ##  distance at perigee:      356889.7 +- 0.5 km
    ## 
    ## supermoon #97
    ##  full moon:                28-Dec-2088 00:57:03 +- 00:01:50 UTC
    ##  perigee:                  28-Dec-2088 06:31:50 +- 00:32:29 UTC
    ##  distance at perigee:      356502 +- 6 km
    ## 
    ## supermoon #98
    ##  full moon:                14-Feb-2090 13:39:02 +- 00:01:50 UTC
    ##  perigee:                  14-Feb-2090 17:12:55 +- 00:32:30 UTC
    ##  distance at perigee:      356620.4 +- 0.6 km
    ## 
    ## supermoon #99
    ##  full moon:                04-Apr-2091 01:31:02 +- 00:01:50 UTC
    ##  perigee:                  04-Apr-2091 03:25:59 +- 00:32:31 UTC
    ##  distance at perigee:      356784 +- 4 km
    ## 
    ## supermoon #100
    ##  full moon:                21-May-2092 09:59:48 +- 00:01:50 UTC
    ##  perigee:                  21-May-2092 11:10:29 +- 00:31:46 UTC
    ##  distance at perigee:      356854 +- 3 km
    ## 
    ## supermoon #101
    ##  full moon:                08-Jul-2093 17:13:31 +- 00:01:50 UTC
    ##  perigee:                  08-Jul-2093 18:16:58 +- 00:32:34 UTC
    ##  distance at perigee:      357098 +- 0.5 km
    ## 
    ## supermoon #102
    ##  full moon:                26-Aug-2094 00:51:25 +- 00:01:50 UTC
    ##  perigee:                  26-Aug-2094 01:39:18 +- 00:31:48 UTC
    ##  distance at perigee:      356867 +- 2 km
    ## 
    ## supermoon #103
    ##  full moon:                13-Oct-2095 09:30:14 +- 00:01:50 UTC
    ##  perigee:                  13-Oct-2095 09:25:10 +- 00:32:36 UTC
    ##  distance at perigee:      356687 +- 5 km
    ## 
    ## supermoon #104
    ##  full moon:                29-Nov-2096 21:33:48 +- 00:02:00 UTC
    ##  perigee:                  29-Nov-2096 19:43:27 +- 00:32:38 UTC
    ##  distance at perigee:      356609 +- 2 km
    ## 
    ## supermoon #105
    ##  full moon:                17-Jan-2098 10:35:40 +- 00:02:00 UTC
    ##  perigee:                  17-Jan-2098 06:41:26 +- 00:32:39 UTC
    ##  distance at perigee:      356437 +- 3 km
    ## 
    ## supermoon #106
    ##  full moon:                06-Mar-2099 22:59:22 +- 00:02:00 UTC
    ##  perigee:                  06-Mar-2099 16:59:32 +- 00:31:50 UTC
    ##  distance at perigee:      356797 +- 3 km
    ## 
    ## supermoon #107
    ##  full moon:                24-Apr-2100 09:43:26 +- 00:02:00 UTC
    ##  perigee:                  24-Apr-2100 02:13:30 +- 00:31:00 UTC
    ##  distance at perigee:      357012 +- 1 km
    ## 
    ## ---------------------------------------------
