21st Century Supermoon Estimation in R
================

``` r
knitr::opts_chunk$set(cache = TRUE, results = "hold") # Group output by chunk.

suppressMessages(library(extrafont))    # Extra fonts, used for graphs.
suppressMessages(library(ggplot2))      # Graph plotting library.
suppressMessages(library(lubridate))    # Date-time manipulation utilities.
suppressMessages(library(pkgcond))      # Suppress inefficiency warnings when converting mpfr values to numerics.
suppressMessages(library(Rmpfr))        # Arbitrary-precision arithmetic.
suppressMessages(library(scales))       # For formatting date axes in graphs.


AngleDegreesToRadians <- function(angle.degrees) {

   # Convert angles in degrees to radians.
   #
   # Input:
   #    angle.degrees   (an angle in degrees).
   #
   # N.B. Short-circuit order optimised for 21st century supermoon processing.

   if ((angle.degrees >= G.360.mpfr) || (angle.degrees < G.zero.mpfr)) {
      # N.B. The %% operator has higher precedence than *.
      return(angle.degrees %% G.360.mpfr * G.pi.on.180.mpfr)
   } else {
      return(angle.degrees * G.pi.on.180.mpfr)
   }
}


AngleRadiansToDegrees <- function(angle.radians) {

   # Convert angles in radians to degrees.
   #
   # Input:
   #    angle.radians   (an angle in radians).

   two.pi <- Const("pi", prec = 128) * 2

   # Convert angles to (0 <= theta < 2*Pi) radians, if necessary, by taking the angle modulo 2*Pi radians.
   if ((angle.radians < G.zero.mpfr) || (angle.radians >= two.pi)) {
      angle.radians <- angle.radians %% two.pi
   }

   return(angle.radians / Const("pi", prec = 128) * 180)
}


CalculateWithUncertainty <- function(function.name, value, uncertainty, receive.uncertainty) {

   # Call another function, propagating uncertainty.
   #    N.B. It assumes the function called also returns a value with an associated uncertainty.
   #    N.B. If time-dates are being manipulated, pass the uncertainty as a lubridate duration.
   #
   # Input:
   #    function.name         the function to invoke
   #    value                 the value to pass to the function
   #    uncertainty           the value's uncertainty to pass to the function
   #    receive.uncertainty   a boolean indicating if the invoked function itself returns a value with an uncertainty

   result.low <- function.name(value - uncertainty)
   result.mid <- function.name(value)
   result.high <- function.name(value + uncertainty)

   if (receive.uncertainty) {
      result.vector <- c(result.low$value - result.low$uncertainty, result.low$value + result.low$uncertainty,
                         result.mid$value - result.mid$uncertainty, result.mid$value + result.mid$uncertainty,
                         result.high$value - result.high$uncertainty, result.high$value + result.high$uncertainty)
   } else {
      result.vector <- c(result.low, result.mid, result.high)
   }

   result.min <- min(result.vector)
   result.max <- max(result.vector)

   return(list("value" = mean(c(result.min, result.max)), "uncertainty" = {result.max - result.min} * 0.5))
}


GetDelT <- function(date.target) {
   
   # Calculate the difference delT, in seconds, between dynamical time and universal time (UT1).
   #    See Five Millennium Canon of Solar Eclipses: -1999 to +3000 (2000 BCE to 3000 CE).
   #
   # Input:
   #    date.target   (date should be universal time (UT1).
   #                  (It's assumed the target date will be in the 21st century.)
   #
   # N.B. An exact decimal year is used here.
   #      The reference material uses a less accurate value rounded to the middle of months.
   # N.B. Where relevant, a lunar secular acceleration of -26 arcsec/cy^2 is assumed.

   # Convert input date to decimal years. E.g. 1987.25 is the end of March, 1987.
   date.target <- G.zero.mpfr + decimal_date(date.target)

   if (date.target < G.threshold.GetDelT) {
      # As known delT values exist until 31-Dec-2020 23:59:59, detT shouldn't be estimated on or prior to that.
      # Allow a 900 millisecond buffer for processing the UT1-UTC uncertainty (always <= 0.9 s).
      # So, G.threshold.GetDelT is (01-01-2021 00:00:00) - 0.9 s.

      cat("\nError in function GetDelT...")
      cat("\nInvalid date processing; date.target:\n")
      print(date.target, digits=20)
      stopifnot(FALSE)

   } else if (date.target < mpfr("2050", 128)) {
      t <- date.target - 2000
      delT <- t * {mpfr("0.005589", 128) * t + mpfr("0.32217", 128)} + mpfr("62.92", 128)

   } else {
      # delT <- -20.0 - 0.5628 * (2150.0 - date.target) + 32.0 * (((date.target - 1820.0) * 0.01) ^ 2)
      # Or more simply...
      delT <- date.target * {mpfr("0.0032", 128) * date.target - mpfr("11.0852", 128)} + mpfr("9369.66", 128)
   }

   # N.B. Here uncertainty is in seconds; algorithm's intrinsic uncertainty.
   N <- date.target - 2005
   return(list("value" = delT,
               "uncertainty" = 365.25 * N * sqrt(N * G.uncertainty.constant.GetDelT * {1 + N * mpfr("0.0004", 128)}) * mpfr("0.001", 128)))
}


GetEccentricAnomaly <- function(time.target) {

   # Get the eccentric anomaly, for a given time (in Julian centuries from Epoch J2000.0).
   #
   # Input:
   #   time.target   Target time, in Julian centuries from Epoch J2000.0.

   return(1 - time.target * {mpfr("0.002516", 128) + time.target * mpfr("0.0000074", 128)})
}


GetListEstimates <- function(function.name, date.start, date.end, interval) {

   # Get a list of dates (in JDE format) for nominated apsides or phases over a selected date range.
   #
   #   date.start   (start of the date range, in Coordinated Universal Time).
   #   date.end     (end of the date range, in Coordinated Universal Time).
   #   interval     (decimal specifying how often, in days, to check for the nearest apside or phase, e.g. 29.4.)
   #                (N.B. For apsides, ensure it's < the anomalistic month over the date range, else apsides will be missed.)
   #                (N.B. For phases, ensure it's < the synodic month over the date range, else phases will be missed.)
   #                (N.B. An mpfr value is recommended.)

   # Total list elements needed is ceiling({date.end - date.start} / ddays(interval)).
   interval.days <- ddays(interval)
   list.estimates <- lapply(0:{ceiling({date.end - date.start} / interval.days) - 1},
         function(x) {
            target.dynamical <- TimeUniversalToDynamical(date.start + x * interval.days)
            CalculateWithUncertainty(function.name, target.dynamical$value, target.dynamical$uncertainty, TRUE)
         }
   )

   # Remove duplicates.
   if (G.verbosity) {
      cat("\nList length, initial: ", length(list.estimates), sep="")
      list.estimates <- unique(list.estimates)
      cat("\nList length, final: ", length(list.estimates), sep="")
   } else {
      list.estimates <- unique(list.estimates)
   }

   return(list.estimates)
}


GetLunarApsisJDE <- function(date.target) {

   # Calculate lunar apsis dates.
   #    See Astronomical Algorithms, Meeus, chapter 50.
   #
   # Input:
   #    date.target   (date should be a Dynamical Time, as a date/time value).

   # Determine new moon coefficient "k".
   #    k = 0 corresponds to the perigee of 22nd Dec, 1999.
   #    Negative values of k give lunar apsides before that perigee.
   #    - For perigee, round to the nearest integer.
   #    - For apogee, round to the nearest integer + 0.5.
   #    Assumes target apsis is a global variable.
   #
   #    N.B. Ensure x.5 values are rounded towards zero if negative, but away from zero if positive.
   #    N.B. Ensure odd and even numbers are rounded consistently (R doesn't do this by default).
   #    N.B. Ensure target dynamical time is converted to decimal years before use.
   #
   k <- {decimal_date(date.target) - mpfr("1999.97", 128)} * mpfr("13.2555", 128)
   #
   if (G.lunar.apsis == "perigee") {
      k <- floor(0.5 + k)

   } else if (G.lunar.apsis == "apogee") {
      k <- 0.5 + floor(k)

   } else {
      cat("\nError in function GetLunarApsisJDE...")
      cat("\nInvalid lunar apsis; G.lunar.apsis:\n")
      print(G.lunar.apsis)
      stopifnot(FALSE)
   }

   # Determine the approximate time "T" in Julian centuries from epoch 2000.
   T <- k / mpfr("1325.55", 128)
   T.x.T <- T * T

   # Determine the mean perigee or apogee, in Julian Ephemeris Days "JDE".
   # Specify the calculation using Horner's Method.
   mean.apsis <- T.x.T * {mpfr("-0.0006691", 128) + T * {T * mpfr("0.0000000052", 128) - mpfr("0.000001098", 128)}} + mpfr("27.55454989", 128) * k + mpfr("2451534.6698", 128)

   # Calculate the Moon's mean elongation, in radians.
   D <- AngleDegreesToRadians(T.x.T * {T * {T * mpfr("0.000000055", 128) - mpfr("0.00001156", 128)} - mpfr("0.0100383", 128)} + mpfr("171.9179", 128) + mpfr("335.9106046", 128) * k)

   # Calculate the Sun's mean anomaly, in radians.
   M <- AngleDegreesToRadians(T.x.T * {-T * mpfr("0.0000010", 128) - mpfr("0.0008130", 128)} + mpfr("347.3477", 128) + mpfr("27.1577721", 128) * k)

   # Calculate the Moon's argument of latitude, in radians.
   F <- AngleDegreesToRadians(T.x.T * {-T * mpfr("0.0000148", 128) - mpfr("0.0125053", 128)} + mpfr("316.6109", 128) + mpfr("364.5287911", 128) * k)
   F.x.2 <- 2 * F

   # Calculate the first set of periodic corrections, in days.
   if (G.lunar.apsis == "perigee") {
      # Set apsis-specific uncertainty from Astronomical Algorithms, Meeus, chapter 50.
      # N.B. The worst-case uncertainty, not mean uncertainty, is taken here.
      # N.B. Values here are minutes converted into JDE days, using 1 day = 1440 minutes.
      apsis.uncertainty <- 31 / mpfr("1440", 128)

      # Calculate periodic corrections for perigee.

      # Perigee timing periodic calculation coefficients.
      vec.sin <-  c(        2 * D,                 4 * D,             6 * D,               8 * D,
                        2 * D - M,                     M,            10 * D,           4 * D - M,
                        6 * D - M,                12 * D,                 D,           8 * D - M,
                           14 * D,                 F.x.2,             3 * D,          10 * D - M,
                           16 * D,            12 * D - M,             5 * D,       2 * D + F.x.2,
                           18 * D,            14 * D - M,             7 * D,           2 * D + M,
                           20 * D,                 D + M,        16 * D - M,           4 * D + M,
                            9 * D,         4 * D + F.x.2,       2 * {D - M},       4 * D - 2 * M,
                    6 * D - 2 * M,                22 * D,        18 * D - M,           6 * D + M,
                           11 * D,             8 * D + M,     4 * D - F.x.2,       6 * D + F.x.2,
                        3 * D + M,             5 * D + M,            13 * D,          20 * D - M,
                    3 * D + 2 * M, 4 * D + F.x.2 - 2 * M,         D + 2 * M,          22 * D - M,
                            4 * F,         6 * D - F.x.2, 2 * D - F.x.2 + M,               2 * M,
                        F.x.2 - M,       2 * {D + F.x.2},     F.x.2 - 2 * M, 2 * {D + M} - F.x.2,
                           24 * D,           4 * {D - F},       2 * {D + M},               D - M)

      corrections.periodic <- sum({G.vec.perigee.coeff + G.vec.perigee.subtrahend * T} * sin(vec.sin))

   } else {
      # Set apsis-specific uncertainty.
      # N.B. Values here are minutes converted into JDE days, using 1 day = 1440 minutes.
      apsis.uncertainty <- 3 / mpfr("1440", 128)

      # Calculate periodic corrections for apogee.
      vec.sin <- c(        2 * D,         4 * D,                 M,     2 * D - M,
                           F.x.2,             D,             6 * D,     4 * D - M,
                   2 * D + F.x.2,         D + M,             8 * D,     6 * D - M,
                   2 * D - F.x.2,   2 * {D - M},             3 * D, 4 * D + F.x.2,
                       8 * D - M, 4 * D - 2 * M,            10 * D,     3 * D + M,
                           2 * M,     2 * D + M,       2 * {D + M}, 6 * D + F.x.2,
                   6 * D - 2 * M,    10 * D - M,             5 * D, 4 * D - F.x.2,
                       F.x.2 + M,        12 * D, 2 * D + F.x.2 - M,         D - M)

      corrections.periodic <- sum({G.vec.apogee.coeff + G.vec.apogee.subtrahend * T} * sin(vec.sin))
   }

   # Diagnostic information.
   if (G.verbosity) {
      cat("   GetLunarApsisJDE settings...\n", sep="")
      cat("\n   k:")
      print(k, digits = 20)
      cat("\n   T:")
      print(T, digits = 20)
      cat("\n   D (degrees):")
      print(AngleRadiansToDegrees(D), digits = 20)
      cat("\n   M (degrees):")
      print(AngleRadiansToDegrees(M), digits = 20)
      cat("\n   F (degrees):")
      print(AngleRadiansToDegrees(F), digits = 20)
      cat("\n   mean.apsis: ")
      print(mean.apsis, digits = 20)
      cat("\n   corrections.periodic: ")
      print(corrections.periodic, digits = 20)
   }

   return(list("value" = mean.apsis + corrections.periodic, "uncertainty" = apsis.uncertainty))
}


GetLunarDistance <- function(JDE.mpfr) {

   # Get the lunar distance, in kilometres, for a given JDE (Julian Ephemeris Day).
   #    See Astronomical Algorithms, Meeus, chapter 47.
   #
   # Input:
   #   JDE.mpfr   Target date, in Julian Ephemeris Days, best passed in mpfr format with, e.g. 128-bit accuracy.
   #              E.g. R stores a passed value of 2443259.9 as 2443259.8999999999 by default.
   #                   This causes errors of +- 1 second easily.
   
   # Input must be >= 0.
   if (JDE.mpfr < 0) {
      cat("\nError in function GetLunarDistance...")
      cat("\nInput dynamical time; JDE.mpfr:\n")
      print(JDE.mpfr, digits=20)
      cat("\nNeed dynamical time >= 0.\n")
      stopifnot(FALSE)
   }

   # Express the target date as "T" in Julian centuries from epoch 2000.
   # N.B. See Astronomical Algorithms, Meeus, chapter 22, p. 143.
   T <- {JDE.mpfr - 2451545} / 36525

   # Calculate E, to correct for the (changing) eccentricity of the Earth's orbit around the Sun.
   E <- GetEccentricAnomaly(T)
   E.x.E <- E * E

   # Calculate the Moon's mean elongation, in radians.
   D <- AngleDegreesToRadians(mpfr("297.8501921", 128) + T * {mpfr("445267.1114034", 128) - T * {mpfr("0.0018819", 128) - T * {G.lunar.mean.elongation.constant - T / 113065000}}})

   # Calculate the Sun's mean anomaly, in radians.
   M <- AngleDegreesToRadians(mpfr("357.5291092", 128) + T * {mpfr("35999.0502909", 128) - T * {mpfr("0.0001536", 128) - T / 24490000}})

   # Calculate the Moon's mean anomaly, in radians.
   M.prime <- AngleDegreesToRadians(mpfr("134.9633964", 128) + T * {mpfr("477198.8675055", 128) + T * {mpfr("0.0087414", 128) + T * {G.lunar.mean.anomaly.constant - T / 14712000}}})

   # Calculate the Moon's argument of latitude, in radians.
   F.x.2 <- 2 * AngleDegreesToRadians(mpfr("93.272095", 128) + T * {mpfr("483202.0175233", 128) - T * {mpfr("0.0036539", 128) + T * {G.lunar.argument.latitude.constant - T / 863310000}}})

   # Calculate periodic distance terms.
   vec.E <-     c( G.1.mpfr,                G.1.mpfr,                G.1.mpfr,                G.1.mpfr,
                          E,                G.1.mpfr,                G.1.mpfr,                       E,
                   G.1.mpfr,                       E,                       E,                G.1.mpfr,
                          E,                G.1.mpfr,                G.1.mpfr,                G.1.mpfr,
                   G.1.mpfr,                G.1.mpfr,                       E,                       E,
                   G.1.mpfr,                       E,                       E,                G.1.mpfr,
                   G.1.mpfr,                G.1.mpfr,                       E,                       E,
                   G.1.mpfr,                   E.x.E,                       E,                   E.x.E,
                   G.1.mpfr,                       E,                G.1.mpfr,                       E,
                          E,                   E.x.E,                   E.x.E,                G.1.mpfr,
                   G.1.mpfr,                       E,                G.1.mpfr,                G.1.mpfr,               
                      E.x.E,                G.1.mpfr)

   distance.value <- mpfr("385000.56", 128) +
                     sum(G.vec.ld.coeff * vec.E *
                        cos(D * G.mask.distance.D +
                            F.x.2 * G.mask.distance.F.x.2 +
                            M * G.mask.distance.M +
                            M.prime * G.mask.distance.M.prime)) *
                     mpfr("0.001", 128)

   # Diagnostic information.
   if (G.verbosity) {
      cat("   GetLunarDistance...\n", sep="")
      cat("\n   T:")
      print(T, digits = 20)
      cat("\n   E:")
      print(E, digits = 20)
      cat("\n   D (degrees):")
      print(AngleRadiansToDegrees(D), digits = 20)
      cat("\n   M (degrees):")
      print(AngleRadiansToDegrees(M), digits = 20)
      cat("\n   M.prime (degrees):")
      print(AngleRadiansToDegrees(M.prime), digits = 20)
      cat("\n   F (degrees):")
      print(AngleRadiansToDegrees(F.x.2 * 0.5), digits = 20)
      cat("\n   distance.periodic: ")
      print({distance.value - mpfr("385000.56", 128)} * mpfr("1000", 128), digits = 20)
   }

   # Determine uncertainty by comparison to the complete ELP-2000/82 theory.
   # Load a Fortran shared library for ELP-2000/82.
   # N.B. Ideally, this shouldn't be loaded/unloaded each time it's needed, but implementation is unclear.
   dyn.load(G.library.path)
   stopifnot(is.loaded("elp82b_2"))
   result <- .Fortran("ELP82B_2", as.double(JDE.mpfr), 0, as.integer(1), double(3), as.integer(-1))
   dyn.unload(G.library.path)

   # ELP-2000/82 distance is the length derived from its three Cartesian coordinate vectors.
   return(list("value" = distance.value,
               "uncertainty" = abs(sqrt(result[[4]][1] * result[[4]][1] +
                                        result[[4]][2] * result[[4]][2] +
                                        result[[4]][3] * result[[4]][3]) - distance.value)))
}


GetLunarPhaseJDE <- function(date.target) {

   # Calculate lunar phase dates.
   #    See Astronomical Algorithms, Meeus, chapter 49.
   #
   # Input:
   #    date.target   (date should be a Dynamical Time, as a date/time value).

   # Determine new moon coefficient "k".
   #    k = 0 corresponds to the New Moon of 6th Jan, 2000.
   #    Negative values of k give lunar phases before the year 2000.
   #    - If new moon, round to the nearest integer.
   #    - If first quarter moon, round to the nearest integer + 0.25.
   #    - If full moon, round to the nearest integer + 0.5.
   #    - If last quarter moon, round to the nearest integer + 0.75.
   #    Assumes target phase is a global variable.
   #
   #    N.B. Ensure x.5 values are rounded towards zero if negative, but away from zero if positive.
   #    N.B. Ensure odd and even numbers are rounded consistently (R doesn't do this by default).
   #    N.B. Ensure target dynamical time is converted to decimal years before use.
   #
   k <- {decimal_date(date.target) - 2000} * mpfr("12.3685", 128)
   #
   if (G.lunar.phase == "full moon") {
      k <- mpfr("0.5", 128) + floor(k)
      
      # Set phase-specific uncertainty from Astronomical Algorithms, Meeus, chapter 49, p. 354.
      # N.B. The worst-case uncertainty, not mean uncertainty, is taken here.
      # N.B. Values here are seconds converted into JDE days, using 1 day = 86,400 seconds.
      phase.uncertainty <- mpfr("17.4", 128) / 86400

   } else if (G.lunar.phase == "last quarter") {
      k <- sign(k) * mpfr("0.75", 128) + floor(mpfr("0.5", 128) - sign(k) * mpfr("0.75", 128) + k)
      phase.uncertainty <- mpfr("13", 128) / 86400

   } else if (G.lunar.phase == "new moon") {
      k <- floor(mpfr("0.5", 128) + k)
      phase.uncertainty <- mpfr("16.4", 128) / 86400

   } else if (G.lunar.phase == "first quarter") {
      k <- sign(k) * mpfr("0.25", 128) + floor(mpfr("0.5", 128) - sign(k) * mpfr("0.25", 128) + k)
      phase.uncertainty <- mpfr("15.3", 128) / 86400

   } else {
      cat("\nError in function GetLunarPhaseJDE...")
      cat("\nInvalid lunar phase; G.lunar.phase:\n")
      print(G.lunar.phase)
      stopifnot(FALSE)
   }

   # Determine the approximate time "T" in Julian centuries from epoch 2000.
   T <- k / mpfr("1236.85", 128)
   T.x.T <- T * T

   # Determine the mean lunar phase, in Julian Ephemeris Days "JDE".
   # Specify the calculation using Horner's Method.
   mean.phase <- T.x.T * {0.00015437 + T * {T * 0.00000000073 - 0.00000015}} + 29.530588861 * k + 2451550.09766

   # Calculate E, to correct for the (changing) eccentricity of the Earth's orbit around the Sun.
   E <- GetEccentricAnomaly(T)
   E.x.E <- E * E

   # Calculate the Sun's mean anomaly, in radians.
   M <- AngleDegreesToRadians(T.x.T * {-T * 0.00000011 - 0.0000014} + 2.5534 + 29.1053567 * k)

   # Calculate the Moon's mean anomaly, in radians.
   M.prime <- AngleDegreesToRadians(T.x.T * {-T.x.T * 0.000000058 + T * 0.00001238 + 0.0107582} + 201.5643 + 385.81693528 * k)

   # Calculate the Moon's argument of latitude, in radians.
   F.x.2 <- 2 * AngleDegreesToRadians(T.x.T * {T * {T * 0.000000011 - 0.00000227} - 0.0016118} + 160.7108 + 390.67050284 * k)

   # Calculate the longitude of ascending node of the lunar orbit, in radians.
   omega <- AngleDegreesToRadians(T.x.T * {T * 0.00000215 + 0.0020672} + 124.7746 - 1.56375588 * k)

   # Calculate the first set of periodic corrections, in days.
   if (G.lunar.phase == "full moon" || G.lunar.phase == "new moon") {

      vec.ph.E <- c(G.1.mpfr,
                    E,
                    G.1.mpfr,
                    G.1.mpfr,
                    E,
                    E,
                    E.x.E,
                    G.1.mpfr,
                    G.1.mpfr,
                    E,
                    G.1.mpfr,
                    E,
                    E,
                    E,
                    G.1.mpfr,
                    G.1.mpfr,
                    G.1.mpfr,
                    G.1.mpfr,
                    G.1.mpfr,
                    G.1.mpfr,
                    G.1.mpfr,
                    G.1.mpfr,
                    G.1.mpfr,
                    G.1.mpfr,
                    G.1.mpfr)

      vec.ph.qrt.sin <- c(M.prime,
                          M,
                          2 * M.prime,
                          F.x.2,
                          M.prime - M,
                          M.prime + M,
                          2 * M,
                          M.prime - F.x.2,
                          M.prime + F.x.2,
                          2 * M.prime + M,
                          3 * M.prime,
                          M + F.x.2,
                          M - F.x.2,
                          2 * M.prime - M,
                          omega,
                          M.prime + 2 * M,
                          2 * M.prime - F.x.2,
                          3 * M,
                          M.prime + M - F.x.2,
                          2 * M.prime + F.x.2,
                          M.prime + M + F.x.2,
                          M.prime - M + F.x.2,
                          M.prime - M - F.x.2,
                          3 * M.prime + M,
                          4 * M.prime)

      if (G.lunar.phase == "full moon") {
         # A full moon is being processed.
         corrections.periodic <- sum(G.vec.ph.full.coeff * vec.ph.E * sin(vec.ph.qrt.sin))
      } else {
         # Assume a new moon is being processed.
         corrections.periodic <- sum(G.vec.ph.new.coeff * vec.ph.E * sin(vec.ph.qrt.sin))
      }
      
   } else {
      # Assume a quarter phase is being processed.

      vec.ph.qrt.E <- c(G.1.mpfr,
                        E,
                        E,
                        G.1.mpfr,
                        G.1.mpfr,
                        E,
                        E.x.E,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr,
                        E,
                        E,
                        E,
                        E.x.E,
                        E,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr)

      vec.ph.qrt.sin <- c(M.prime,
                          M,
                          M.prime + M,
                          2 * M.prime,
                          F.x.2,
                          M.prime - M,
                          2 * M,
                          M.prime - F.x.2,
                          M.prime + F.x.2,
                          3 * M.prime,
                          2 * M.prime - M,
                          M + F.x.2,
                          M - F.x.2,
                          M.prime + 2 * M,
                          2 * M.prime + M,
                          omega,
                          M.prime - M - F.x.2,
                          2 * M.prime + F.x.2,
                          M.prime + M + F.x.2,
                          M.prime - 2 * M,
                          M.prime + M - F.x.2,
                          3 * M,
                          2 * M.prime - F.x.2,
                          M.prime - M + F.x.2,
                          3 * M.prime + M)

      corrections.periodic <- sum(G.vec.ph.qrt.coeff * vec.ph.qrt.E * sin(vec.ph.qrt.sin))

      # Supplemental quarter phase corrections.
      vec.ph.sup.E <- c(E,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr,
                        G.1.mpfr)

      vec.ph.sup.cos <- c(M,
                          M.prime,
                          M.prime - M,
                          M.prime + M,
                          F.x.2)

      W <- 0.00306 + sum(G.vec.ph.sup.coeff * vec.ph.sup.E  * cos(vec.ph.sup.cos))

      if (G.lunar.phase == "first quarter") {
         # Apply supplemental first quarter corrections.
         corrections.periodic <- corrections.periodic + W

      } else {
         # Apply supplemental last quarter corrections.
         corrections.periodic <- corrections.periodic - W
      }

      # Diagnostic information.
      if (G.verbosity) {
         cat("   GetLunarPhaseJDE settings (quarter-specific)...\n", sep="")
         cat("\n   W:")
         print(W, digits = 20)
      }
   }

   # Calculate additional corrections, based on "planetary arguments" in days.
   corrections.additional <- sum(G.vec.ph.add.coeff.outer * sin(AngleDegreesToRadians(G.vec.ph.add.addend + G.vec.ph.add.coeff.inner * k - G.vec.ph.add.subtrahend * T.x.T)))

   # Diagnostic information.
   if (G.verbosity) {
      cat("   GetLunarPhaseJDE settings...\n", sep="")
      cat("\n   k:")
      print(k, digits = 20)
      cat("\n   T:")
      print(T, digits = 20)
      cat("\n   E:")
      print(E, digits = 20)
      cat("\n   mean.phase: ")
      print(mean.phase, digits = 20)
      cat("\n   corrections.periodic: ")
      print(corrections.periodic, digits = 20)
      cat("\n   corrections.additional: ")
      print(corrections.additional, digits = 20)
   }

   return(list("value" = corrections.additional + corrections.periodic + mean.phase,
               "uncertainty" = phase.uncertainty))
}


GetUTC <- function(DT) {

   # Estimate UTC from an input Dynamical Time.
   # Effectively reverses the operation in GetDelT (except it returns UTC, not UT1).
   #
   # Input:
   #    date.target   (date should be a Dynamical Time).
   #                  (It's assumed the recovered UTC date will be in the 21st century.)

   # Ensure a valid date range has been entered.
   # Known delT values exist until 31-Dec-2020 23:59:59.
   # Thus, UTC should be determined exactly on or prior to that.
   # N.B. Assume 01-Jan-2021 00:00:00.000 UTC = 01-01-2021 00:01:12.150 DT.
   #
   if (DT < dmy_hms("01-01-2021 00:01:12.150", tz="UTC")) {
      cat("\nError in function GetUTC...")
      cat("\nInput date:  ", format(DT, "%d-%b-%Y %H:%M:%OS"), " (assumed Dynamical Time)")
      cat("\nNeed date >= 01-01-2021 00:01:12.150 Dynamical Time", sep="")
      stopifnot(FALSE)
   }

   # Determine how many seconds were in the UTC year, get a inverse scaling factor.
   # N.B. Remember, delT is a function of the UTC year (not DT year).
   if (IsUTCLeapYear(DT)) {
      # 31622400 seconds = 366 days.
      scaling <- G.366.day.scaling
   } else {
      # 31536000 seconds = 365 days.
      scaling <- G.365.day.scaling
   }

   # Estimate value.
   # Remember, 0 = UTC (decimal year) + (ay^2 + by + c) (seconds) - DT (decimal year).
   DT.decimal <- G.zero.mpfr + decimal_date(DT)
   small.diff <- G.zero.mpfr + decimal_date(DT - G.5.minutes)

   if (DT < dmy_hms("01-01-2050 00:01:33.000", tz="UTC")) {
      # Assume 01-Jan-2050 00:00:00.000 UTC = 01-Jan-2050 00:01:33.000 DT.
      f1 <- function (x) {
               t <- x - 2000
               x + scaling * {t * {0.005589 * t + 0.32217} + 62.92} - DT.decimal
            }

   } else {
      f1 <- function (x) {
               x + scaling * {x * {0.0032 * x - 11.0852} + 9369.66} - DT.decimal
            }
   }
   
   # Return the UTC value as a date-time.
   #    N.B. By default unirootR warns if the root finding calculation fails to converge.
   #    N.B. Using tol = 1e-20 is generally fine, but we use the R numerical mantissa for added rigor.
   UTC <- unirootR(f1, interval = c(small.diff, DT.decimal), verbose = FALSE, tol = 1e-36)

   UTC.decimal <- date_decimal(asNumeric(UTC$root), tz = "UTC")

   return(list("value" = UTC.decimal,
               "uncertainty" = difftime(date_decimal(asNumeric(UTC$root + UTC$estim.prec * 0.5), tz = "UTC"),
                                        UTC.decimal, tz = "UTC", units = c("auto"))))
}


IsUTCLeapYear <- function(date.target) {

   # Test if a Dynamical Time is equivalent to a UTC leap year.
   # This information is used to conversion of DT to UTC date-times > 31-Dec-2020 23:59:59.999 (UTC).
   #
   # Input:
   #   date.target   (date should be a UTC time).
   #   N.B. An input date equivalent to a 21st century UTC date is assumed.

   # Convert input date to decimal format.
   date.target.decimal <- G.zero.mpfr + decimal_date(date.target)

   # Ensure only dates >= 01-Jan-2021 UTC are being processed.
   # N.B. Assume 01-Jan-2021 00:00:00.000 UTC = 01-Jan-2021 00:01:12.150 Dynamical Time.
   #                                           ~2021.000002287871666339920000000000000000.
   #      Technically, here, 01-Jan-2021 00:00:00.000 UTC = 01-Jan-2021 00:01:12.150 +- 00:00:03.261 Dynamical Time.
   #
   date.threshold.decimal <- mpfr("2021.000002287871666339920000000000000000", 128)
   if (date.target.decimal < date.threshold.decimal) {
      cat("\nWarning in function IsUTCLeapYear...")
      cat("\n   Input date: ", formatMpfr(date.target.decimal, NULL), sep="")
      cat("\n   Expected date >= ", formatMpfr(date.threshold.decimal, NULL), sep="")
   }

   # Get the year of the input Dynamical Time.
   year.target <- year(date.target)

   # Allow for Dynamical Times falling shortly after New Year's Eve of a UTC leap year.
   #
   #    E.g. 2024 is a leap year.
   #    Assume 01-Jan-2025 00:00:00.000 UTC --> 01-Jan-2025 00:01:14.467 +- 4.561s Dynamical Time.
   #                                           ~2025.000002361344968448980000000000000003.
   #    So 2025 Dynamical Time dates < 2025.000002361344968448980000000000000003 actually started in 2024 UTC.
   #

   # Get the index of the input year in the vector vec.post.leap.years, if it's there.
   # N.B. R indices start from 1, not 0.
   x <- match(year.target, G.seq.leap.year.check, nomatch = 0)
   if (x > 0) {
      if (date.target.decimal < G.vec.DT.transition[x]) {
         year.target <- year.target - 1
      }
   }

   return(leap_year(year.target))
}


IsValidDate <- function(date.target) {

   # This program is designed to process 21st century dates.
   #    I.e. 01-Jan-2001 to 31-Dec-2100.
   #
   # Input:
   #    date.target   (date should be a UTC time).

   if (date.target >= dmy_hms("01-01-2001 00:00:00", tz="UTC")) {
      if (date.target < dmy_hms("01-01-2101 00:00:00", tz="UTC")) {
         return(TRUE)
      } else {
         return(FALSE)
      }
   } else {
      return(FALSE)
   }
}


TimeDynamicalToJDE <- function(DT) {

   # Convert Dynamical Time to JDE (Julian Ephemeris Day).
   #    See Astronomical Algorithms, Meeus, chapter 7.
   #    N.B. Assumes a Gregorian calendar.
   #
   # Input:
   #   DT   (a date-time, in Dynamical Time format).

   Y <- G.zero.mpfr + year(DT)
   M <- G.zero.mpfr + month(DT)
   D <- G.zero.mpfr + day(DT)

   # Add the date's fractional part to D.
   D <- D + {DT - dmy_hms(paste(sprintf("%02d", as.integer(D)), "-", sprintf("%02d", as.integer(M)), "-", as.integer(Y), " 00:00:00", sep=""), tz="UTC")} / G.1.day

   if (M < 3) {
      Y <- Y - 1
      M <- M + 12
   }

   A <- floor(Y * 0.01)
   B <- 2 - A + floor(A * 0.25)

   return(floor(365.25 * {Y + 4716}) + floor(30.6001 * {M + 1}) + D + B - 1524.5)
}


TimeDynamicalToUniversal <- function(DT) {

   # Convert Dynamical Time to Coordinated Universal Time.
   #
   #    N.B. This function reverses the process in TimeUniversalToDynamical.
   #    N.B. It's assumed a 21st century date will be processed.
   #         I.e. 01-Jan-2001 to 31-Dec-2100.
   #         Limited support for other dates, for testing purposes, has been provided.
   #
   # Input:
   #   DT   (a date-time, in Dynamical Time format).
   # Output:
   #   UTC = DT - offset.

   # These values are given in UTC, to avoid accidental assignment to other time zones.
   # If the UTC label was omitted, they would be valid dynamical time values.
   #
   if (DT >= dmy_hms("01-01-2021 00:01:09.18", tz="UTC")) {
      # Get UTC via algorithmic estimation.
      #    N.B. GetUTC returns a date-time.
      #         GetUTC returns uncertainty as a time difference.
      return(GetUTC(DT))

   } else if (DT >= dmy_hms("01-01-2017 00:01:08.18", tz="UTC")) {
      # Offset = 42.18 + 27 seconds.
      return(list("value" = DT - G.69.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-2015 00:01:07.18", tz="UTC")) {
      # Offset = 42.18 + 26 seconds.
      return(list("value" = DT - G.68.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-2012 00:01:06.18", tz="UTC")) {
      # Offset = 42.18 + 25 seconds.
      return(list("value" = DT - G.67.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-2009 00:01:05.18", tz="UTC")) {
      # Offset = 42.18 + 24 seconds.
      return(list("value" = DT - G.66.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-2006 00:01:04.18", tz="UTC")) {
      # Offset = 42.18 + 23 seconds.
      return(list("value" = DT - G.65.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1999 00:01:03.18", tz="UTC")) {
      # Offset = 42.18 + 22 seconds.
      return(list("value" = DT - G.64.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-1997 00:01:02.18", tz="UTC")) {
      # Offset = 42.18 + 21 seconds.
      return(list("value" = DT - G.63.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1996 00:01:01.18", tz="UTC")) {
      # Offset = 42.18 + 20 seconds.
      return(list("value" = DT - G.62.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-1994 00:01:00.18", tz="UTC")) {
      # Offset = 42.18 + 19 seconds.
      return(list("value" = DT - G.61.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-1993 00:00:59.18", tz="UTC")) {
      # Offset = 42.18 + 18 seconds.
      return(list("value" = DT - G.60.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-1992 00:00:58.18", tz="UTC")) {
      # Offset = 42.18 + 17 seconds.
      return(list("value" = DT - G.59.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1991 00:00:57.18", tz="UTC")) {
      # Offset = 42.18 + 16 seconds.
      return(list("value" = DT - G.58.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1990 00:00:56.18", tz="UTC")) {
      # Offset = 42.18 + 15 seconds.
      return(list("value" = DT - G.57.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1988 00:00:55.18", tz="UTC")) {
      # Offset = 42.18 + 14 seconds.
      return(list("value" = DT - G.56.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-1985 00:00:54.18", tz="UTC")) {
      # Offset = 42.18 + 13 seconds.
      return(list("value" = DT - G.55.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-1983 00:00:53.18", tz="UTC")) {
      # Offset = 42.18 + 12 seconds.
      return(list("value" = DT - G.54.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-1982 00:00:52.18", tz="UTC")) {
      # Offset = 42.18 + 11 seconds.
      return(list("value" = DT - G.53.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-1981 00:00:51.18", tz="UTC")) {
      # Offset = 42.18 + 10 seconds.
      return(list("value" = DT - G.52.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1980 00:00:50.18", tz="UTC")) {
      # Offset = 42.18 + 9 seconds.
      return(list("value" = DT - G.51.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1979 00:00:49.18", tz="UTC")) {
      # Offset = 42.18 + 8 seconds.
      return(list("value" = DT - G.50.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1978 00:00:48.18", tz="UTC")) {
      # Offset = 42.18 + 7 seconds.
      return(list("value" = DT - G.49.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1977 00:00:47.18", tz="UTC")) {
      # Offset = 42.18 + 6 seconds.
      return(list("value" = DT - G.48.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1976 00:00:46.18", tz="UTC")) {
      # Offset = 42.18 + 5 seconds.
      return(list("value" = DT - G.47.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1975 00:00:45.18", tz="UTC")) {
      # Offset = 42.18 + 4 seconds.
      return(list("value" = DT - G.46.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1974 00:00:44.18", tz="UTC")) {
      # Offset = 42.18 + 3 seconds.
      return(list("value" = DT - G.45.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-01-1973 00:00:43.18", tz="UTC")) {
      # Offset = 42.18 + 2 seconds.
      return(list("value" = DT - G.44.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (DT >= dmy_hms("01-07-1972 00:00:42.18", tz="UTC")) {
      # Offset = 42.18 + 1 seconds.
      return(list("value" = DT - G.43.18.seconds, "uncertainty" = G.zero.seconds))

   } else {
      return(list("value" = DT - G.42.18.seconds, "uncertainty" = G.zero.seconds))
   }
}


TimeJDEToDynamical <- function(JDE.mpfr) {

   # Convert JDE (Julian Ephemeris Day) to Dynamical Time.
   #    See Astronomical Algorithms, Meeus, chapter 7.
   #    N.B. Error on p. 63: "2291 161" should be "2299 161".
   #
   # Input:
   #   JDE.mpfr   Target date, in Julian Ephemeris Days, best passed in mpfr format with, e.g. 128-bit accuracy.
   #              E.g. R stores a passed value of 2443259.9 as 2443259.8999999999 by default.
   #                   This causes errors of +- 1 second easily.

   # Input must be >= 0.
   if (JDE.mpfr < 0) {
      cat("\nError in function TimeJDEToDynamical...")
      cat("\nInput dynamical time; JDE.mpfr:\n")
      print(JDE.mpfr, digits=20)
      cat("\nNeed dynamical time >= 0.\n")
      stopifnot(FALSE)
   }

   JDE <- 0.5 + JDE.mpfr
   Z <- trunc(JDE)   # Integer part.
   F <- JDE - Z      # Fractional part.

   if (Z < 2299161) {
      A <- Z
   } else {
      alpha <- floor({Z - 1867216.25} / 36524.25)
      A <- Z + 1 + alpha - floor(alpha * 0.25)
   }

   B <- A + 1524
   C <- floor({B - 122.1} / 365.25)
   D <- floor(365.25 * C)
   E <- floor({B - D} / 30.6001)

   # Day, as a decimal.
   day <- B - D - floor(30.6001 * E) + F

   # Month, as an integer.
   if (E < 14) {
      month <- E - 1
   } else {
      month <- E - 13
   }

   # Year, as an integer.
   if (month > 2) {
      year <- C - 4716
   } else {
      year <- C - 4715
   }

   # The Dynamical Time returned here is stored as a UTC date-time format.
   # If the time zone label is neglected however, the date and time are the correct Dynamical Time.
   # N.B. Don't forget to add any fractional part of the day to the output.

   return(dmy_hms(paste(sprintf("%02d", as.integer(trunc(day))), "-", sprintf("%02d", as.integer(month)), "-", as.integer(year), " 00:00:00", sep=""), tz="UTC") + ddays(day - trunc(day)))
}


TimeUniversalToDynamical <- function(UTC) {

   # Convert Coordinated Universal Time to Dynamical Time.
   #    N.B. It's assumed a 21st century date will be processed.
   #         I.e. 01-Jan-2001 to 31-Dec-2100.
   #         Limited support for other dates, for testing purposes, has been provided.
   #
   # Input:
   #   UTC   (target date, in Coordinated Universal Time).
   #
   # N.B. The output value is stored as a UTC date-time format.
   #      If the time zone label is neglected however, the date and time are the correct Dynamical Time.
   #
   #    DT = UTC + 42.18 s + # leap seconds
   #
   # Ref: Polynomial approximations to Delta T, 1620-2000 AD (Meeus and Simons).
   #
   # N.B. Dynamical time can be determined exactly for dates from 1972 to the present, estimated after.
   # N.B. Number of leap seconds added from 1972-2020 = 27.
   # N.B. Timings are available from Wikipedia, the Astronomical Almanac, etc.

   if (UTC >= G.2021.start) {
      # Get delta T via algorithmic estimation.
      # N.B. The algorithm finds delta T for a UT (or UT1).
      # N.B. Thus, allow for uncertainty in UTC vs UT1 (always <= 0.9 s; as maintained by leap seconds).
      DelT <- CalculateWithUncertainty(GetDelT, UTC, G.900.milliseconds, TRUE)
      return(list("value" = UTC + dseconds(DelT$value), "uncertainty" = dseconds(DelT$uncertainty)))

   } else if (UTC >= G.2017.start) {
      # offset = 42.18 seconds + 27 leap seconds.
      # No forecast change to the end of 2021; keep leap seconds unchanged.
      return(list("value" = UTC + G.69.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.2015.mid) {
      # offset = 42.18 seconds + 26 leap seconds.
      return(list("value" = UTC + G.68.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.2012.mid) {
      # offset = 42.18 seconds + 25 leap seconds.
      return(list("value" = UTC + G.67.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.2009.start) {
      # offset = 42.18 seconds + 24 leap seconds.
      return(list("value" = UTC + G.66.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.2006.start) {
      # offset = 42.18 seconds + 23 leap seconds.
      return(list("value" = UTC + G.65.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1999.start) {
      # offset = 42.18 seconds + 22 leap seconds.
      return(list("value" = UTC + G.64.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1997.mid) {
      # offset = 42.18 seconds + 21 leap seconds.
      return(list("value" = UTC + G.63.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1996.start) {
      # offset = 42.18 seconds + 20 leap seconds.
      return(list("value" = UTC + G.62.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1994.mid) {
      # offset = 42.18 seconds + 19 leap seconds.
      return(list("value" = UTC + G.61.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1993.mid) {
      # offset = 42.18 seconds + 18 leap seconds.
      return(list("value" = UTC + G.60.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1992.mid) {
      # offset = 42.18 seconds + 17 leap seconds.
      return(list("value" = UTC + G.59.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1991.start) {
      # offset = 42.18 seconds + 16 leap seconds.
      return(list("value" = UTC + G.58.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1990.start) {
      # offset = 42.18 seconds + 15 leap seconds.
      return(list("value" = UTC + G.57.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1988.start) {
      # offset = 42.18 seconds + 14 leap seconds.
      return(list("value" = UTC + G.56.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1985.mid) {
      # offset = 42.18 seconds + 13 leap seconds.
      return(list("value" = UTC + G.55.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1983.mid) {
      # offset = 42.18 seconds + 12 leap seconds.
      return(list("value" = UTC + G.54.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1982.mid) {
      # offset = 42.18 seconds + 11 leap seconds.
      return(list("value" = UTC + G.53.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1981.mid) {
      # offset = 42.18 seconds + 10 leap seconds.
      return(list("value" = UTC + G.52.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1980.start) {
      # offset = 42.18 seconds + 9 leap seconds.
      return(list("value" = UTC + G.51.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1979.start) {
      # offset = 42.18 seconds + 8 leap seconds.
      return(list("value" = UTC + G.50.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1978.start) {
      # offset = 42.18 seconds + 7 leap seconds.
      return(list("value" = UTC + G.49.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1977.start) {
      # offset = 42.18 seconds + 6 leap seconds.
      return(list("value" = UTC + G.48.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1976.start) {
      # offset = 42.18 seconds + 5 leap seconds.
      return(list("value" = UTC + G.47.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1975.start) {
      # offset = 42.18 seconds + 4 leap seconds.
      return(list("value" = UTC + G.46.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1974.start) {
      # offset = 42.18 seconds + 3 leap seconds.
      return(list("value" = UTC + G.45.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1973.start) {
      # offset = 42.18 seconds + 2 leap seconds.
      return(list("value" = UTC + G.44.18.seconds, "uncertainty" = G.zero.seconds))

   } else if (UTC >= G.1972.mid) {
      # offset = 42.18 seconds + 1 leap seconds.
      return(list("value" = UTC + G.43.18.seconds, "uncertainty" = G.zero.seconds))

   } else {
      # offset = 42.18 seconds + 0 leap seconds.
      return(list("value" = UTC + G.42.18.seconds, "uncertainty" = G.zero.seconds))
   }
}


#######################
# Set global variables.

# Set output verbosity.
G.verbosity <- FALSE

# Set lunar apsis.
# G.lunar.apsis <- "apogee"
G.lunar.apsis <- "perigee"

# Set lunar phase.
# G.lunar.phase <- "new moon"
# G.lunar.phase <- "first quarter"
G.lunar.phase <- "full moon"
# G.lunar.phase <- "last quarter"

# Set the path to a Fortran shared library for ELP-2000/82.
# Designate the shared library extension in platform-independent manner.
G.library.path <- file.path(paste("elp82b_2", .Platform$dynlib.ext, sep=""))

# Control extra fonts installation.
G.extra.font.installed <- TRUE

# Set general purpose numeric globals.
G.max.indent <- 27
G.pi.on.180.mpfr <- Const("pi", prec = 128) / mpfr("180", 128)
G.360.mpfr <- mpfr("360", 128)
G.1.mpfr <- mpfr("1", 128)
G.zero.mpfr <- mpfr("0", 128)

# Set durations.
G.zero.seconds <- dseconds(G.zero.mpfr)
G.900.milliseconds <- dmilliseconds(mpfr("900", 128))
G.42.18.seconds <- dseconds(mpfr("42.18", 128))
G.43.18.seconds <- dseconds(mpfr("43.18", 128))
G.44.18.seconds <- dseconds(mpfr("44.18", 128))
G.45.18.seconds <- dseconds(mpfr("45.18", 128))
G.46.18.seconds <- dseconds(mpfr("46.18", 128))
G.47.18.seconds <- dseconds(mpfr("47.18", 128))
G.48.18.seconds <- dseconds(mpfr("48.18", 128))
G.49.18.seconds <- dseconds(mpfr("49.18", 128))
G.50.18.seconds <- dseconds(mpfr("50.18", 128))
G.51.18.seconds <- dseconds(mpfr("51.18", 128))
G.52.18.seconds <- dseconds(mpfr("52.18", 128))
G.53.18.seconds <- dseconds(mpfr("53.18", 128))
G.54.18.seconds <- dseconds(mpfr("54.18", 128))
G.55.18.seconds <- dseconds(mpfr("55.18", 128))
G.56.18.seconds <- dseconds(mpfr("56.18", 128))
G.57.18.seconds <- dseconds(mpfr("57.18", 128))
G.58.18.seconds <- dseconds(mpfr("58.18", 128))
G.59.18.seconds <- dseconds(mpfr("59.18", 128))
G.60.18.seconds <- dseconds(mpfr("60.18", 128))
G.61.18.seconds <- dseconds(mpfr("61.18", 128))
G.62.18.seconds <- dseconds(mpfr("62.18", 128))
G.63.18.seconds <- dseconds(mpfr("63.18", 128))
G.64.18.seconds <- dseconds(mpfr("64.18", 128))
G.65.18.seconds <- dseconds(mpfr("65.18", 128))
G.66.18.seconds <- dseconds(mpfr("66.18", 128))
G.67.18.seconds <- dseconds(mpfr("67.18", 128))
G.68.18.seconds <- dseconds(mpfr("68.18", 128))
G.69.18.seconds <- dseconds(mpfr("69.18", 128))
G.5.minutes <- dminutes(mpfr("5", 128))
G.1.day <- ddays(G.1.mpfr)

# Set time conversion constants.
G.2021.start <- dmy_hms("01-01-2021 00:00:00", tz="UTC")
G.2017.start <- dmy_hms("01-01-2017 00:00:00", tz="UTC")
G.2015.mid <- dmy_hms("01-07-2015 00:00:00", tz="UTC")
G.2012.mid <- dmy_hms("01-07-2012 00:00:00", tz="UTC")
G.2009.start <- dmy_hms("01-01-2009 00:00:00", tz="UTC")
G.2006.start <- dmy_hms("01-01-2006 00:00:00", tz="UTC")
G.1999.start <- dmy_hms("01-01-1999 00:00:00", tz="UTC")
G.1997.mid <- dmy_hms("01-07-1997 00:00:00", tz="UTC")
G.1996.start <- dmy_hms("01-01-1996 00:00:00", tz="UTC")
G.1994.mid <- dmy_hms("01-07-1994 00:00:00", tz="UTC")
G.1993.mid <- dmy_hms("01-07-1993 00:00:00", tz="UTC")
G.1992.mid <- dmy_hms("01-07-1992 00:00:00", tz="UTC")
G.1991.start <- dmy_hms("01-01-1991 00:00:00", tz="UTC")
G.1990.start <- dmy_hms("01-01-1990 00:00:00", tz="UTC")
G.1988.start <- dmy_hms("01-01-1988 00:00:00", tz="UTC")
G.1985.mid <- dmy_hms("01-07-1985 00:00:00", tz="UTC")
G.1983.mid <- dmy_hms("01-07-1983 00:00:00", tz="UTC")
G.1982.mid <- dmy_hms("01-07-1982 00:00:00", tz="UTC")
G.1981.mid <- dmy_hms("01-07-1981 00:00:00", tz="UTC")
G.1980.start <- dmy_hms("01-01-1980 00:00:00", tz="UTC")
G.1979.start <- dmy_hms("01-01-1979 00:00:00", tz="UTC")
G.1978.start <- dmy_hms("01-01-1978 00:00:00", tz="UTC")
G.1977.start <- dmy_hms("01-01-1977 00:00:00", tz="UTC")
G.1976.start <- dmy_hms("01-01-1976 00:00:00", tz="UTC")
G.1975.start <- dmy_hms("01-01-1975 00:00:00", tz="UTC")
G.1974.start <- dmy_hms("01-01-1974 00:00:00", tz="UTC")
G.1973.start <- dmy_hms("01-01-1973 00:00:00", tz="UTC")
G.1972.mid <- dmy_hms("01-07-1972 00:00:00", tz="UTC")

# Set thresholds and miscellaneous constants.
G.365.day.scaling <- 1 / mpfr("31536000", 128)
G.366.day.scaling <- 1 / mpfr("31622400", 128)
G.lunar.argument.latitude.constant <- 1 / mpfr("3526000", 128)
G.lunar.mean.anomaly.constant <- 1 / mpfr("69699", 128)
G.lunar.mean.elongation.constant <- 1 / mpfr("545868", 128)
G.threshold.GetDelT <- G.zero.mpfr + decimal_date(dmy_hms("01-01-2021 00:00:00", tz="UTC") - G.900.milliseconds)
G.uncertainty.constant.GetDelT <- mpfr("0.058", 128) / 3

# Set leap year check transition constants.
G.vec.DT.transition <- c(mpfr("2025.000002361344968448980000000000000003", 128),
                         mpfr("2029.000002440489652144610000000000000000", 128),
                         mpfr("2033.000002525305490053139999999999999999", 128),
                         mpfr("2037.000002615792482174580000000000000000", 128),
                         mpfr("2041.000002711950855882609999999999999999", 128),
                         mpfr("2045.000002813780383803530000000000000003", 128),
                         mpfr("2049.000002921281065937360000000000000003", 128),
                         mpfr("2053.000003143493358948040000000000000003", 128),
                         mpfr("2057.000003405644292797659999999999999996", 128),
                         mpfr("2061.000003671042577479970000000000000000", 128),
                         mpfr("2065.000003939687758247599999999999999994", 128),
                         mpfr("2069.000004211580289847920000000000000000", 128),
                         mpfr("2073.000004486719717533559999999999999995", 128),
                         mpfr("2077.000004765106496051890000000000000002", 128),
                         mpfr("2081.000005046740170655539999999999999999", 128),
                         mpfr("2085.000005331621196091869999999999999999", 128),
                         mpfr("2089.000005619748662866189999999999999998", 128),
                         mpfr("2093.000005911123935220530000000000000000", 128),
                         mpfr("2097.000006205745648912850000000000000004", 128))

# Set leap year check sequence.
G.seq.leap.year.check <- seq(2025, 2097, by=4)

# Lunar apogee calculation coefficients.
G.vec.apogee.coeff <- c(mpfr( "0.4392", 128),
                        mpfr( "0.0684", 128),
                        mpfr( "0.0456", 128),
                        mpfr( "0.0426", 128),
                        mpfr( "0.0212", 128),
                        mpfr("-0.0189", 128),
                        mpfr( "0.0144", 128),
                        mpfr( "0.0113", 128),
                        mpfr( "0.0047", 128),
                        mpfr( "0.0036", 128),
                        mpfr( "0.0035", 128),
                        mpfr( "0.0034", 128),
                        mpfr("-0.0034", 128),
                        mpfr( "0.0022", 128),
                        mpfr("-0.0017", 128),
                        mpfr( "0.0013", 128),
                        mpfr( "0.0011", 128),
                        mpfr( "0.0010", 128),
                        mpfr( "0.0009", 128),
                        mpfr( "0.0007", 128),
                        mpfr( "0.0006", 128),
                        mpfr( "0.0005", 128),
                        mpfr( "0.0005", 128),
                        mpfr( "0.0004", 128),
                        mpfr( "0.0004", 128),
                        mpfr( "0.0004", 128),
                        mpfr("-0.0004", 128),
                        mpfr("-0.0004", 128),
                        mpfr( "0.0003", 128),
                        mpfr( "0.0003", 128),
                        mpfr( "0.0003", 128),
                        mpfr("-0.0003", 128))
#
G.vec.apogee.subtrahend <- c(G.zero.mpfr,
                             G.zero.mpfr,
                             mpfr("-0.00011", 128),
                             mpfr("-0.00011", 128),
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr,
                             G.zero.mpfr)

# Lunar perigee calculation coefficients.
G.vec.perigee.coeff <- c(mpfr("-1.6769", 128),
                         mpfr( "0.4589", 128),
                         mpfr("-0.1856", 128),
                         mpfr( "0.0883", 128),
                         mpfr("-0.0773", 128),
                         mpfr( "0.0502", 128),
                         mpfr("-0.0460", 128),
                         mpfr( "0.0422", 128),
                         mpfr("-0.0256", 128),
                         mpfr( "0.0253", 128),
                         mpfr( "0.0237", 128),
                         mpfr( "0.0162", 128),
                         mpfr("-0.0145", 128),
                         mpfr( "0.0129", 128),
                         mpfr("-0.0112", 128),
                         mpfr("-0.0104", 128),
                         mpfr( "0.0086", 128),
                         mpfr( "0.0069", 128),
                         mpfr( "0.0066", 128),
                         mpfr("-0.0053", 128),
                         mpfr("-0.0052", 128),
                         mpfr("-0.0046", 128),
                         mpfr("-0.0041", 128),
                         mpfr( "0.0040", 128),
                         mpfr( "0.0032", 128),
                         mpfr("-0.0032", 128),
                         mpfr( "0.0031", 128),
                         mpfr("-0.0029", 128),
                         mpfr( "0.0027", 128),
                         mpfr( "0.0027", 128),
                         mpfr("-0.0027", 128),
                         mpfr( "0.0024", 128),
                         mpfr("-0.0021", 128),
                         mpfr("-0.0021", 128),
                         mpfr("-0.0021", 128),
                         mpfr( "0.0019", 128),
                         mpfr("-0.0018", 128),
                         mpfr("-0.0014", 128),
                         mpfr("-0.0014", 128),
                         mpfr("-0.0014", 128),
                         mpfr( "0.0014", 128),
                         mpfr("-0.0014", 128),
                         mpfr( "0.0013", 128),
                         mpfr( "0.0013", 128),
                         mpfr( "0.0011", 128),
                         mpfr("-0.0011", 128),
                         mpfr("-0.0010", 128),
                         mpfr("-0.0009", 128),
                         mpfr("-0.0008", 128),
                         mpfr( "0.0008", 128),
                         mpfr( "0.0008", 128),
                         mpfr( "0.0007", 128),
                         mpfr( "0.0007", 128),
                         mpfr( "0.0007", 128),
                         mpfr("-0.0006", 128),
                         mpfr("-0.0006", 128),
                         mpfr( "0.0006", 128),
                         mpfr( "0.0005", 128),
                         mpfr( "0.0005", 128),
                         mpfr("-0.0004", 128))
#
G.vec.perigee.subtrahend <- c(G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              mpfr( "0.00019", 128),
                              mpfr("-0.00013", 128),
                              G.zero.mpfr,
                              mpfr("-0.00011", 128),
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr)

# Lunar distance periodic calculation coefficients.
G.vec.ld.coeff <- c(-20905355, -3699111, -2955968, -569925,
                        48888,    -3149,   246158, -152138,
                      -170733,  -204586,  -129620,  108743,
                       104755,    10321,    79661,  -34782,
                       -23210,   -21636,    24208,   30824,
                        -8379,   -16675,   -12831,  -10445,
                       -11650,    14403,    -7003,   10056,
                         6322,    -9884,     5751,   -4950,
                         4130,    -3958,     3258,    2616,
                        -1897,    -2117,     2354,   -1423,
                        -1117,    -1571,    -1739,   -4421,
                         1165,     8752)

# Set lunar distance masks for vectorised operations.
G.mask.distance.D <- c(0,
                       2,
                       2,
                       0,
                       0,
                       0,
                       2,
                       2,
                       2,
                       2,
                       0,
                       1,
                       0,
                       2,
                       0,
                       4,
                       0,
                       4,
                       2,
                       2,
                       1,
                       1,
                       2,
                       2,
                       4,
                       2,
                       0,
                       2,
                       1,
                       2,
                       0,
                       2,
                       2,
                       4,
                       3,
                       2,
                       4,
                       0,
                       2,
                       4,
                       0,
                       4,
                       1,
                       0,
                       0,
                       2)
#
G.mask.distance.F.x.2 <- c( 0,
                            0,
                            0,
                            0,
                            0,
                            1,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                           -1,
                           -1,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                           -1,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                           -1,
                            0,
                           -1)
#
G.mask.distance.M <- c( 0,
                        0,
                        0,
                        0,
                        1,
                        0,
                        0,
                       -1,
                        0,
                       -1,
                        1,
                        0,
                        1,
                        0,
                        0,
                        0,
                        0,
                        0,
                        1,
                        1,
                        0,
                        1,
                       -1,
                        0,
                        0,
                        0,
                        1,
                       -1,
                        0,
                       -2,
                        1,
                       -2,
                        0,
                       -1,
                        0,
                        1,
                       -1,
                        2,
                        2,
                        0,
                        0,
                       -1,
                        0,
                        0,
                        2,
                        0)
#
G.mask.distance.M.prime <- c( 1,
                             -1,
                              0,
                              2,
                              0,
                              0,
                             -2,
                             -1,
                              1,
                              0,
                             -1,
                              0,
                              1,
                              0,
                              1,
                             -1,
                              3,
                             -2,
                             -1,
                              0,
                             -1,
                              0,
                              1,
                              2,
                              0,
                             -3,
                             -2,
                             -2,
                              1,
                              0,
                              2,
                             -1,
                              1,
                             -1,
                             -1,
                              1,
                             -2,
                             -1,
                             -1,
                              1,
                              4,
                              0,
                             -2,
                              2,
                              1,
                             -1)

# Lunar phase (full moon) calculation coefficients.
G.vec.ph.full.coeff <- c(mpfr("-0.40614", 128),
                         mpfr( "0.17302", 128),
                         mpfr( "0.01614", 128),
                         mpfr( "0.01043", 128),
                         mpfr( "0.00734", 128),
                         mpfr("-0.00515", 128),
                         mpfr( "0.00209", 128),
                         mpfr("-0.00111", 128),
                         mpfr("-0.00057", 128),
                         mpfr( "0.00056", 128),
                         mpfr("-0.00042", 128),
                         mpfr( "0.00042", 128),
                         mpfr( "0.00038", 128),
                         mpfr("-0.00024", 128),
                         mpfr("-0.00017", 128),
                         mpfr("-0.00007", 128),
                         mpfr( "0.00004", 128),
                         mpfr( "0.00004", 128),
                         mpfr( "0.00003", 128),
                         mpfr( "0.00003", 128),
                         mpfr("-0.00003", 128),
                         mpfr( "0.00003", 128),
                         mpfr("-0.00002", 128),
                         mpfr("-0.00002", 128),
                         mpfr( "0.00002", 128))

# Lunar phase (new moon) calculation coefficients.
G.vec.ph.new.coeff <- c(mpfr("-0.40720", 128),
                        mpfr( "0.17241", 128),
                        mpfr( "0.01608", 128),
                        mpfr( "0.01039", 128),
                        mpfr( "0.00739", 128),
                        mpfr("-0.00514", 128),
                        mpfr( "0.00208", 128),
                        mpfr("-0.00111", 128),
                        mpfr("-0.00057", 128),
                        mpfr( "0.00056", 128),
                        mpfr("-0.00042", 128),
                        mpfr( "0.00042", 128),
                        mpfr( "0.00038", 128),
                        mpfr("-0.00024", 128),
                        mpfr("-0.00017", 128),
                        mpfr("-0.00007", 128),
                        mpfr( "0.00004", 128),
                        mpfr( "0.00004", 128),
                        mpfr( "0.00003", 128),
                        mpfr( "0.00003", 128),
                        mpfr("-0.00003", 128),
                        mpfr( "0.00003", 128),
                        mpfr("-0.00002", 128),
                        mpfr("-0.00002", 128),
                        mpfr( "0.00002", 128))

# Lunar phase (quarter phase) calculation coefficients.
G.vec.ph.qrt.coeff <- c(mpfr("-0.62801", 128),
                        mpfr( "0.17172", 128),
                        mpfr("-0.01183", 128),
                        mpfr( "0.00862", 128),
                        mpfr( "0.00804", 128),
                        mpfr( "0.00454", 128),
                        mpfr( "0.00204", 128),
                        mpfr("-0.00180", 128),
                        mpfr("-0.00070", 128),
                        mpfr("-0.00040", 128),
                        mpfr("-0.00034", 128),
                        mpfr( "0.00032", 128),
                        mpfr( "0.00032", 128),
                        mpfr("-0.00028", 128),
                        mpfr( "0.00027", 128),
                        mpfr("-0.00017", 128),
                        mpfr("-0.00005", 128),
                        mpfr( "0.00004", 128),
                        mpfr("-0.00004", 128),
                        mpfr( "0.00004", 128),
                        mpfr( "0.00003", 128),
                        mpfr( "0.00003", 128),
                        mpfr( "0.00002", 128),
                        mpfr( "0.00002", 128),
                        mpfr("-0.00002", 128))

# Lunar phase (quarter phase) supplemental corrections.
G.vec.ph.sup.coeff <- c(mpfr("-0.00038", 128),
                        mpfr(" 0.00026", 128),
                        mpfr("-0.00002", 128),
                        mpfr(" 0.00002", 128),
                        mpfr(" 0.00002", 128))

# Lunar phase additional corrections coefficients.
G.vec.ph.add.coeff.outer <- c(mpfr("0.000325", 128),
                              mpfr("0.000165", 128),
                              mpfr("0.000164", 128),
                              mpfr("0.000126", 128),
                              mpfr("0.000110", 128),
                              mpfr("0.000062", 128),
                              mpfr("0.000060", 128),
                              mpfr("0.000056", 128),
                              mpfr("0.000047", 128),
                              mpfr("0.000042", 128),
                              mpfr("0.000040", 128),
                              mpfr("0.000037", 128),
                              mpfr("0.000035", 128),
                              mpfr("0.000023", 128))
#
G.vec.ph.add.addend      <- c(mpfr("299.77", 128),
                              mpfr("251.88", 128),
                              mpfr("251.83", 128),
                              mpfr("349.42", 128),
                              mpfr( "84.66", 128),
                              mpfr("141.74", 128),
                              mpfr("207.14", 128),
                              mpfr("154.84", 128),
                              mpfr( "34.52", 128),
                              mpfr("207.19", 128),
                              mpfr("291.34", 128),
                              mpfr("161.72", 128),
                              mpfr("239.56", 128),
                              mpfr("331.55", 128))
#
G.vec.ph.add.coeff.inner <- c(mpfr( "0.107408", 128),
                              mpfr( "0.016321", 128),
                              mpfr("26.651886", 128),
                              mpfr("36.412478", 128),
                              mpfr("18.206239", 128),
                              mpfr("53.303771", 128),
                              mpfr( "2.453732", 128),
                              mpfr( "7.306860", 128),
                              mpfr("27.261239", 128),
                              mpfr( "0.121824", 128),
                              mpfr( "1.844379", 128),
                              mpfr("24.198154", 128),
                              mpfr("25.513099", 128),
                              mpfr( "3.592518", 128))
#
G.vec.ph.add.subtrahend  <- c(mpfr("0.009173", 128),
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr,
                              G.zero.mpfr)
```

``` r
# Precision Testing
# -----------------
#
# Quantify the effect of using R's default double (53-bit) precision vs mpfr 128-bit precision.
# Here, converting a date-time in JDE (Julian Ephemeris Day) to Dynamical Time is the test case.

G.lunar.apsis <- "perigee"
G.lunar.phase <- "full moon"

JDE.A <- 2443259.9
JDE.B <- mpfr("2443259.9", 128)

result.A <- TimeJDEToDynamical(JDE.A)
result.B <- TimeJDEToDynamical(JDE.B)

cat("Precision Testing\n-----------------\n\n", sep="")

cat("JDE, low precision...\n", sep="")
print(JDE.A, digits = 17)
cat("\nJDE, high precision...\n", sep="")
print(JDE.B, digits = 17)

cat("\n\nJDE --> Dynamical time, low precision (top) vs high precision (bottom)...\n", sep="")
print(result.A, digits = 17)
print(result.B, digits = 17)

cat("\nWhat's the corresponding difference in lunar distance?")
lunar.distance.A <- GetLunarDistance(JDE.A)
lunar.distance.B <- GetLunarDistance(JDE.B)
difference <- abs(lunar.distance.A$value - lunar.distance.B$value) * mpfr("1000000", 128)
cat("\n\U2248", formatMpfr(difference, 8), " mm\n", sep="")
```

    ## Precision Testing
    ## -----------------
    ## 
    ## JDE, low precision...
    ## [1] 2443259.8999999999
    ## 
    ## JDE, high precision...
    ## 1 'mpfr' number of precision  128   bits 
    ## [1] 2443259.9
    ## 
    ## 
    ## JDE --> Dynamical time, low precision (top) vs high precision (bottom)...
    ## [1] "1977-04-26 09:35:59 UTC"
    ## [1] "1977-04-26 09:36:00 UTC"
    ## 
    ## What's the corresponding difference in lunar distance?
    ## ˜0.46645042 mm

``` r
# Testing
# -------

cat("Testing\n-------\n\n", sep="")
cat("Test the system's functions against examples in 'Astronomical Algorithms'.\n")
cat("Ref: Meeus J., Astronomical Algorithms, 2nd Ed, Willmann-Bell, Richmond, VA, 1998.\n")

######
cat("\n-------\nTesting: New Moon Feb 1977 (page 353)\n\n")
G.lunar.phase <- "new moon"
date.target.UTC <- date_decimal(1977.13)
date.target.dynamical <- TimeUniversalToDynamical(date.target.UTC)
phase.date <- CalculateWithUncertainty(GetLunarPhaseJDE, date.target.dynamical$value, date.target.dynamical$uncertainty, TRUE)
label <- paste("Next ", G.lunar.phase, ":", sep="")
cat(label, replicate(G.max.indent - nchar(label), " "), formatMpfr(phase.date$value, 11), " +- ", formatMpfr(phase.date$uncertainty, 1), " JDE\n", sep="")
phase.date.dynamical <- CalculateWithUncertainty(TimeJDEToDynamical, phase.date$value, phase.date$uncertainty, FALSE)
cat(replicate(G.max.indent - 2, " "), "= ", format(round_date(phase.date.dynamical$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ", format(hms::as_hms(signif(as.numeric(phase.date.dynamical$uncertainty, units = "secs"), 2)), "%H:%M:%S"), " Dynamical Time\n", sep="")
phase.date.UTC <- CalculateWithUncertainty(TimeDynamicalToUniversal, phase.date.dynamical$value, phase.date.dynamical$uncertainty, TRUE)
cat(replicate(G.max.indent - 2, " "), "= ", format(round_date(phase.date.UTC$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ", format(hms::as_hms(signif(as.numeric(phase.date.UTC$uncertainty, units = "secs"), 2)), "%H:%M:%S"), " UTC\n\n", sep="")

label <- paste("Expected:")
cat(label, replicate(G.max.indent - nchar(label), " "), "18-Feb-1977 03:37 UT.\n", sep="")
cat(replicate(G.max.indent, " "), "N.B. UT = UTC +- 0.9 s.\n", sep="")
######


######
cat("\n-------\nTesting: Last Quarter Moon Jan 2044 (page 353)\n\n")
G.lunar.phase <- "last quarter"
date.target.UTC <- date_decimal(2044.026)
date.target.dynamical <- TimeUniversalToDynamical(date.target.UTC)
phase.date <- CalculateWithUncertainty(GetLunarPhaseJDE, date.target.dynamical$value, date.target.dynamical$uncertainty, TRUE)
label <- paste("Next ", G.lunar.phase, ":", sep="")
cat(label, replicate(G.max.indent - nchar(label), " "), formatMpfr(phase.date$value, 11), " +- ", formatMpfr(phase.date$uncertainty, 1), " JDE\n", sep="")
phase.date.dynamical <- CalculateWithUncertainty(TimeJDEToDynamical, phase.date$value, phase.date$uncertainty, FALSE)
cat(replicate(G.max.indent - 2, " "), "= ", format(round_date(phase.date.dynamical$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ", format(hms::as_hms(signif(as.numeric(phase.date.dynamical$uncertainty, units = "secs"), 2)), "%H:%M:%S"), " Dynamical Time\n\n", sep="")

label <- paste("Expected:")
cat(label, replicate(G.max.indent - nchar(label), " "), "21-Jan-2044 23:48:17 Dynamical Time.\n", sep="")
######


######
cat("\n-------\nTesting: Lunar Apogee Oct 1988 (page 357)\n\n")
G.lunar.apsis <- "apogee"
date.target.UTC <- date_decimal(1988.75)
date.target.dynamical <- TimeUniversalToDynamical(date.target.UTC)
apsis.date <- CalculateWithUncertainty(GetLunarApsisJDE, date.target.dynamical$value, date.target.dynamical$uncertainty, TRUE)
label <- paste("Next ", G.lunar.apsis, ":", sep="")
cat(label, replicate(G.max.indent - nchar(label), " "), formatMpfr(apsis.date$value, 10), " +- ", formatMpfr(apsis.date$uncertainty, 1), " JDE\n", sep="")
apsis.date.dynamical <- CalculateWithUncertainty(TimeJDEToDynamical, apsis.date$value, apsis.date$uncertainty, FALSE)
cat(replicate(G.max.indent - 2, " "), "= ", format(round_date(apsis.date.dynamical$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ", format(hms::as_hms(signif(as.numeric(apsis.date.dynamical$uncertainty, units = "secs"), 3)), "%H:%M:%S"), " Dynamical Time\n\n", sep="")

label <- paste("Expected:")
cat(label, replicate(G.max.indent - nchar(label), " "), "07-Oct-1988 20:29 Dynamical Time.\n", sep="")
######


######
cat("\n-------\nTesting: Lunar Distance for 12 April 1992, at 00:00:00 Dynamical Time (pages 342-343).\n")
lunar.distance <- CalculateWithUncertainty(GetLunarDistance, mpfr("2448724.5", 128), G.zero.mpfr, TRUE)
label <- paste("Estimated lunar distance:", sep="")
cat("\n", label, replicate(G.max.indent - nchar(label), " "), signif(asNumeric(lunar.distance$value), 6), " +- ", signif(asNumeric(lunar.distance$uncertainty), 1), " km\n\n", sep="")

label <- paste("Expected:")
cat(label, replicate(G.max.indent - nchar(label), " "), "368409.7 km.\n", sep="")
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
    ##                          = 07-Oct-1988 20:30:12 +- 00:03:00 Dynamical Time
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
# Compare this system vs astropixels.com and fourmilab.com.

G.lunar.apsis <- "perigee"
G.lunar.phase <- "full moon"

label.date <- "Target date:"
chars.label.date <- nchar(label.date)
label.distance <- "Lunar distance"

vec.dates.astropixels <- c(dmy_hms("09-01-2001 20:24:00", tz="UTC"),
                           dmy_hms("04-07-2050 18:51:00", tz="UTC"),
                           dmy_hms("23-05-2100 17:25:00", tz="UTC"))

vec.dates.fourmilab   <- c(dmy_hms("10-01-2001 09:00:00", tz="UTC"),
                           dmy_hms("07-07-2050 02:26:00", tz="UTC"),
                           dmy_hms("22-05-2100 11:08:00", tz="UTC"))

vec.distances.astropixels <- c(357406, 367058, 360904)
vec.distances.fourmilab   <- c(357131, 363255, 359497)

# Get lunar distance at target dates, for astropixels dates.
# N.B. Convert dates from UTC --> Dynamical Time--> JDE.
vec.calculated.distances <- lapply(
                               lapply(
                                  lapply(
                                     vec.dates.astropixels, TimeUniversalToDynamical),
                                  function(x) {CalculateWithUncertainty(TimeDynamicalToJDE, x$value, x$uncertainty, FALSE)}
                               ),
                               function(x) {CalculateWithUncertainty(GetLunarDistance, x$value, x$uncertainty, TRUE)}
                             )

# Display astropixels.com data.
cat("Lunar distance benchmarking: astropixels.com\n\n")
for (i in 1:length(vec.dates.astropixels)) {
   
   # Display target date.
   cat("------\n", label.date, replicate(G.max.indent - chars.label.date, " "), format(vec.dates.astropixels[i], "%d-%b-%Y %H:%M:%S %Z"), "\n", sep="")

   # Display astropixels.com's estimate.
   cat("\n", label.distance, sep="")
   cat("\nEstimate, astropixels.com:", replicate(G.max.indent - nchar("Estimate, astropixels.com:"), " "), vec.distances.astropixels[i], " +- <unspecified> km", sep="")

   # Display lunar distance at target date.
   cat("\nEstimate, R:", replicate(G.max.indent - nchar("Estimate, R:"), " "), signif(asNumeric(vec.calculated.distances[[i]]$value), 6), " +- ", signif(asNumeric(vec.calculated.distances[[i]]$uncertainty), 1), " km\n\n", sep="")
}


# Get lunar distance at target dates, for fourmilab dates.
# N.B. That source states that its times are only accurate to +- 2 minutes.
#      That error is incorporated, below.
vec.calculated.distances <- lapply(
                               lapply(
                                  lapply(
                                     vec.dates.fourmilab, function(x) {CalculateWithUncertainty(TimeUniversalToDynamical, x, dminutes(mpfr("2", 128)), TRUE)}),
                                  function(x) {CalculateWithUncertainty(TimeDynamicalToJDE, x$value, x$uncertainty, FALSE)}
                               ),
                               function(x) {CalculateWithUncertainty(GetLunarDistance, x$value, x$uncertainty, TRUE)}
                             )

# Display fourmilab.com data.
# N.B. The first test distance needs a different number of decimal places than the rest.
cat(replicate(G.max.indent, "-"), sep="")
cat("\n\nLunar distance benchmarking: www.fourmilab.ch\n\n")
#
# Display target date.
cat("------\n", label.date, replicate(G.max.indent - chars.label.date, " "), format(vec.dates.fourmilab[1], "%d-%b-%Y %H:%M:%S"), " +- 00:00:02 UTC\n", sep="")

# Display fourmilab.com's estimate.
cat("\n", label.distance, sep="")
cat("\nEstimate, fourmilab.com:", replicate(G.max.indent - nchar("Estimate, fourmilab.com:"), " "), format(vec.distances.fourmilab[1], nsmall=1), " +- ", signif(asNumeric(vec.calculated.distances[[1]]$uncertainty), 1), " km (assumed uncertainty, from ELP 2000-82)", sep="")

# Display lunar distance at target date.
cat("\nEstimate, R:", replicate(G.max.indent - nchar("Estimate, R:"), " "), format(signif(asNumeric(vec.calculated.distances[[1]]$value), 7), nsmall=1), " +- ", signif(asNumeric(vec.calculated.distances[[1]]$uncertainty), 1), " km\n\n", sep="")
#
for (i in 2:length(vec.dates.fourmilab)) {

   # Display target date.
   cat("------\n", label.date, replicate(G.max.indent - chars.label.date, " "), format(vec.dates.fourmilab[i], "%d-%b-%Y %H:%M:%S"), " +- 00:00:02 UTC\n", sep="")

   # Display fourmilab.com's estimate.
   cat("\n", label.distance, sep="")
   cat("\nEstimate, fourmilab.com:", replicate(G.max.indent - nchar("Estimate, fourmilab.com:"), " "), vec.distances.fourmilab[i], " +- ", signif(asNumeric(vec.calculated.distances[[i]]$uncertainty), 1), " km (assumed uncertainty, from ELP 2000-82)", sep="")

   # Display lunar distance at target date.
   cat("\nEstimate, R:", replicate(G.max.indent - nchar("Estimate, R:"), " "), signif(asNumeric(vec.calculated.distances[[i]]$value), 6), " +- ", signif(asNumeric(vec.calculated.distances[[i]]$uncertainty), 1), " km\n\n", sep="")
}
```

    ## Lunar distance benchmarking: astropixels.com
    ## 
    ## ------
    ## Target date:               09-Jan-2001 20:24:00 UTC
    ## 
    ## Lunar distance
    ## Estimate, astropixels.com: 357406 +- <unspecified> km
    ## Estimate, R:               357407 +- 1 km
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

cat("Example 1\n---------\n\n", sep="")
cat("For a target date, determine the time of the next....\n")
cat("a) full moon,\n")
cat("b) perigee and associated lunar distance.\n\n")

G.lunar.apsis <- "perigee"
G.lunar.phase <- "full moon"

# Set target date.
# Assume format is dd-mm-yyyy hh:mm:ss, with UTC time zone.
date.target.UTC <- dmy_hms("18-02-2029 00:00:00", tz="UTC")
label <- "Target date:"
cat(label, replicate(G.max.indent - nchar(label), " "), format(date.target.UTC, "%d-%b-%Y %H:%M:%S %Z"), "\n", sep="")

# Ensure a 21st century date is being processed.
if (!IsValidDate(date.target.UTC)) {
   cat("\nInvalid date; date.target.UTC:\n")
   print(date.target.UTC, digits=20)
   cat("\nNeed 21st century dates.")
   stopifnot(FALSE)
}

# Get the dynamical date equivalent of the target date.
date.target.dynamical <- TimeUniversalToDynamical(date.target.UTC)
cat(replicate(G.max.indent - 2, " "), "= ", format(round_date(date.target.dynamical$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ", format(signif(date.target.dynamical$uncertainty, 1), "%H:%M:%S"), " Dynamical Time\n\n", sep="")

# Get the date of the next nominated lunar phase (set to full moon by default).
phase.date <- CalculateWithUncertainty(GetLunarPhaseJDE, date.target.dynamical$value, date.target.dynamical$uncertainty, TRUE)
label <- paste("Next ", G.lunar.phase, ":", sep="")
cat(label, replicate(G.max.indent - nchar(label), " "), format(signif(asNumeric(phase.date$value), 11), nsmall=4), " +- ", signif(asNumeric(phase.date$uncertainty), 1), " JDE\n", sep="")

# Convert JDE to dynamical time.
# N.B. Handle uncertainty differently than for target date; it used a decimal uncertainty, this uses a difftime value.
phase.date.dynamical <- CalculateWithUncertainty(TimeJDEToDynamical, phase.date$value, phase.date$uncertainty, FALSE)
cat(replicate(G.max.indent - 2, " "), "= ", format(round_date(phase.date.dynamical$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ", format(hms::as_hms(signif(as.numeric(phase.date.dynamical$uncertainty, units = "secs"), 2)), "%H:%M:%S"), " Dynamical Time\n", sep="")

# Convert dynamical time to UTC.
# N.B. As for the JDE to dynamical time conversion, the uncertainty here uses a difftime value.
phase.date.UTC <- CalculateWithUncertainty(TimeDynamicalToUniversal, phase.date.dynamical$value, phase.date.dynamical$uncertainty, TRUE)
cat(replicate(G.max.indent - 2, " "), "= ", format(round_date(phase.date.UTC$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ", format(hms::as_hms(signif(as.numeric(phase.date.UTC$uncertainty, units = "secs"), 2)), "%H:%M:%S"), " UTC\n", sep="")


# Get the perigee estimation.
apsis.date <- CalculateWithUncertainty(GetLunarApsisJDE, date.target.dynamical$value, date.target.dynamical$uncertainty, TRUE)
label <- paste("Next ", G.lunar.apsis, ":", sep="")
cat("\n", label, replicate(G.max.indent - nchar(label), " "), format(signif(asNumeric(apsis.date$value), 9), nsmall=2), " +- ", signif(asNumeric(apsis.date$uncertainty), 1), " JDE\n", sep="")

# Convert JDE to dynamical time.
# N.B. Handle uncertainty differently than for target date; it used a decimal uncertainty, this uses a difftime value.
apsis.date.dynamical <- CalculateWithUncertainty(TimeJDEToDynamical, apsis.date$value, apsis.date$uncertainty, FALSE)
cat(replicate(G.max.indent - 2, " "), "= ", format(round_date(apsis.date.dynamical$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ", format(hms::as_hms(signif(as.numeric(apsis.date.dynamical$uncertainty, units = "secs"), 4)), "%H:%M:%S"), " Dynamical Time\n", sep="")

# Convert dynamical time to UTC.
# N.B. As for the JDE to dynamical time conversion, the uncertainty here uses a difftime value.
apsis.date.UTC <- CalculateWithUncertainty(TimeDynamicalToUniversal, apsis.date.dynamical$value, apsis.date.dynamical$uncertainty, TRUE)
cat(replicate(G.max.indent - 2, " "), "= ", format(round_date(apsis.date.UTC$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ", format(hms::as_hms(signif(as.numeric(apsis.date.UTC$uncertainty, units = "secs"), 4)), "%H:%M:%S"), " UTC\n", sep="")


# Get lunar distance at target apsis.
lunar.distance <- CalculateWithUncertainty(GetLunarDistance, apsis.date$value, apsis.date$uncertainty, TRUE)
label <- paste("Lunar distance at ", G.lunar.apsis, ":", sep="")
cat("\n", label, replicate(G.max.indent - nchar(label), " "), signif(asNumeric(lunar.distance$value), 6), " +- ", signif(asNumeric(lunar.distance$uncertainty), 1), " km\n", sep="")
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
    ##                          = 01-Mar-2029 18:29:55 +- 00:31:19 UTC
    ## 
    ## Lunar distance at perigee: 358630 +- 3 km

``` r
# Example 2a
# ----------
#
# Graph variation in lunar distance for 2019-2020.
# Overlay full and new moons.
# First step: generate required data.

cat("Example 2\n---------\n\n", sep="")
cat("Graph variation in lunar distance for 2019-2020.\n")
cat("Overlay full and new moons.\n")

# Prepare poster data.
date.start <- dmy_hms("01-01-2019 00:00:00", tz="UTC")
date.end <- dmy_hms("31-12-2020 00:00:00", tz="UTC")

# Series 1: Get all full moons and their distances.
G.lunar.phase <- "full moon"
list.full.moon <- GetListEstimates(GetLunarPhaseJDE, date.start, date.end, mpfr("29.4", 128))
#
list.full.moon.distances <- lapply(list.full.moon,
   function(x) {CalculateWithUncertainty(GetLunarDistance, x$value, x$uncertainty, TRUE)})
#
# Convert full moon JDE dates to UTC.
list.full.moon <- lapply(list.full.moon,
   function(x) {
      target.Dynamical <- CalculateWithUncertainty(TimeJDEToDynamical, x$value, x$uncertainty, FALSE)
      CalculateWithUncertainty(TimeDynamicalToUniversal, target.Dynamical$value, target.Dynamical$uncertainty, TRUE)
})


# Series 2: Get all new moons and their distances.
G.lunar.phase <- "new moon"
list.new.moon <- GetListEstimates(GetLunarPhaseJDE, date.start, date.end, mpfr("29.4", 128))
#
list.new.moon.distances <- lapply(list.new.moon,
   function(x) {CalculateWithUncertainty(GetLunarDistance, x$value, x$uncertainty, TRUE)})
#
# Convert new moon JDE dates to UTC.
list.new.moon <- lapply(list.new.moon,
   function(x) {
      target.Dynamical <- CalculateWithUncertainty(TimeJDEToDynamical, x$value, x$uncertainty, FALSE)
      CalculateWithUncertainty(TimeDynamicalToUniversal, target.Dynamical$value, target.Dynamical$uncertainty, TRUE)
})


# Series 3: Get enough other lunar distances to graph. Try a point every 0.5 days.
interval.days <- ddays(mpfr("0.5", 128))
list.filler.dates <- 0:{ceiling({date.end - date.start} / interval.days) - 1}
list.filler.dates <- date.start + list.filler.dates * interval.days
list.filler.distances <- lapply(list.filler.dates,
   function(x) {
      target.Dynamical <- TimeUniversalToDynamical(x)
      target.JDE <- CalculateWithUncertainty(TimeDynamicalToJDE, target.Dynamical$value, target.Dynamical$uncertainty, FALSE)
      CalculateWithUncertainty(GetLunarDistance, target.JDE$value, target.JDE$uncertainty, TRUE)
})


# Convert lists to dataframes.
#
# Full moon UTC dates (values and uncertainties).
list.full.moon.dates <- lapply(list.full.moon, function(x) {x$value})
list.full.moon.dates.uncertainties <- lapply(list.full.moon, function(x) {x$uncertainty})
#
# Full moon distances (values and uncertainties).
list.full.moon.distances.approx <- lapply(list.full.moon.distances, function(x) {suppress_warnings(asNumeric(x$value), "*")})
list.full.moon.distances.uncertainties.approx <- lapply(list.full.moon.distances, function(x) {suppress_warnings(asNumeric(x$uncertainty), "*")})
#
# New moon UTC dates (values and uncertainties).
list.new.moon.dates <- lapply(list.new.moon, function(x) {x$value})
list.new.moon.dates.uncertainties <- lapply(list.new.moon, function(x) {x$uncertainty})
#
# New moon distances (values and uncertainties).
list.new.moon.distances.approx <- lapply(list.new.moon.distances, function(x) {suppress_warnings(asNumeric(x$value), "*")})
list.new.moon.distances.uncertainties.approx <- lapply(list.new.moon.distances, function(x) {suppress_warnings(asNumeric(x$uncertainty), "*")})
#
# Filler distances (values and uncertainties).
list.filler.distances.approx <- lapply(list.filler.distances, function(x) {suppress_warnings(asNumeric(x$value), "*")})
list.filler.distances.uncertainties.approx <- lapply(list.filler.distances, function(x) {suppress_warnings(asNumeric(x$uncertainty), "*")})
#
lunar.data <- cbind(data.frame(day = as_datetime(unlist(list.filler.dates)),
                          distance = unlist(list.filler.distances.approx),
              distance.uncertainty = unlist(list.filler.distances.uncertainties.approx))
                   )

full.moon.data <- cbind(data.frame(day.full.moon = as_datetime(unlist(list.full.moon.dates)),
                       day.full.moon.uncertainty = as_datetime(unlist(list.full.moon.dates.uncertainties)),
                              distance.full.moon = unlist(list.full.moon.distances.approx),
                  distance.full.moon.uncertainty = unlist(list.full.moon.distances.uncertainties.approx))
                       )

new.moon.data <- cbind(data.frame(day.new.moon = as_datetime(unlist(list.new.moon.dates)),
                      day.new.moon.uncertainty = as_datetime(unlist(list.new.moon.dates.uncertainties)),
                             distance.new.moon = unlist(list.new.moon.distances.approx),
                 distance.new.moon.uncertainty = unlist(list.new.moon.distances.uncertainties.approx))
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

# Load needed font, if needed, but only once.
# N.B. Toggle G.extra.font.installed as needed.
# N.B. It may need console-based confirmation.

# Import Candara font, if not already available.
G.extra.font.installed <- TRUE
if (!G.extra.font.installed) {
   font_import(pattern = "Candara.ttf")
   loadfonts(device = "win")
   fonts()

   G.extra.font.installed <- TRUE
}

# Compose graph.
p <- ggplot() +
        geom_line(data = lunar.data, aes(x = day, y = distance, color = "distance")) +
        geom_point(data = full.moon.data, aes(x = day.full.moon, y = distance.full.moon, color = "full moon"), shape = 16, size = 4) +
        geom_point(data = new.moon.data, aes(x = day.new.moon, y = distance.new.moon, color = "new moon"), shape = 4, size = 4, stroke = 2) +
        ggtitle("Lunar Distance vs. Date", subtitle = "(with full and new moons overlaid)") + 
        guides(color = guide_legend(override.aes = list(shape = c(NA, 16, 4)))) + 
        scale_color_manual(values = c("#18677E", "#18677E", "#18677E"), name = "Legend") + 
        scale_x_datetime(labels = date_format("%d-%b-%Y")) + 
        scale_y_continuous(labels = comma) +
        theme_bw() + 
        theme_classic() + 
        theme(axis.text.x = element_text(color="#000000"),
              axis.text.y = element_text(color="#000000"),
              legend.background = element_rect(color = NA, fill = "transparent"),
              legend.box.background = element_rect(color = NA, fill = "transparent"),
              panel.background = element_rect(fill = "transparent"),
              plot.background = element_rect(fill = "transparent", color = NA),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              text = element_text(color = "#717171", family = "Candara", size = 23.76097441)) +
        xlab("\ndate (UTC)") +
        ylab("distance (km)\n")

# Save graph, large format for poster use.
ggsave("supermoon_graph.png", bg = "transparent", width = 11.9, height = 8.499, units = "in")
```

``` r
# Example 2c
# ----------
#
# Graph variation in lunar distance for 2019-2020.
# Overlay full and new moons.
# Third step: display research paper-style graph.

# Compose graph.
p <- ggplot() +
        geom_line(data = lunar.data, aes(x = day, y = distance, color = "distance")) +
        geom_point(data = full.moon.data, aes(x = day.full.moon, y = distance.full.moon, color = "full moon"), shape = 16) +
        geom_point(data = new.moon.data, aes(x = day.new.moon, y = distance.new.moon, color = "new moon"), shape = 4) + 
        ggtitle("Lunar Distance vs. Date", subtitle = "(with full and new moons overlaid)") + 
        guides(color = guide_legend(override.aes = list(shape = c(NA, 16, 4)))) + 
        scale_color_manual(values = c("#000000", "#000000", "#000000"), name = "Legend") + 
        scale_x_datetime(labels = date_format("%d-%b-%Y")) + 
        scale_y_continuous(labels = comma) +
        theme_bw() + 
        theme_classic() + 
        theme(axis.text.x = element_text(color="#000000"),
          axis.text.y = element_text(color="#000000"),
              legend.background = element_rect(color = NA, fill = "transparent"),
              legend.box.background = element_rect(color = NA, fill = "transparent"),
              panel.background = element_rect(fill = "transparent"),
              plot.background = element_rect(fill = "transparent", color = NA),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              text = element_text(color = "#000000", family = "Times New Roman")) +
        xlab("\ndate (UTC)") +
        ylab("distance (km)\n")

# Display graph.
p
```

![](graph_distance_vs_date.png)<!-- -->

``` r
# Save graph, small format for research paper use.
ggsave(file = "supermoon_graph.eps", device = cairo_ps, family = "Times", scale = 0.7)
```

    ## Saving 4.9 x 3.5 in image

``` r
# Example 3
# ---------
#
# Get all 21st century supermoons.
# Define "supermoon" as a full moon and lunar perigee separated by <= 12 h.

cat("Example 3\n---------\n\n", sep="")
cat("Get all 21st century supermoons.\n")
cat("Define 'supermoon' as a full moon and lunar perigee separated by <= 12 h.\n")

G.lunar.apsis <- "perigee"
G.lunar.phase <- "full moon"

# Get the 21st century full moons.
list.full.moons <- GetListEstimates(GetLunarPhaseJDE, dmy_hms("01-01-2001 00:00:00", tz="UTC"), dmy_hms("31-12-2100 00:00:00", tz="UTC"), mpfr("29.4", 128))

# Get the 21st century lunar perigees.
list.lunar.perigees <- GetListEstimates(GetLunarApsisJDE, dmy_hms("01-01-2001 00:00:00", tz="UTC"), dmy_hms("31-12-2100 00:00:00", tz="UTC"), mpfr("27.4", 128))

# Build the supermoons list.
# Each list element will hold: (date-time of full moon, date-time of perigee, distance).
# N.B. A 12 h UTC or DT time difference = 0.5 JDE.
list.supermoons <- list()
num.full.moons <- length(list.full.moons)
num.lunar.perigees <- length(list.lunar.perigees)
match.found <- FALSE
#
for (x in 1:num.full.moons) {
   for (y in 1:num.lunar.perigees) {

      # Only assess, in detail, relatively close data points.
      if (abs(list.full.moons[[x]]$value - list.lunar.perigees[[y]]$value) <= 2) {

         if (list.full.moons[[x]]$value < list.lunar.perigees[[y]]$value) {
            if ({list.lunar.perigees[[y]]$value - list.lunar.perigees[[y]]$uncertainty} - {list.full.moons[[x]]$value + list.full.moons[[x]]$uncertainty} <= 0.5) {
               match.found <- TRUE
               break
            }

         } else{
            if ({list.full.moons[[x]]$value - list.full.moons[[x]]$uncertainty} - {list.lunar.perigees[[y]]$value + list.lunar.perigees[[y]]$uncertainty} <= 0.5) {
               match.found <- TRUE
               break
            }
         }

      }
   }

   if (match.found) {
      # a) Convert time of full moon from JDE to UTC.
      # b) Convert time of perigee from JDE to UTC.
      # c) Get lunar distance at perigee.

      phase.date.dynamical <- CalculateWithUncertainty(TimeJDEToDynamical, list.full.moons[[x]]$value, list.full.moons[[x]]$uncertainty, FALSE)

      apsis.date.dynamical <- CalculateWithUncertainty(TimeJDEToDynamical, list.lunar.perigees[[y]]$value, list.lunar.perigees[[y]]$uncertainty, FALSE)

      list.supermoons[[length(list.supermoons) + 1]] <- list(
         "full.moon" = CalculateWithUncertainty(TimeDynamicalToUniversal, phase.date.dynamical$value, phase.date.dynamical$uncertainty, TRUE),
         "perigee" = CalculateWithUncertainty(TimeDynamicalToUniversal, apsis.date.dynamical$value, apsis.date.dynamical$uncertainty, TRUE),
         "distance" = CalculateWithUncertainty(GetLunarDistance, list.lunar.perigees[[y]]$value, list.lunar.perigees[[y]]$uncertainty, TRUE))

      # Delete the perigee from the perigee list.
      list.lunar.perigees[- y]
      num.lunar.perigees <- length(list.lunar.perigees)

      match.found <- FALSE
   }
}


# Display results.
cat("\n---------------------------------------------\n")
#
for (x in 1:length(list.supermoons)) {

   cat("supermoon #", x, "\n", sep="")

   label <- paste(" ", G.lunar.phase, ":", sep="")
   if (list.supermoons[[x]]$full.moon$uncertainty < 60) {
      significant.figures <- 2
   } else {
      significant.figures <- 3
   }
   cat(label, replicate(G.max.indent - nchar(label), " "), format(round_date(list.supermoons[[x]]$full.moon$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ",
       format(hms::as_hms(signif(as.numeric(list.supermoons[[x]]$full.moon$uncertainty, units = "secs"), significant.figures)), "%H:%M:%S"), " UTC\n", sep="")

   label <- paste(" ", G.lunar.apsis, ":", sep="")
   cat(label, replicate(G.max.indent - nchar(label), " "), format(round_date(list.supermoons[[x]]$perigee$value, unit = "seconds"), "%d-%b-%Y %H:%M:%S"), " +- ",  
       format(hms::as_hms(signif(as.numeric(list.supermoons[[x]]$perigee$uncertainty, units = "secs"), 4)), "%H:%M:%S"), " UTC\n", sep="")

   label <- paste(" distance at ", G.lunar.apsis, ":", sep="")
   if (list.supermoons[[x]]$distance$uncertainty < 1) {
      significant.figures <- 7
   } else {
      significant.figures <- 6
   }
   cat(label, replicate(G.max.indent - nchar(label), " "), signif(asNumeric(list.supermoons[[x]]$distance$value), significant.figures), " +- ", signif(asNumeric(list.supermoons[[x]]$distance$uncertainty), 1), " km\n\n", sep="")
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
    ##  full moon:                27-Apr-2021 03:31:44 +- 00:00:36 UTC
    ##  perigee:                  27-Apr-2021 15:24:11 +- 00:31:36 UTC
    ##  distance at perigee:      357377.5 +- 1 km
    ## 
    ## supermoon #24
    ##  full moon:                26-May-2021 11:14:01 +- 00:00:45 UTC
    ##  perigee:                  26-May-2021 01:51:53 +- 00:31:00 UTC
    ##  distance at perigee:      357312 +- 1 km
    ## 
    ## supermoon #25
    ##  full moon:                14-Jun-2022 11:51:17 +- 00:00:36 UTC
    ##  perigee:                  14-Jun-2022 23:21:22 +- 00:31:36 UTC
    ##  distance at perigee:      357435 +- 2 km
    ## 
    ## supermoon #26
    ##  full moon:                13-Jul-2022 18:37:36 +- 00:00:45 UTC
    ##  perigee:                  13-Jul-2022 09:07:41 +- 00:31:18 UTC
    ##  distance at perigee:      357263.4 +- 0.8 km
    ## 
    ## supermoon #27
    ##  full moon:                01-Aug-2023 18:31:28 +- 00:00:54 UTC
    ##  perigee:                  02-Aug-2023 05:52:04 +- 00:31:18 UTC
    ##  distance at perigee:      357312 +- 1 km
    ## 
    ## supermoon #28
    ##  full moon:                31-Aug-2023 01:35:27 +- 00:00:17 UTC
    ##  perigee:                  30-Aug-2023 15:51:02 +- 00:31:00 UTC
    ##  distance at perigee:      357185 +- 4 km
    ## 
    ## supermoon #29
    ##  full moon:                18-Sep-2024 02:34:19 +- 00:00:55 UTC
    ##  perigee:                  18-Sep-2024 13:26:06 +- 00:31:19 UTC
    ##  distance at perigee:      357286 +- 0.5 km
    ## 
    ## supermoon #30
    ##  full moon:                17-Oct-2024 11:26:24 +- 00:00:55 UTC
    ##  perigee:                  17-Oct-2024 00:45:12 +- 00:31:19 UTC
    ##  distance at perigee:      357173 +- 2 km
    ## 
    ## supermoon #31
    ##  full moon:                05-Nov-2025 13:19:15 +- 00:00:55 UTC
    ##  perigee:                  05-Nov-2025 22:28:52 +- 00:31:37 UTC
    ##  distance at perigee:      356833.2 +- 0.4 km
    ## 
    ## supermoon #32
    ##  full moon:                04-Dec-2025 23:14:01 +- 00:00:38 UTC
    ##  perigee:                  04-Dec-2025 11:05:36 +- 00:31:19 UTC
    ##  distance at perigee:      356965 +- 3 km
    ## 
    ## supermoon #33
    ##  full moon:                24-Dec-2026 01:28:17 +- 00:00:47 UTC
    ##  perigee:                  24-Dec-2026 08:29:51 +- 00:31:19 UTC
    ##  distance at perigee:      356652 +- 1 km
    ## 
    ## supermoon #34
    ##  full moon:                10-Feb-2028 15:03:27 +- 00:00:47 UTC
    ##  perigee:                  10-Feb-2028 19:53:16 +- 00:31:38 UTC
    ##  distance at perigee:      356679.7 +- 0.6 km
    ## 
    ## supermoon #35
    ##  full moon:                30-Mar-2029 02:26:20 +- 00:00:39 UTC
    ##  perigee:                  30-Mar-2029 05:39:48 +- 00:31:39 UTC
    ##  distance at perigee:      356666 +- 2 km
    ## 
    ## supermoon #36
    ##  full moon:                17-May-2030 11:19:12 +- 00:00:48 UTC
    ##  perigee:                  17-May-2030 13:45:20 +- 00:31:39 UTC
    ##  distance at perigee:      357017 +- 2 km
    ## 
    ## supermoon #37
    ##  full moon:                04-Jul-2031 19:01:11 +- 00:00:57 UTC
    ##  perigee:                  04-Jul-2031 21:13:53 +- 00:31:20 UTC
    ##  distance at perigee:      357009.2 +- 0.9 km
    ## 
    ## supermoon #38
    ##  full moon:                21-Aug-2032 01:46:37 +- 00:00:57 UTC
    ##  perigee:                  21-Aug-2032 03:51:38 +- 00:31:00 UTC
    ##  distance at perigee:      356881 +- 2 km
    ## 
    ## supermoon #39
    ##  full moon:                08-Oct-2033 10:58:10 +- 00:00:49 UTC
    ##  perigee:                  08-Oct-2033 12:11:20 +- 00:31:40 UTC
    ##  distance at perigee:      356824 +- 2 km
    ## 
    ## supermoon #40
    ##  full moon:                25-Nov-2034 22:31:51 +- 00:00:49 UTC
    ##  perigee:                  25-Nov-2034 22:06:15 +- 00:31:40 UTC
    ##  distance at perigee:      356448 +- 3 km
    ## 
    ## supermoon #41
    ##  full moon:                13-Jan-2036 11:15:50 +- 00:00:50 UTC
    ##  perigee:                  13-Jan-2036 08:46:57 +- 00:31:20 UTC
    ##  distance at perigee:      356519 +- 2 km
    ## 
    ## supermoon #42
    ##  full moon:                02-Mar-2037 00:27:56 +- 00:00:59 UTC
    ##  perigee:                  01-Mar-2037 19:47:36 +- 00:31:21 UTC
    ##  distance at perigee:      356711 +- 3 km
    ## 
    ## supermoon #43
    ##  full moon:                19-Apr-2038 10:35:49 +- 00:00:59 UTC
    ##  perigee:                  19-Apr-2038 04:30:40 +- 00:31:21 UTC
    ##  distance at perigee:      356843 +- 6 km
    ## 
    ## supermoon #44
    ##  full moon:                06-Jun-2039 18:47:44 +- 00:00:51 UTC
    ##  perigee:                  06-Jun-2039 12:00:54 +- 00:31:21 UTC
    ##  distance at perigee:      357207 +- 2 km
    ## 
    ## supermoon #45
    ##  full moon:                24-Jul-2040 02:05:14 +- 00:00:51 UTC
    ##  perigee:                  23-Jul-2040 19:15:19 +- 00:31:43 UTC
    ##  distance at perigee:      357112 +- 2 km
    ## 
    ## supermoon #46
    ##  full moon:                10-Sep-2041 09:23:36 +- 00:01:00 UTC
    ##  perigee:                  10-Sep-2041 02:11:53 +- 00:31:43 UTC
    ##  distance at perigee:      357006 +- 6 km
    ## 
    ## supermoon #47
    ##  full moon:                28-Oct-2042 19:48:32 +- 00:00:43 UTC
    ##  perigee:                  28-Oct-2042 11:27:14 +- 00:31:43 UTC
    ##  distance at perigee:      356973 +- 2 km
    ## 
    ## supermoon #48
    ##  full moon:                16-Nov-2043 21:52:28 +- 00:00:53 UTC
    ##  perigee:                  17-Nov-2043 09:11:26 +- 00:31:22 UTC
    ##  distance at perigee:      356948 +- 1 km
    ## 
    ## supermoon #49
    ##  full moon:                16-Dec-2043 08:01:57 +- 00:00:53 UTC
    ##  perigee:                  15-Dec-2043 22:00:39 +- 00:31:00 UTC
    ##  distance at perigee:      356772 +- 4 km
    ## 
    ## supermoon #50
    ##  full moon:                03-Jan-2045 10:20:13 +- 00:00:17 UTC
    ##  perigee:                  03-Jan-2045 19:24:35 +- 00:31:00 UTC
    ##  distance at perigee:      356774 +- 1 km
    ## 
    ## supermoon #51
    ##  full moon:                01-Feb-2045 21:05:13 +- 00:01:02 UTC
    ##  perigee:                  01-Feb-2045 08:42:32 +- 00:31:22 UTC
    ##  distance at perigee:      357105 +- 1 km
    ## 
    ## supermoon #52
    ##  full moon:                20-Feb-2046 23:43:59 +- 00:01:02 UTC
    ##  perigee:                  21-Feb-2046 06:43:17 +- 00:31:45 UTC
    ##  distance at perigee:      356806 +- 2 km
    ## 
    ## supermoon #53
    ##  full moon:                10-Apr-2047 10:35:03 +- 00:01:03 UTC
    ##  perigee:                  10-Apr-2047 16:08:00 +- 00:31:45 UTC
    ##  distance at perigee:      356790 +- 4 km
    ## 
    ## supermoon #54
    ##  full moon:                27-May-2048 18:56:46 +- 00:00:46 UTC
    ##  perigee:                  27-May-2048 23:56:00 +- 00:31:46 UTC
    ##  distance at perigee:      357115 +- 1 km
    ## 
    ## supermoon #55
    ##  full moon:                15-Jul-2049 02:29:15 +- 00:01:04 UTC
    ##  perigee:                  15-Jul-2049 07:17:56 +- 00:31:46 UTC
    ##  distance at perigee:      357062 +- 1 km
    ## 
    ## supermoon #56
    ##  full moon:                01-Sep-2050 09:30:26 +- 00:00:47 UTC
    ##  perigee:                  01-Sep-2050 14:02:31 +- 00:31:24 UTC
    ##  distance at perigee:      356899 +- 4 km
    ## 
    ## supermoon #57
    ##  full moon:                19-Oct-2051 19:12:50 +- 00:01:06 UTC
    ##  perigee:                  19-Oct-2051 22:40:51 +- 00:31:48 UTC
    ##  distance at perigee:      356809 +- 2 km
    ## 
    ## supermoon #58
    ##  full moon:                06-Dec-2052 07:17:32 +- 00:00:58 UTC
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
    ##  perigee:                  13-Mar-2055 06:25:07 +- 00:31:52 UTC
    ##  distance at perigee:      356698 +- 4 km
    ## 
    ## supermoon #61
    ##  full moon:                29-Apr-2056 18:30:49 +- 00:01:10 UTC
    ##  perigee:                  29-Apr-2056 14:47:49 +- 00:31:53 UTC
    ##  distance at perigee:      356811 +- 5 km
    ## 
    ## supermoon #62
    ##  full moon:                17-Jun-2057 02:18:27 +- 00:01:03 UTC
    ##  perigee:                  16-Jun-2057 22:08:24 +- 00:31:27 UTC
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
    ##  full moon:                08-Nov-2060 04:17:14 +- 00:00:58 UTC
    ##  perigee:                  07-Nov-2060 22:10:49 +- 00:31:29 UTC
    ##  distance at perigee:      356812 +- 3 km
    ## 
    ## supermoon #66
    ##  full moon:                26-Dec-2061 16:52:34 +- 00:01:16 UTC
    ##  perigee:                  26-Dec-2061 08:55:32 +- 00:31:59 UTC
    ##  distance at perigee:      356619 +- 3 km
    ## 
    ## supermoon #67
    ##  full moon:                14-Jan-2063 19:11:13 +- 00:01:00 UTC
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
    ##  perigee:                  03-Mar-2064 17:31:27 +- 00:32:01 UTC
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
    ##  perigee:                  12-Sep-2068 00:17:10 +- 00:32:06 UTC
    ##  distance at perigee:      356953 +- 4 km
    ## 
    ## supermoon #75
    ##  full moon:                30-Oct-2069 03:35:07 +- 00:01:25 UTC
    ##  perigee:                  30-Oct-2069 09:14:21 +- 00:31:34 UTC
    ##  distance at perigee:      356831 +- 1 km
    ## 
    ## supermoon #76
    ##  full moon:                17-Dec-2070 16:05:22 +- 00:01:26 UTC
    ##  perigee:                  17-Dec-2070 19:40:29 +- 00:32:09 UTC
    ##  distance at perigee:      356443 +- 4 km
    ## 
    ## supermoon #77
    ##  full moon:                04-Feb-2072 04:55:16 +- 00:01:27 UTC
    ##  perigee:                  04-Feb-2072 06:25:57 +- 00:32:10 UTC
    ##  distance at perigee:      356545.9 +- 0.6 km
    ## 
    ## supermoon #78
    ##  full moon:                23-Mar-2073 17:16:51 +- 00:01:28 UTC
    ##  perigee:                  23-Mar-2073 16:57:26 +- 00:32:11 UTC
    ##  distance at perigee:      356722 +- 3 km
    ## 
    ## supermoon #79
    ##  full moon:                11-May-2074 02:17:19 +- 00:01:21 UTC
    ##  perigee:                  11-May-2074 01:01:08 +- 00:32:12 UTC
    ##  distance at perigee:      356815 +- 4 km
    ## 
    ## supermoon #80
    ##  full moon:                28-Jun-2075 09:46:14 +- 00:01:31 UTC
    ##  perigee:                  28-Jun-2075 08:12:49 +- 00:32:13 UTC
    ##  distance at perigee:      357100 +- 1 km
    ## 
    ## supermoon #81
    ##  full moon:                14-Aug-2076 17:11:36 +- 00:01:23 UTC
    ##  perigee:                  14-Aug-2076 15:29:27 +- 00:32:15 UTC
    ##  distance at perigee:      356914 +- 3 km
    ## 
    ## supermoon #82
    ##  full moon:                02-Oct-2077 01:20:21 +- 00:01:33 UTC
    ##  perigee:                  01-Oct-2077 22:58:31 +- 00:32:16 UTC
    ##  distance at perigee:      356757 +- 7 km
    ## 
    ## supermoon #83
    ##  full moon:                19-Nov-2078 12:52:16 +- 00:01:35 UTC
    ##  perigee:                  19-Nov-2078 08:56:51 +- 00:32:17 UTC
    ##  distance at perigee:      356690 +- 1 km
    ## 
    ## supermoon #84
    ##  full moon:                07-Jan-2080 01:44:33 +- 00:01:36 UTC
    ##  perigee:                  06-Jan-2080 19:50:10 +- 00:31:39 UTC
    ##  distance at perigee:      356508 +- 4 km
    ## 
    ## supermoon #85
    ##  full moon:                23-Feb-2081 14:26:57 +- 00:01:37 UTC
    ##  perigee:                  23-Feb-2081 06:17:54 +- 00:32:20 UTC
    ##  distance at perigee:      356861 +- 2 km
    ## 
    ## supermoon #86
    ##  full moon:                14-Mar-2082 16:44:37 +- 00:01:21 UTC
    ##  perigee:                  15-Mar-2082 04:16:12 +- 00:31:40 UTC
    ##  distance at perigee:      357176 +- 1 km
    ## 
    ## supermoon #87
    ##  full moon:                13-Apr-2082 01:45:08 +- 00:01:38 UTC
    ##  perigee:                  12-Apr-2082 15:53:20 +- 00:32:21 UTC
    ##  distance at perigee:      357106 +- 2 km
    ## 
    ## supermoon #88
    ##  full moon:                02-May-2083 02:29:28 +- 00:01:40 UTC
    ##  perigee:                  02-May-2083 12:57:27 +- 00:31:41 UTC
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
    ##  full moon:                17-Jul-2084 17:01:24 +- 00:01:32 UTC
    ##  perigee:                  17-Jul-2084 06:12:13 +- 00:32:24 UTC
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
    ##  full moon:                22-Oct-2086 09:55:40 +- 00:01:40 UTC
    ##  perigee:                  21-Oct-2086 22:00:07 +- 00:32:26 UTC
    ##  distance at perigee:      357178 +- 5 km
    ## 
    ## supermoon #96
    ##  full moon:                10-Nov-2087 12:04:45 +- 00:01:40 UTC
    ##  perigee:                  10-Nov-2087 19:53:15 +- 00:31:44 UTC
    ##  distance at perigee:      356889.7 +- 0.5 km
    ## 
    ## supermoon #97
    ##  full moon:                28-Dec-2088 00:57:03 +- 00:01:50 UTC
    ##  perigee:                  28-Dec-2088 06:32:34 +- 00:31:44 UTC
    ##  distance at perigee:      356502 +- 6 km
    ## 
    ## supermoon #98
    ##  full moon:                14-Feb-2090 13:39:02 +- 00:01:50 UTC
    ##  perigee:                  14-Feb-2090 17:12:55 +- 00:32:30 UTC
    ##  distance at perigee:      356620.4 +- 0.6 km
    ## 
    ## supermoon #99
    ##  full moon:                04-Apr-2091 01:31:11 +- 00:01:40 UTC
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
    ##  perigee:                  08-Jul-2093 18:17:45 +- 00:31:47 UTC
    ##  distance at perigee:      357098 +- 0.5 km
    ## 
    ## supermoon #102
    ##  full moon:                26-Aug-2094 00:51:25 +- 00:01:50 UTC
    ##  perigee:                  26-Aug-2094 01:38:31 +- 00:32:35 UTC
    ##  distance at perigee:      356867 +- 2 km
    ## 
    ## supermoon #103
    ##  full moon:                13-Oct-2095 09:30:32 +- 00:01:36 UTC
    ##  perigee:                  13-Oct-2095 09:25:59 +- 00:31:48 UTC
    ##  distance at perigee:      356687 +- 5 km
    ## 
    ## supermoon #104
    ##  full moon:                29-Nov-2096 21:33:39 +- 00:01:50 UTC
    ##  perigee:                  29-Nov-2096 19:44:15 +- 00:31:49 UTC
    ##  distance at perigee:      356609 +- 2 km
    ## 
    ## supermoon #105
    ##  full moon:                17-Jan-2098 10:35:40 +- 00:02:00 UTC
    ##  perigee:                  17-Jan-2098 06:41:26 +- 00:32:39 UTC
    ##  distance at perigee:      356437 +- 3 km
    ## 
    ## supermoon #106
    ##  full moon:                06-Mar-2099 22:59:22 +- 00:02:00 UTC
    ##  perigee:                  06-Mar-2099 16:58:42 +- 00:32:40 UTC
    ##  distance at perigee:      356797 +- 3 km
    ## 
    ## supermoon #107
    ##  full moon:                24-Apr-2100 09:43:26 +- 00:02:00 UTC
    ##  perigee:                  24-Apr-2100 02:12:39 +- 00:31:51 UTC
    ##  distance at perigee:      357012 +- 1 km
    ## 
    ## ---------------------------------------------
