# This R script downloads selected files comprising the third-party lunar
# solution ELP 2000-82B (Éphéméride Lunaire Parisienne 2000-82B).
# The files are downloaded from Laboratoire Temps Espace
# (formerly Institut de Mécanique Céleste et de Calcul des Éphémérides).
#
# The R function download.file is called to download files.
# R's default download timeout is 60 seconds if network connectivity is lost.

multiplatform_download <- function(target_url, destination_file) {

  # Download a remote file to the local file system.
  # The function is designed to be as multiplatform-compatible as possible in R.
  #
  # Input:
  #   target_url        String. A target URL to download from.
  #   destination_file  String. A name to give to the downloaded file.
  #
  # Output:
  #   <unnamed>         Integer. A code specifying download status.
  #                        Zero if successful, non-zero otherwise.

  tryCatch({
    download.file(
      target_url,
      destination_file,
      method = "auto",
      mode = "w"
    )
  },
  error = function(e) {
    -1  # Download failure status code.
  }
  )
}

process_download_status <- function(input_result, num_downloads) {

  # Test if a download operation return value indicates success,
  # and return the count of successful of downloads.
  #
  # Input:
  #   input_result   Integer. Code specifying download operation status.
  #   num_downloads  Integer. Old number of successful downloads.
  #
  # Output:
  #   num_downloads  Integer. New number of successful downloads.

  # The multiplatform_download function returns 0 if successful.
  if (input_result == 0) {
    num_downloads + 1
  } else {
    num_downloads
  }
}

# Initialise download tracking and configuration variables.
files_target <- 37  # Number of files that need to be downloaded.
files_downloaded <- 0 # Number of files actually downloaded.
#
url_prefix <- "https://ftp.imcce.fr/pub/ephem/moon/elp82b/"
local_path_prefix <- paste0(this.path::this.dir(), "/")
# Avoiding underscore use in file names for portability.
fortran_file <- "elp82b_2.f"
cat("Downloading files...", fill = TRUE)

# Download the Fortran program.
# N.B. The program is legacy, fixed-form Fortran, not modern, free-form Fortran.
download_result <- multiplatform_download(
  target_url = paste0(url_prefix, "elp82b_2"),
  destination_file = paste0(local_path_prefix, fortran_file)
)
files_downloaded <- process_download_status(download_result, files_downloaded)

# Lint the Fortran program.
if (files_downloaded == 1) {

  # For portability, use a lowercase, underscore-free subroutine name.
  # Assume two instances of the name need modification.
  # N.B. Using "fixed = TRUE" matches a literal string, not a regex.
  updated_code <- sub(
    "ELP82B_2", "elp82b2", readLines(fortran_file), fixed = TRUE
  )
  updated_code <- sub("ELP82B_2", "elp82b2", updated_code, fixed = TRUE)

  # Annotate the file to leave a linting record.
  updated_code <- sub(
    "2000-82B.",
    paste0(
      "2000-82B.\n*\n",
      "*     N.B. This file is a linted derivative of the original.\n",
      "*          Linting date: ", format(Sys.Date(), "%d-%b-%Y")
    ),
    updated_code, fixed = TRUE
  )

  # Explicitly declare variable types.
  # Ensure Fortran uses 4-byte integers to match R's 4-byte integers.
  updated_code <- sub(
    "fich*60",
    "fich*60\n      integer*4 i,ideb,ierr,ific,ir,itab,iv,iz,j,k,nt,nul",
    updated_code, fixed = TRUE
  )
  updated_code <- sub("ilu(4),ipla(11)", "zone(6)", updated_code, fixed = TRUE)
  updated_code <- sub(
    "dimension nterm(3,12),nrang(3,12),zone(6)",
    "integer*4 nterm(3,12),nrang(3,12),ilu(4),ipla(11)",
    updated_code, fixed = TRUE
  )

  # Remove an unused array element from the program.
  updated_code <- sub(",t(0:4)", ",t(1:4)", updated_code, fixed = TRUE)
  updated_code <- sub("t/1.d0,4*0.d0/", "t/4*0.d0/", updated_code, fixed = TRUE)

  # Remove an unused Fortran "parameter" (a constant) from the program.
  updated_code <- sub("cpi2=2.d0*cpi,", "", updated_code, fixed = TRUE)

  # Store large arrays on the heap, not the stack.
  updated_code <- sub(
    "dimension per1(3,19537),per2(3,6766),per3(3,8924)",
    "double precision, allocatable :: per1(:,:),per2(:,:),per3(:,:)",
    updated_code, fixed = TRUE
  )
  updated_code <- sub(
    "4*0.d0/",
    "4*0.d0/\n      allocate (per1(3,19537),per2(3,6766),per3(3,8924))",
    updated_code, fixed = TRUE
  )
  # Ensure large arrays are removed from the heap after use.
  updated_code <- sub(
    "*     Change",
    "      deallocate(per1, per2, per3)\n*\n*     Change",
    updated_code, fixed = TRUE
  )

  # Ensure variables are initialised before they're read.
  updated_code <- sub(
    "t(3)*t(1)",
    paste0(
      "t(3)*t(1)\n      per1(:,:)=0\n      per2(:,:)=0\n      per3(:,:)=0\n",
      "      nrang(:,:)=0"
    ),
    updated_code, fixed = TRUE
  )

  # Remove multi-session precision data caching.
  updated_code <- sub("prec0/-1.d0/,", "", updated_code, fixed = TRUE)
  updated_code <- sub(
    "nrang(:,:)=0",
    "nrang(:,:)=0\n      prec0=-1.d0",
    updated_code, fixed = TRUE
  )

  # Compare "REAL(8)" (double-precision floating point) values reliably.
  # Assume that if the absolute value of their difference is greater than the
  # machine epsilon for "REAL(8)" values, they're not equal.
  updated_code <- sub(
    "prec.ne.prec0", "abs(prec - prec0) > epsilon(1.0_8)", updated_code,
    fixed = TRUE
  )

  writeLines(updated_code, fortran_file)
  rm(updated_code)
}

# Download the associated data files.
for (i in 1:36) {
  download_result <- multiplatform_download(
    target_url = paste0(url_prefix, "ELP", i),
    destination_file = paste0(local_path_prefix, "ELP", i)
  )
  files_downloaded <- process_download_status(download_result, files_downloaded)
}

if (files_downloaded == files_target) {
  cat("Download complete.", fill = TRUE)
} else {
  cat(
    "Download finished incomplete. ",
    files_downloaded, "/", files_target, " files downloaded.",
    sep = "",
    fill = TRUE
  )
}
