# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}


# kd should be a list (or frame) with an $x, $y, and $z component.
#
# $x and $y are the x an y coordinates
#
# $z is a matrix with values coresponding to each x/y place on the grid.
#
# we will be using these $z values to find the peaks in the contours
#

#' @export
find_peaks <- function(kd, combine_within=1, top_peak_pct=0.05) {

  if ("x" %in% names(kd) && "y" %in% names(kd) && "z" %in% names(kd)) {
    # print("2d")
    return(find_peaks_2d(kd, combine_within, top_peak_pct))
  }
  if ("x" %in% names(kd) && "y" %in% names(kd)) {
    # print("1d")
    return(find_peaks_1d(kd, combine_within, top_peak_pct))
  }

  stop("Bad density object: ", kd)
}

find_peaks_1d <- function(kd, combine_within=1, top_peak_pct=0.05) {

  flow.m <- rep(0, length(kd$x))

  for (i in 1:length(kd$x)) {
      # print(paste0("i=",i,", j=",j))
      left = NA
      me = kd$y[i]
      right = NA

      if (i > 1) {
        left = kd$y[i-1]
      }
      if (i < length(kd$x)) {
        right = kd$y[i+1]
      }

      flow.m[i] <- which(c(left, me, right)==max(c(left, me, right), na.rm = T))[1]

  }


  ## Find the peaks in order of magnitude
  ##
  ## Follow the flow established above to find the peaks
  ## and then assign each position in the grid to a peak
  ##


  tmp.m <- kd$y

  # peaks will have a flow.m value of 5. Everything else is ignored for this part
  tmp.m[flow.m!=2] <- 0
  peak.m <- rep(0, length(kd$x))

  peak_nums <- c()
  peak_cols <- c()
  n <- 1

  # Identify peaks in order of peak height.
  while (max(tmp.m) > 0) {
    max_peak <- max(tmp.m)
    peak_i <- which(tmp.m == max_peak)[1]
    tmp.m[peak_i] <- 0
    peak.m[peak_i] <- n

    # look for adjacent peaks
    #
    # because we are going in order of peak height, we will merge
    # any adjacent peaks to this one.

    # print(paste0("Found peak at: [",peak_i,",",peak_j,"] ", sum(tmp.m > 0)))

    for (i in max(1, peak_i-combine_within): min(length(kd$y), peak_i+combine_within)) {
      if (tmp.m[i] > 0) {
        tmp.m[i] <- 0
        peak.m[i] <- n
      }
    }

    peak_nums <- c(peak_nums, n)
    peak_cols <- c(peak_cols, peak_i)

    n <- n + 1
  }
  # print ("1")
  ## assign each position to a peak
  while(length(which(peak.m==0))>0) {
    i <- which(peak.m==0)[1]
    peak.m[i] <- recur_peak_1d(flow.m,i,peak.m)
  }
  # print ("2")

  peak_size <- rep(0,length(peak_nums))
  for (i in 1:length(kd$x)) {
      peak <- peak.m[i]
      peak_size[peak] <- peak_size[peak] + kd$y[i]
  }
  # print ("3")

  peak_vals <- c()
  for (i in 1:length(peak_cols)) {
    peak_vals <- c(peak_vals, kd$y[peak_cols[i]])
  }
  # print ("4")

  return(list(
    num=peak_nums,
    cols=peak_cols,
    x=kd$x[peak_cols],
    size=peak_size,
    y=peak_vals,
    pct=peak_size/sum(peak_size),
    peak.m=peak.m,
    top_peaks=sum(peak_size/sum(peak_size) > top_peak_pct)
  ))
}


find_peaks_2d <- function(kd, combine_within=1, top_peak_pct=0.05) {

  flow.m <- matrix(0, nrow=nrow(kd$z), ncol=ncol(kd$z))

  for (i in 1:nrow(kd$z)) {
    for (j in 1:ncol(kd$z)) {
    # print(paste0("i=",i,", j=",j))
      ul = NA
      top = NA
      ur = NA
      left = NA
      me = kd$z[i,j]
      right = NA
      bl = NA
      bottom = NA
      br = NA

      if (i > 1) {
        if (j > 1) {
          ul = kd$z[i-1,j-1]
        }
        top = kd$z[i-1,j]
        if (j < ncol(kd$z)) {
          ur = kd$z[i-1,j+1]
        }
      }
      if (i < nrow(kd$z)) {
        if (j > 1) {
          bl = kd$z[i+1,j-1]
        }
        bottom = kd$z[i+1,j]
        if (j < ncol(kd$z)) {
          br = kd$z[i+1,j+1]
        }
      }
      if (j > 1) {
        left = kd$z[i,j-1]
      }
      if (j < ncol(kd$z)) {
        right = kd$z[i,j+1]
      }

      flow.m[i,j] <- which(c(ul, top, ur,left, me, right, bl, bottom, br)==max(c(ul, top, ur,left, me, right, bl, bottom, br), na.rm = T))[1]

    }
  }


  ## Find the peaks in order of magnitude
  ##
  ## Follow the flow established above to find the peaks
  ## and then assign each position in the grid to a peak
  ##


  tmp.m <- kd$z

  # peaks will have a flow.m value of 5. Everything else is ignored for this part
  tmp.m[flow.m!=5] <- 0
  peak.m <- matrix(0, nrow=nrow(kd$z), ncol=ncol(kd$z))

  peak_nums <- c()
  peak_rows <- c()
  peak_cols <- c()
  n <- 1

  # Identify peaks in order of peak height.
  while (max(tmp.m) > 0) {
    max_peak <- max(tmp.m)
    peak_xy <- which(tmp.m == max_peak, arr.ind = TRUE)
    #print(peak_xy)
    peak_i <- peak_xy[1,1]
    peak_j <- peak_xy[1,2]
    tmp.m[peak_i, peak_j] <- 0
    peak.m[peak_i, peak_j] <- n

    # if (length(peak_rows) > 0 && peak_i == peak_rows[length(peak_rows)] && peak_j == peak_cols[length(peak_cols)] ) {
    #   print(tmp.m)
    # }

    # look for adjacent peaks
    #
    # because we are going in order of peak height, we will merge
    # any adjacent peaks to this one.

    # print(paste0("Found peak at: [",peak_i,",",peak_j,"] ", sum(tmp.m > 0), " max: ", max_peak, ", tmp.m[]=",tmp.m[peak_i, peak_j]))

    for (i in max(1, peak_i-combine_within): min(nrow(kd$z), peak_i+combine_within)) {
      for (j in max(1, peak_j-combine_within): min(ncol(kd$z), peak_j+combine_within)) {
        if (tmp.m[i,j] > 0) {
          # print(paste0("Combining peak[",i,",",j,"] with existing peak [",peak_i,",",peak_j,"]"))
          tmp.m[i,j] <- 0
          peak.m[i,j] <- n
        }
      }
    }

    peak_nums <- c(peak_nums, n)
    peak_rows <- c(peak_rows, peak_i)
    peak_cols <- c(peak_cols, peak_j)

    n <- n + 1
  }

  # print("1")

  ## assign each position to a peak
  while(length(which(peak.m==0, arr.ind = T))>0) {
    i <- which(peak.m==0, arr.ind = T)[1,1]
    j <- which(peak.m==0, arr.ind = T)[1,2]
    peak.m[i,j] <- recur_peak_2d(flow.m,i,j,peak.m)
  }
  # print("2")

  peak_size <- rep(0,length(peak_nums))
  for (i in 1:nrow(kd$z)) {
    for (j in 1:ncol(kd$z)) {
      peak <- peak.m[i,j]
      peak_size[peak] <- peak_size[peak] + kd$z[i,j]
    }
  }
  # print("3")

  peak_vals <- c()
  for (i in 1:length(peak_rows)) {
    peak_vals <- c(peak_vals, kd$z[peak_rows[i], peak_cols[i]])
  }
  # print("4")

  return(list(
    num=peak_nums,
    rows=peak_rows,
    cols=peak_cols,
    x=kd$x[peak_rows],
    y=kd$y[peak_cols],
    size=peak_size,
    z=peak_vals,
    pct=peak_size/sum(peak_size),
    peak.m=peak.m,
    top_peaks=sum(peak_size/sum(peak_size) > top_peak_pct)
  ))
}

recur_peak_2d <- function(mat, row, col, ret.m) {
  if (ret.m[row,col] != 0) {
    return (ret.m[row,col])
  }
  if (mat[row,col] == 1) {
    v <- recur_peak_2d(mat, row-1, col-1, ret.m)
    return (v)
  }
  if (mat[row,col] == 2) {
    v <- recur_peak_2d(mat, row-1, col, ret.m)
    return (v)
  }
  if (mat[row,col] == 3) {
    v <- recur_peak_2d(mat, row-1, col+1, ret.m)
    return (v)
  }
  if (mat[row,col] == 4) {
    v <- recur_peak_2d(mat, row, col-1, ret.m)
    return (v)
  }
  if (mat[row,col] == 5) {
    return (ret.m[row,col])
  }
  if (mat[row,col] == 6) {
    v <- recur_peak_2d(mat, row, col+1, ret.m)
    return (v)
  }
  if (mat[row,col] == 7) {
    v <- recur_peak_2d(mat, row+1, col-1, ret.m)
    return (v)
  }
  if (mat[row,col] == 8) {
    v <- recur_peak_2d(mat, row+1, col, ret.m)
    return (v)
  }
  if (mat[row,col] == 9) {
    v <- recur_peak_2d(mat, row+1, col+1, ret.m)
    return (v)
  }
}



recur_peak_1d <- function(v, i, ret.m) {
  if (ret.m[i] != 0) {
    return (ret.m[i])
  }
  if (v[i] == 1) {
    v <- recur_peak_1d(v, i-1, ret.m)
    return (v)
  }
  if (v[i] == 2) {
    return(ret.m[i])
  }
  if (v[i] == 3) {
    v <- recur_peak_1d(v, i+1, ret.m)
    return (v)
  }
}

