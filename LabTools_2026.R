# =============================================================================
# JLC Lab Tools
# =============================================================================
# General-purpose utilities for data import, processing, plotting, and
# quantification used across the JLC/Howorka lab.
#
# Conventions used throughout:
#   - Data frames passed as the primary argument are always named `input_df`
#   - Boolean-like parameters use TRUE/FALSE or NULL (not "y"/"n" strings)
#   - File paths constructed with file.path() where possible
#   - message() used for progress; warning() for recoverable issues
# =============================================================================


# ── Section 0: Package Setup ──────────────────────────────────────────────────

.jlc_pkgs <- c(
  "ggplot2", "grid", "gridExtra", "cowplot", "zoo",
  "readxl", "scales", "openxlsx", "ggthemes", "ggside",
  "patchwork", "drc"
)
.not_installed <- .jlc_pkgs[!(.jlc_pkgs %in% installed.packages()[, "Package"])]
if (length(.not_installed)) install.packages(.not_installed)
invisible(lapply(.jlc_pkgs, require, character.only = TRUE))


# ── Section 1: Global Theme Elements ─────────────────────────────────────────

#' ggplot2 theme modifier: remove all facet strip labels
#'
#' Add to any ggplot to strip facet labels while keeping panel borders.
#' Usage: `my_plot + strip_blank`
lab_rmv <- ggplot2::theme(
  strip.background   = ggplot2::element_blank(),
  strip.text.x       = ggplot2::element_blank(),
  strip.text.y       = ggplot2::element_blank(),
  strip.background.x = ggplot2::element_blank(),
  strip.background.y = ggplot2::element_blank()
)

#' ggplot2 theme modifier: clean up ggside margin panels
#'
#' Removes gridlines, backgrounds, borders, and axis decorations from
#' ggside marginal panels. Add to plots that use geom_xsidedensity etc.
#' Usage: `my_plot + side_panel_clean`
side_panel_clean <- ggplot2::theme(
  ggside.panel.grid       = ggplot2::element_blank(),
  ggside.panel.background = ggplot2::element_blank(),
  ggside.panel.border     = ggplot2::element_blank(),
  ggside.axis.text        = ggplot2::element_blank(),
  ggside.axis.ticks       = ggplot2::element_blank()
)

# Global font size defaults — change these to rescale all plots uniformly
.txt_size <- 12
.txt_fam  <- "sans"


#' Lab ggplot2 theme (bw base)
#'
#' Clean black-and-white theme with consistent text sizing. Intended to be
#' called as a function so txt_size can be overridden per-plot.
#'
#' @param base_size  numeric  Base font size passed to theme_bw (default: 12).
#' @param txt_size   numeric  Text size for all labels (default: 12).
#' @param base_family character  Font family (default: "sans").
#'
#' @return ggplot2 theme object.
My_Theme <- function(base_size = 12, txt_size = 12, base_family = "sans") {
  
  txt      <- ggplot2::element_text(size = txt_size, colour = "black", face = "plain")
  bold_txt <- ggplot2::element_text(size = txt_size, colour = "black", face = "bold")
  
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      strip.background    = ggplot2::element_rect(fill = "white", color = "black"),
      plot.title          = ggplot2::element_text(size = txt_size, family = base_family),
      plot.margin         = ggplot2::unit(c(0.6, 0.6, 0.6, 0.6), "cm"),
      panel.background    = ggplot2::element_rect(fill = "white", colour = "white"),
      panel.spacing       = ggplot2::unit(2, "lines"),
      text                = txt,
      axis.title.x        = txt,
      axis.title.y        = txt,
      axis.text.x         = txt,
      axis.text.y         = txt,
      strip.text          = txt,
      strip.text.x        = txt,
      strip.text.y        = txt,
      legend.title        = bold_txt,
      legend.text         = txt
    )
}


#' Lab ggplot2 theme (ggthemes::theme_few base)
#'
#' Alternative minimal theme. Otherwise identical signature to My_Theme.
#'
#' @param base_size   numeric  Base font size (default: 12).
#' @param txt_size    numeric  Text size for all labels (default: 12).
#' @param base_family character  Font family (default: "sans").
#'
#' @return ggplot2 theme object.
stefTheme <- function(base_size = 12, txt_size = 12, base_family = "sans") {
  
  txt      <- ggplot2::element_text(size = txt_size, colour = "black", face = "plain")
  bold_txt <- ggplot2::element_text(size = txt_size, colour = "black", face = "bold")
  
  ggthemes::theme_few(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", color = "black"),
      plot.title       = ggplot2::element_text(size = txt_size, family = base_family),
      plot.margin      = ggplot2::unit(c(0.6, 0.6, 0.6, 0.6), "cm"),
      panel.background = ggplot2::element_rect(fill = "white", colour = "white"),
      panel.spacing    = ggplot2::unit(2, "lines"),
      text             = txt,
      axis.text.x      = txt,
      axis.text.y      = txt,
      axis.title.x     = txt,
      axis.title.y     = txt,
      strip.text       = txt,
      legend.title     = bold_txt,
      legend.text      = txt
    )
}


# ── Section 2: Utility Functions ─────────────────────────────────────────────

#' Min-max normalise a numeric vector to [0, 1]
#'
#' @param x  numeric  Vector to normalise.
#'
#' @return numeric vector scaled to [0, 1]. Returns NaN for constant input.
norm <- function(x) {
  lo <- min(x, na.rm = TRUE)
  hi <- max(x, na.rm = TRUE)
  (x - lo) / (hi - lo)
}


#' Row-wise mean and SD across a list of equal-length columns
#'
#' @param col_list  list  Columns of the same length to average row-wise.
#' @param id        character  Label appended as the `id` column.
#'
#' @return data.frame with columns: y (row mean), stdev (row SD), id.
avg <- function(col_list, id) {
  df        <- data.frame(col_list)
  out       <- data.frame(
    y     = rowMeans(df, na.rm = TRUE),
    stdev = apply(df, 1, sd, na.rm = TRUE),
    id    = id
  )
  return(out)
}


#' Rolling mean and SD on the second column of a data frame
#'
#' Computes a block rolling mean and SD, then decimates to one row per block
#' so the output has ~1/block as many rows as the input.
#'
#' @param input_df  data.frame  Two-column (x, y) input.
#' @param block     integer  Rolling window size (default: 5).
#' @param id        character  Label appended as the `id` column.
#'
#' @return data.frame with columns: x, y (rolling mean), stdev (rolling SD), id.
run_avg <- function(input_df, block = 5, id = "id_missing") {
  
  # Core computation (unchanged)
  .compute <- function(df, id_val) {
    rm_vals   <- zoo::rollmean(df[, 2], k = block, na.pad = TRUE)
    sd_vals   <- zoo::rollapply(df[, 2], width = block, FUN = sd, fill = NA)
    out       <- data.frame(df[, 1], rm_vals, sd_vals)
    out       <- out[seq(2, nrow(out), block), ]
    names(out) <- c("x", "y", "stdev")
    out$id    <- id_val
    return(out)
  }
  
  # Auto-detect id column (looks for a column literally named "id")
  if ("id" %in% tolower(names(input_df))) {
    id_col  <- which(tolower(names(input_df)) == "id")
    id_vals <- unique(input_df[, id_col])
    
    results <- lapply(id_vals, function(i) {
      subset_df <- input_df[input_df[, id_col] == i, -id_col]
      .compute(subset_df, id_val = i)
    })
    
    return(do.call(rbind, results))
    
  } else {
    return(.compute(input_df, id_val = id))
  }
}


#' Add a trailing slash to a directory string if absent
#'
#' @param path  character  Directory path.
#'
#' @return character  Path guaranteed to end with "/".
endSlash <- function(path) {
  if (!endsWith(path, "/")) path <- paste0(path, "/")
  return(path)
}


#' Row-bind two data frames with different columns (fills missing with NA)
#'
#' @param x  data.frame
#' @param y  data.frame
#'
#' @return data.frame combining all columns from both, with NA where absent.
rbind_fill <- function(x, y) {
  all_cols <- union(names(x), names(y))
  x[setdiff(all_cols, names(x))] <- NA
  y[setdiff(all_cols, names(y))] <- NA
  rbind(x[all_cols], y[all_cols])
}


#' Generate axis breaks and labels in Howorka-group publication style
#'
#' Produces 5 breaks (min, Q1, mid, Q3, max) for linear axes, or decade
#' breaks for log axes. Every other label is blanked to reduce clutter.
#'
#' @param axis_limits  numeric vector  c(min, max) axis limits.
#' @param log          character/NULL  Non-NULL activates log-decade spacing.
#'
#' @return data.frame with columns: breaks, labels.
stefBreaks <- function(axis_limits, log = NULL) {
  
  lo <- axis_limits[1]
  hi <- axis_limits[2]
  
  if (!is.null(log)) {
    # Log-decade breaks between lo and hi
    if (hi > lo) {
      n_steps <- ceiling(log10(hi / lo))
      breaks  <- lo * 10^(seq(0, n_steps))
    } else {
      n_steps <- ceiling(log10(lo / hi))
      breaks  <- hi * 10^(seq(0, n_steps))
    }
  } else {
    span   <- hi - lo
    mid    <- lo + span / 2
    breaks <- c(lo, lo + span / 4, mid, hi - span / 4, hi)
  }
  
  labels <- as.character(breaks)
  labels[seq(2, length(labels), 2)] <- ""
  
  return(data.frame(breaks = breaks, labels = labels))
}


#' Convenience wrapper: apply stefBreaks to both axes simultaneously
#'
#' Returns a list of two ggplot2 scale calls ready to be added to a plot.
#'
#' @param xlim  numeric vector  c(xmin, xmax).
#' @param ylim  numeric vector  c(ymin, ymax).
#'
#' @return list of two ggplot2 scale_*_continuous() calls.
stefXY <- function(xlim, ylim) {
  list(
    ggplot2::scale_x_continuous(
      breaks = stefBreaks(xlim)[[1]],
      labels = stefBreaks(xlim)[[2]],
      limits = xlim
    ),
    ggplot2::scale_y_continuous(
      breaks = stefBreaks(ylim)[[1]],
      labels = stefBreaks(ylim)[[2]],
      limits = ylim
    )
  )
}


#' Round up to the nearest power of ten
#'
#' @param x  numeric  Value to round up.
#'
#' @return numeric  Smallest power of 10 ≥ x.
roundUp <- function(x) 10^ceiling(log10(x))


#' Generate a colour ramp between two colours
#'
#' Thin wrapper around colorRampPalette for convenience.
#'
#' @param col1    character  Start colour (name or hex).
#' @param col2    character  End colour (name or hex).
#' @param n_cols  integer  Number of colours to generate.
#'
#' @return character vector of `n_cols` hex colour codes.
colRampR <- function(col1, col2, n_cols) {
  colorRampPalette(c(col1, col2))(n_cols)
}


#' Continuous modal average via kernel density estimation
#'
#' Returns the x-value at the peak of the density estimate. Falls back to
#' the median if density estimation fails or there are too few values.
#'
#' @param x  numeric  Input vector.
#'
#' @return numeric scalar: estimated mode.
modalAvg <- function(x) {
  
  x <- stats::na.omit(x)
  
  if (length(x) < 2) {
    warning("modalAvg: fewer than 2 values — returning median.")
    return(stats::median(x, na.rm = TRUE))
  }
  
  d <- tryCatch(
    stats::density(x, bw = "nrd0"),
    error = function(e) {
      warning("modalAvg: density estimation failed — returning median.")
      NULL
    }
  )
  
  if (is.null(d)) stats::median(x, na.rm = TRUE)
  else            d$x[which.max(d$y)]
}


#' Add log-scale tick marks to a faceted ggplot
#'
#' Created by neuron & mne1 on StackOverflow. Useful because the standard
#' annotation_logticks() cannot be targeted at specific facets.
#'
#' @param base       numeric  Log base (default: 10).
#' @param sides      character  Which axes: "b" bottom, "l" left, "bl" both.
#' @param scaled     logical  Whether to scale ticks with the axis.
#' @param tick_size  numeric  Relative tick length multiplier (default: 1).
#' @param colour     character  Tick colour (default: "black").
#' @param linetype   integer  Line type (default: 1).
#' @param alpha      numeric  Opacity (default: 1).
#' @param data       data.frame  Used to target specific facets, e.g.
#'                    data.frame(x = NA, id = "my_facet_label").
#' @param ...        Additional arguments forwarded to the layer.
#'
#' @return ggplot2 layer.
#'
#' @examples
#' # Add log ticks to only the facet labelled "ctrl":
#' my_plot + add_logticks(sides = "bl", data = data.frame(x = NA, id = "ctrl"))
add_logticks <- function(base      = 10,
                         sides     = "bl",
                         scaled    = TRUE,
                         tick_size = 1,
                         colour    = "black",
                         linetype  = 1,
                         alpha     = 1,
                         data      = data.frame(x = NA),
                         ...) {
  short <- ggplot2::unit(0.1 * tick_size, "cm")
  mid   <- ggplot2::unit(0.2 * tick_size, "cm")
  long  <- ggplot2::unit(0.3 * tick_size, "cm")
  
  ggplot2::layer(
    geom     = "logticks",
    params   = list(base = base, sides = sides, scaled = scaled,
                    short = short, mid = mid, long = long,
                    colour = colour, linetype = linetype, alpha = alpha,
                    outside = TRUE, ...),
    stat     = "identity",
    data     = data,
    mapping  = NULL,
    inherit.aes  = FALSE,
    position = "identity",
    show.legend  = FALSE
  )
}


# ── Section 3: Data Import ────────────────────────────────────────────────────

#' Import a single two-column (x, y) CSV or text file
#'
#' Reads the file, selects the requested columns, coerces to numeric, removes
#' NA rows, and optionally normalises or zero-subtracts the axes.
#'
#' @param file_name  character  Filename (extension optional if supplied in `ext`).
#' @param direc      character  Directory containing the file.
#' @param ext        character  File extension including the dot (default: ".csv").
#' @param normalise  character/NULL  Normalisation mode:
#'                     "x"    → subtract first x value (zero x axis),
#'                     "zero" → subtract minimum y value,
#'                     "yes"  → 0-1 normalise y via norm(),
#'                     NULL   → no normalisation.
#' @param precut     numeric  Discard rows where x ≤ this value (default: 0).
#' @param columns    integer vector  Column indices to import; must be an even
#'                     number (pairs of x, y). Default: c(1, 2).
#'
#' @return data.frame with columns: x, y (repeated for multi-column imports),
#'   id (= file_name).
xy_data_importr <- function(file_name,
                            direc,
                            ext       = ".csv",
                            normalise = "x",
                            precut    = 0,
                            columns   = c(1, 2)) {
  
  # Strip extension from file_name if it already contains it
  if (grepl(ext, file_name, fixed = TRUE)) ext <- ""
  
  dat        <- read.csv(file.path(direc, paste0(file_name, ext)), header = FALSE)
  dat        <- dat[, columns, drop = FALSE]
  
  # Name columns as alternating x/y pairs
  names(dat) <- rep(c("x", "y"), length(columns) / 2)
  
  # Coerce to numeric and drop incomplete rows
  dat[]  <- lapply(dat, function(col) as.numeric(as.character(col)))
  dat    <- stats::na.omit(dat)
  dat    <- subset(dat, dat$x > precut)
  
  # Normalisation
  if (!is.null(normalise)) {
    nrm <- tolower(normalise)
    if (nrm == "x") {
      dat$x <- dat$x - dat$x[1]
    } else if (nrm %in% c("zero", "0")) {
      dat$y <- dat$y - min(dat$y, na.rm = TRUE)
    } else if (nrm %in% c("yes", "y")) {
      dat$y <- norm(dat$y)
    }
  }
  
  dat$id <- file_name
  dat    <- stats::na.omit(dat)
  return(dat)
}


#' Import multiple two-column CSV files into one long-format data frame
#'
#' Wrapper around xy_data_importr that reads a list of files and row-binds
#' the results.
#'
#' @param file_list  character vector  Filenames to import.
#' @param direc      character  Directory containing all files.
#' @param ext        character  File extension (default: ".csv").
#' @param normalise  character/NULL  Passed to xy_data_importr (default: NULL).
#' @param precut     numeric  Passed to xy_data_importr (default: 0).
#' @param columns    integer vector  Passed to xy_data_importr (default: c(1,2)).
#'
#' @return data.frame: row-bound output of xy_data_importr for all files.
data_framr <- function(file_list,
                       direc,
                       ext       = ".csv",
                       normalise = NULL,
                       precut    = 0,
                       columns   = c(1, 2)) {
  
  # If any file_name already contains the extension, disable adding it
  ext <- ifelse(any(grepl(ext, file_list, fixed = TRUE)), "", ext)
  if (length(ext) > 1) ext <- ext[1]
  
  do.call(rbind, lapply(
    file_list, xy_data_importr,
    direc     = direc,
    ext       = ext,
    normalise = normalise,
    precut    = precut,
    columns   = columns
  ))
}


#' Import multi-cuvette fluorimeter CSV output to long format
#'
#' Reads a CSV where each cuvette contributes two adjacent columns (x, y).
#' Handles the trailing NA columns some fluorimeters append. Converts to
#' long format with one row per (x, cuvette) observation.
#'
#' @param file_name  character  Filename.
#' @param direc      character  Directory containing the file.
#' @param ext        character  File extension (default: ".csv").
#' @param normalise  character/NULL  Normalisation mode:
#'                     "yes"  → 0-1 normalise × 100,
#'                     "zero" → subtract first y value,
#'                     "min"  → subtract minimum y value,
#'                     NULL   → no normalisation.
#'
#' @return data.frame with columns: x, y, sample (cuvette label), id (file_name).
cuvette_importer <- function(file_name,
                             direc,
                             ext       = ".csv",
                             normalise = NULL) {
  
  if (grepl(ext, file_name, fixed = TRUE)) ext <- ""
  
  raw <- try(read.csv(file.path(direc, paste0(file_name, ext)), header = TRUE))
  if (inherits(raw, "try-error"))
    stop("cuvette_importer: could not read '", file_name, "'")
  
  # Drop the header repeat row that some fluorimeters insert
  raw <- utils::tail(raw, -1)
  # Drop any column that is entirely NA
  raw <- raw[, !apply(is.na(raw), 2, all), drop = FALSE]
  # Coerce to numeric and drop incomplete rows
  raw[] <- lapply(raw, function(col) as.numeric(as.character(col)))
  raw   <- raw[, colSums(is.na(raw)) < nrow(raw), drop = FALSE]
  raw   <- raw[stats::complete.cases(raw), ]
  
  # Every even-positioned column is a y; its preceding column is the x.
  # Column names are the cuvette/sample labels.
  sample_names <- names(raw)[seq(1, ncol(raw), 2)]
  
  out_df <- do.call(rbind, lapply(seq_along(sample_names), function(i) {
    xi  <- (i - 1) * 2 + 1
    yi  <- xi + 1
    tmp <- data.frame(x = raw[, xi], y = raw[, yi], sample = sample_names[i])
    
    if (!is.null(normalise)) {
      nrm <- tolower(normalise)
      if (nrm %in% c("y", "yes")) {
        tmp$y <- norm(tmp$y) * 100
      } else if (nrm %in% c("zero", "0")) {
        tmp$y <- tmp$y - tmp$y[1]
      } else if (nrm == "min") {
        tmp$y <- tmp$y - min(tmp$y, na.rm = TRUE)
      }
    }
    tmp
  }))
  
  out_df    <- stats::na.omit(out_df)
  out_df$id <- file_name
  return(out_df)
}


#' Import multiple multi-cuvette files into one long-format data frame
#'
#' Wrapper around cuvette_importer.
#'
#' @param file_list  character vector  Filenames to import.
#' @param direc      character  Directory containing all files.
#' @param ext        character  File extension (default: ".csv").
#' @param normalise  character/NULL  Passed to cuvette_importer.
#'
#' @return data.frame: row-bound output of cuvette_importer for all files.
cuvette_framer <- function(file_list,
                           direc,
                           ext       = ".csv",
                           normalise = NULL) {
  do.call(rbind, lapply(
    file_list, cuvette_importer,
    direc     = direc,
    ext       = ext,
    normalise = normalise
  ))
}


#' Read a 96-well plate fluorescence file (xlsx) into a tidy data frame
#'
#' Subsets the raw Excel output to the data region and attaches well labels
#' (A1, B2, …). Designed for use with mapply over a list of plate files.
#'
#' @param data_path   character  Full path to the xlsx file.
#' @param data_region integer vector  c(row_start, row_end, col_start, col_end)
#'                      giving the rectangular data region within the sheet
#'                      (default: c(44, 51, 2, 13)).
#' @param id          character  Plate/experiment identifier.
#'
#' @return data.frame with columns: well_pos (e.g. "A1"), fluor_int, id.
#'
#' @examples
#' # Read multiple plates:
#' df <- do.call(rbind, mapply(
#'   read_96well_plate,
#'   data_path = paste0(data_dir, plate_files),
#'   id        = plate_ids,
#'   SIMPLIFY  = FALSE   # essential
#' ))
read_96well_plate <- function(data_path,
                              data_region = c(44, 51, 2, 13),
                              id) {
  
  raw  <- readxl::read_xlsx(data_path)
  raw  <- raw[data_region[1]:data_region[2], data_region[3]:data_region[4]]
  
  # Build standard well labels A1–H12
  well_grid <- expand.grid(
    row = LETTERS[1:8],
    col = 1:12,
    stringsAsFactors = FALSE
  )
  well_grid$label <- paste0(well_grid$row, well_grid$col)
  
  out <- data.frame(
    well_pos  = well_grid$label,
    fluor_int = as.numeric(unlist(raw)),
    id        = as.character(id),
    stringsAsFactors = FALSE
  )
  return(out)
}


# ── Section 4: Plotting ───────────────────────────────────────────────────────
#' 
#' #' General-purpose xy scatter / line plot
#' #'
#' #' Builds a ggplot2 line or point plot from a data frame. Supports error
#' #' ribbons or bars, manual colours, faceting, log axes, and axis limits.
#' #'
#' #' @param input_df    data.frame  Must contain columns named by `xvar` and
#' #'                     `yvar`. For error display it also needs a `stdev` column.
#' #' @param xvar        character  Column name for the x axis (default: "x").
#' #' @param yvar        character  Column name for the y axis (default: "y").
#' #' @param xlab        character/NULL  X axis label. NULL defaults to `xvar`.
#' #' @param ylab        character/NULL  Y axis label. NULL defaults to `yvar`.
#' #' @param colour_by   character  Column used for colour/fill grouping
#' #'                     (default: "id").
#' #' @param leg_loc     character  Legend position: "right", "left", "top",
#' #'                     "bottom", or "none" (default: "right").
#' #' @param error       character  Error display: "full" → full ±SD,
#' #'                     "half" → ±0.5 SD, "none" → no error (default: "none").
#' #' @param error_type  character  "ribbon" (default) or "errorbar".
#' #' @param cols        character/NULL  Manual colour values (recycled if short).
#' #' @param facet       character/NULL  Column name to facet by (rows ~ facet).
#' #' @param marker      character  "line" (default) or "point".
#' #' @param marker_size numeric  Geom size (default: 1).
#' #' @param n_ticks     integer  Approximate number of x-axis breaks (default: 10).
#' #' @param log         character  Log-scale axes: "x", "y", or "none" (default).
#' #' @param sci_notation logical  TRUE enables scientific notation on axes
#' #'                     (default: FALSE — plain number format).
#' #' @param leg_title   character/NULL  Legend title. NULL removes the title.
#' #' @param plot_ylim   numeric/NULL  Y axis limits c(min, max).
#' #' @param plot_xlim   numeric/NULL  X axis limits c(min, max).
#' #' @param txt_size    numeric  Base font size passed to My_Theme (default: 12).
#' #'
#' #' @return ggplot2 object.
#' plotr <- function(input_df,
#'                   xvar         = "x",
#'                   yvar         = "y",
#'                   xlab         = NULL,
#'                   ylab         = NULL,
#'                   colour_by    = "id",
#'                   leg_loc      = "right",
#'                   error        = "none",
#'                   error_type   = "ribbon",
#'                   cols         = NULL,
#'                   facet        = NULL,
#'                   marker       = "line",
#'                   marker_size  = 1,
#'                   n_ticks      = 10,
#'                   log          = "none",
#'                   sci_notation = FALSE,
#'                   leg_title    = NULL,
#'                   plot_ylim    = NULL,
#'                   plot_xlim    = NULL,
#'                   txt_size     = 12) {
#'   
#'   if (is.null(xlab)) xlab <- xvar
#'   if (is.null(ylab)) ylab <- yvar
#'   
#'   base_aes <- ggplot2::aes_string(
#'     x      = xvar,
#'     y      = yvar,
#'     colour = colour_by,
#'     fill   = colour_by
#'   )
#'   
#'   # ── Error display ───────────────────────────────────────────────────────────
#'   err_mode <- tolower(error)
#'   ribbon_layer <- NULL
#'   
#'   if (err_mode %in% c("full", "y", "yes")) {
#'     input_df$ymax <- input_df[[yvar]] + input_df$stdev
#'     input_df$ymin <- input_df[[yvar]] - input_df$stdev
#'   } else if (err_mode %in% c("half", "0.5")) {
#'     input_df$ymax <- input_df[[yvar]] + 0.5 * input_df$stdev
#'     input_df$ymin <- input_df[[yvar]] - 0.5 * input_df$stdev
#'   }
#'   
#'   if (err_mode != "none" && err_mode != "n") {
#'     err_aes <- ggplot2::aes_string(ymax = "ymax", ymin = "ymin")
#'     ribbon_layer <- if (tolower(error_type) == "ribbon") {
#'       ggplot2::geom_ribbon(err_aes, alpha = 0.25, colour = NA)
#'     } else {
#'       ggplot2::geom_errorbar(err_aes, alpha = 0.5, colour = "black",
#'                              width = marker_size/5)
#'     }
#'   }
#'   
#'   # ── Base plot ───────────────────────────────────────────────────────────────
#'   plt <- ggplot2::ggplot(input_df, base_aes) +
#'     ggplot2::theme_bw(base_size = txt_size) +
#'     ggplot2::theme(
#'       text              = ggplot2::element_text(size = txt_size),
#'       legend.position   = leg_loc
#'     ) +
#'     ggplot2::xlab(xlab) +
#'     ggplot2::ylab(ylab) +
#'     My_Theme(txt_size = txt_size) +
#'     ggplot2::guides(
#'       color = ggplot2::guide_legend(override.aes = list(size = 10, pch = 20)),
#'       fill  = "none"
#'     )
#'   
#'   if (!is.null(ribbon_layer)) plt <- plt + ribbon_layer
#'   
#'   # ── Marker ──────────────────────────────────────────────────────────────────
#'   plt <- plt + if (marker == "line") {
#'     ggplot2::geom_path(size = marker_size)
#'   } else {
#'     ggplot2::geom_point(size = marker_size)
#'   }
#'   
#'   # ── Colours ─────────────────────────────────────────────────────────────────
#'   if (!is.null(cols)) {
#'     plt <- plt +
#'       ggplot2::scale_colour_manual(leg_title, values = cols, drop = TRUE) +
#'       ggplot2::scale_fill_manual(leg_title, values = cols, drop = TRUE)
#'   }
#'   
#'   # ── Legend title ─────────────────────────────────────────────────────────────
#'   if (is.null(leg_title)) {
#'     plt <- plt + ggplot2::theme(legend.title = ggplot2::element_blank())
#'   }
#'   
#'   # ── Facet ───────────────────────────────────────────────────────────────────
#'   if (!is.null(facet)) {
#'     plt <- plt + ggplot2::facet_grid(stats::reformulate(facet, "."))
#'   }
#'   
#'   # ── Axis limits (applied before log so limits are in data space) ─────────────
#'   if (!is.null(plot_xlim)) {
#'     plt <- plt +
#'       ggplot2::scale_x_continuous(
#'         breaks = stefBreaks(plot_xlim)[[1]],
#'         labels = stefBreaks(plot_xlim)[[2]]
#'       ) +
#'       ggplot2::coord_cartesian(xlim = plot_xlim)
#'   }
#'   
#'   if (!is.null(plot_ylim)) {
#'     plt <- plt +
#'       ggplot2::scale_y_continuous(
#'         breaks = stefBreaks(plot_ylim)[[1]],
#'         labels = stefBreaks(plot_ylim)[[2]]
#'       ) +
#'       ggplot2::coord_cartesian(ylim = plot_ylim)
#'   }
#'   
#'   # ── Log axes (override limits above if requested) ────────────────────────────
#'   log <- tolower(log)
#'   if (log == "x") plt <- plt + ggplot2::scale_x_log10()
#'   if (log == "y") plt <- plt + ggplot2::scale_y_log10()
#'   
#'   # ── Suppress scientific notation by default ──────────────────────────────────
#'   if (!isTRUE(sci_notation)) {
#'     plt <- plt +
#'       ggplot2::scale_x_continuous(
#'         labels = function(x) format(x, scientific = FALSE)
#'       ) +
#'       ggplot2::scale_y_continuous(
#'         labels = function(y) format(y, scientific = FALSE)
#'       )
#'   }
#'   
#'   return(plt)
#' }

#' General-purpose xy scatter / line plot
#'
#' Builds a ggplot2 line or point plot from a data frame. Supports error
#' ribbons or bars, manual colours, faceting, log axes, and axis limits.
#'
#' @param input_df    data.frame  Must contain columns named by `xvar` and
#'                     `yvar`. For error display it also needs a `stdev` column.
#' @param xvar        character  Column name for the x axis (default: "x").
#' @param yvar        character  Column name for the y axis (default: "y").
#' @param xlab        character/NULL  X axis label. NULL defaults to `xvar`.
#' @param ylab        character/NULL  Y axis label. NULL defaults to `yvar`.
#' @param colour_by   character  Column used for colour/fill grouping
#'                     (default: "id").
#' @param leg_loc     character  Legend position: "right", "left", "top",
#'                     "bottom", or "none" (default: "right").
#' @param error       character  Error display: "full" → full ±SD,
#'                     "half" → ±0.5 SD, "none" → no error (default: "none").
#' @param error_type  character  "ribbon" (default) or "errorbar".
#' @param cols        character/NULL  Manual colour values (recycled if short).
#' @param facet       character/NULL  Column name to facet by (rows ~ facet).
#' @param marker      character  "line" (default) or "point".
#' @param marker_size numeric  Geom size (default: 1).
#' @param n_ticks     integer  Approximate number of x-axis breaks (default: 10).
#' @param log         character  Log-scale axes: "x", "y", or "none" (default).
#' @param sci_notation logical  TRUE enables scientific notation on axes
#'                     (default: FALSE — plain number format).
#' @param leg_title   character/NULL  Legend title. NULL removes the title.
#' @param plot_ylim   numeric/NULL  Y axis limits c(min, max).
#' @param plot_xlim   numeric/NULL  X axis limits c(min, max).
#' @param txt_size    numeric  Base font size passed to My_Theme (default: 12).
#'
#' @return ggplot2 object.
plotr <- function(input_df,
                  xvar         = "x",
                  yvar         = "y",
                  xlab         = NULL,
                  ylab         = NULL,
                  colour_by    = "id",
                  leg_loc      = "right",
                  error        = "none",
                  error_type   = "ribbon",
                  cols         = NULL,
                  facet        = NULL,
                  marker       = "line",
                  marker_size  = 1,
                  n_ticks      = 10,
                  log          = "none",
                  sci_notation = FALSE,
                  leg_title    = NULL,
                  plot_ylim    = NULL,
                  plot_xlim    = NULL,
                  txt_size     = 12) {
  
  if (is.null(xlab)) xlab <- xvar
  if (is.null(ylab)) ylab <- yvar
  
  base_aes <- ggplot2::aes_string(
    x      = xvar,
    y      = yvar,
    colour = colour_by,
    fill   = colour_by
  )
  
  # ── Error display ───────────────────────────────────────────────────────────
  err_mode <- tolower(error)
  ribbon_layer <- NULL
  
  if (err_mode %in% c("full", "y", "yes")) {
    input_df$ymax <- input_df[[yvar]] + input_df$stdev
    input_df$ymin <- input_df[[yvar]] - input_df$stdev
  } else if (err_mode %in% c("half", "0.5")) {
    input_df$ymax <- input_df[[yvar]] + 0.5 * input_df$stdev
    input_df$ymin <- input_df[[yvar]] - 0.5 * input_df$stdev
  }
  
  if (err_mode != "none" && err_mode != "n") {
    err_aes <- ggplot2::aes_string(ymax = "ymax", ymin = "ymin")
    ribbon_layer <- if (tolower(error_type) == "ribbon") {
      ggplot2::geom_ribbon(err_aes, alpha = 0.25, colour = NA)
    } else {
      ggplot2::geom_errorbar(err_aes, alpha = 0.5, colour = "black",
                             width = marker_size / 5)
    }
  }
  
  # ── Base plot ───────────────────────────────────────────────────────────────
  plt <- ggplot2::ggplot(input_df, base_aes) +
    ggplot2::theme_bw(base_size = txt_size) +
    ggplot2::theme(
      text            = ggplot2::element_text(size = txt_size),
      legend.position = leg_loc
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    My_Theme(txt_size = txt_size) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 10, pch = 20)),
      fill  = "none"
    )
  
  if (!is.null(ribbon_layer)) plt <- plt + ribbon_layer
  
  # ── Marker ──────────────────────────────────────────────────────────────────
  plt <- plt + if (marker == "line") {
    ggplot2::geom_path(size = marker_size)
  } else {
    ggplot2::geom_point(size = marker_size)
  }
  
  # ── Colours ─────────────────────────────────────────────────────────────────
  if (!is.null(cols)) {
    plt <- plt +
      ggplot2::scale_colour_manual(leg_title, values = cols, drop = TRUE) +
      ggplot2::scale_fill_manual(leg_title, values = cols, drop = TRUE)
  }
  
  # ── Legend title ────────────────────────────────────────────────────────────
  if (is.null(leg_title)) {
    plt <- plt + ggplot2::theme(legend.title = ggplot2::element_blank())
  }
  
  # ── Facet ───────────────────────────────────────────────────────────────────
  if (!is.null(facet)) {
    plt <- plt + ggplot2::facet_grid(stats::reformulate(facet, "."))
  }
  
  # ── Axis scales (limits, log, sci notation — built once per axis) ───────────
  make_scale <- function(axis, is_log, lims, use_sci) {
    scale_fn <- if (is_log) {
      if (axis == "x") ggplot2::scale_x_log10 else ggplot2::scale_y_log10
    } else {
      if (axis == "x") ggplot2::scale_x_continuous else ggplot2::scale_y_continuous
    }
    
    args <- list()
    
    if (!is_log && !is.null(lims)) {
      sb           <- stefBreaks(lims)
      args$breaks  <- sb[[1]]
      args$labels  <- sb[[2]]
    } else if (!use_sci) {
      args$labels <- function(v) format(v, scientific = FALSE)
    }
    
    do.call(scale_fn, args)
  }
  
  log <- tolower(log)
  plt <- plt +
    make_scale("x", log == "x", plot_xlim, isTRUE(sci_notation)) +
    make_scale("y", log == "y", plot_ylim, isTRUE(sci_notation))
  
  if (!is.null(plot_xlim)) plt <- plt + ggplot2::coord_cartesian(xlim = plot_xlim)
  if (!is.null(plot_ylim)) plt <- plt + ggplot2::coord_cartesian(ylim = plot_ylim)
  
  return(plt)
}

#' Bar chart with optional error bars and value labels
#'
#' @param input_df    data.frame  Requires columns: id (x grouping), y.
#'                     For error bars also needs: stdev.
#' @param labels      character/NULL  Custom x-axis tick labels (in id order).
#' @param cols        character/NULL  Fill colours (one per unique id).
#' @param ylab        character/NULL  Y axis label.
#' @param error       character  "full" → full error bars, "half" → ±0.5 SD,
#'                     "none" (default) → no bars.
#' @param marker_size numeric  Bar border and errorbar line size (default: 1).
#' @param value_labels logical  TRUE adds numeric value labels to bar ends.
#' @param label_hjust numeric  Horizontal adjustment for value labels (default: 5).
#' @param txt_size    numeric  Base font size (default: 12).
#'
#' @return ggplot2 object.
barPlotr <- function(input_df,
                     labels       = NULL,
                     cols         = NULL,
                     ylab         = NULL,
                     error        = "none",
                     marker_size  = 1,
                     value_labels = FALSE,
                     label_hjust  = 5,
                     txt_size     = 12,
                     colour_by    = NULL) {
  
  err_mode <- tolower(error)
  
  colour_col <- .resolve_colour_col(
    input_df,
    colour_by,
    default_order = c("id", "mem", "ins", "pore")
  )
  
  if (err_mode %in% c("full", "y", "yes")) {
    input_df$ymax <- input_df$y + input_df$stdev
    input_df$ymin <- input_df$y
  } else if (err_mode %in% c("half", "0.5")) {
    input_df$ymax <- input_df$y + 0.5 * input_df$stdev
    input_df$ymin <- input_df$y
  }
  
  plt <- ggplot2::ggplot(
    input_df,
    ggplot2::aes(
      x      = id,
      y      = y,
      colour = .data[[colour_col]],
      fill   = .data[[colour_col]]
    )
  )
  
  #pd2 <- ggplot2::position_dodge2(width = 0.9, preserve = "single")
  pd <- ggplot2::position_dodge(width = 0.9, preserve = "single")
  
  
  
  if (err_mode != "none" && err_mode != "n") {
    plt <- plt + ggplot2::geom_errorbar(
      ggplot2::aes(ymax = ymax, ymin = ymin),
      position = pd,
      width    = 0.2,
      colour   = "black",
      size     = marker_size/3
    )
  }
  
  plt <- plt +
    ggplot2::geom_bar(
      stat     = "identity",
      position = pd,
      size     = marker_size
    ) +
    ggplot2::guides(fill = "none", color = "none") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme_bw(base_size = txt_size) +
    ggplot2::xlab("") +
    ggplot2::ylab(ylab) +
    My_Theme(txt_size = txt_size)
  
  if (!is.null(cols)) {
    
    fill_scale <- ggplot2::scale_fill_manual(values = cols)
    colour_scale <- ggplot2::scale_colour_manual(values = cols)
    
    if (!is.null(labels) && colour_col == "id") {
      fill_scale <- ggplot2::scale_fill_manual(
        values = cols,
        labels = labels
      )
    }
    
    plt <- plt + fill_scale + colour_scale
  }
  
  if (isTRUE(value_labels)) {
    plt <- plt + ggplot2::geom_text(
      ggplot2::aes(label = round(y, 1)),
      hjust  = label_hjust,
      colour = "black"
    )
  }
  
  return(plt)
}

#' Resolve which column to use for colour/fill aesthetics.
.resolve_colour_col <- function(df, colour_by, default_order) {
  if (!is.null(colour_by)) {
    if (!colour_by %in% names(df))
      stop("colour_by column '", colour_by, "' not found in data. ",
           "Available columns: ", paste(names(df), collapse = ", "))
    return(colour_by)
  }
  for (col in default_order) {
    if (col %in% names(df)) return(col)
  }
  warning(".resolve_colour_col: none of the default columns (",
          paste(default_order, collapse = ", "),
          ") found in data — using first column as group.")
  return(names(df)[1L])
}

#' Save a ggplot to PNG with a forced white background
#'
#' Defaults match Nature double-column figures (183 × 170 mm, 300 dpi).
#'
#' @param plot    ggplot  Plot object to save.
#' @param direc   character  Directory to save into.
#' @param name    character  Output filename stem (no extension needed).
#' @param width   numeric  Width in mm (default: 183).
#' @param height  numeric  Height in mm (default: 170).
#' @param dpi     numeric  Resolution (default: 300).
pngOut <- function(plot, direc, name, width = 183, height = 170, dpi = 300) {
  
  out <- cowplot::ggdraw(plot) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  ggplot2::ggsave(
    filename = file.path(direc, paste0(name, ".png")),
    plot     = out,
    device   = "png",
    width    = width,
    height   = height,
    units    = "mm",
    dpi      = dpi,
    limitsize = FALSE
  )
}


# ── Section 5: Signal Processing ─────────────────────────────────────────────

#' First derivative of an xy dataset via smoothing spline
#'
#' Fits a smoothing spline to the data and evaluates its first derivative
#' at each original x value.
#'
#' @param x  numeric  x-axis values.
#' @param y  numeric  y-axis values.
#'
#' @return numeric vector of first-derivative values at each x.
derivative1 <- function(x, y) {
  model <- stats::smooth.spline(x = x, y = y)
  stats::predict(model, x = x, deriv = 1)$y
}


#' Detect local maxima in a numeric vector
#'
#' Identifies peaks by locating sign changes in the second difference.
#' Optional smoothing reduces noise sensitivity, and a prominence threshold
#' filters out weak peaks.
#'
#' @param y               numeric  Input signal vector.
#' @param smooth          integer  Moving-average window width. 0 = no smoothing
#'                          (default: 0).
#' @param min_prominence  numeric  Minimum peak prominence on a [0, 1] scale
#'                          relative to the total signal range (default: 0.01).
#'
#' @return data.frame with columns: index, value, smoothed_value, prominence.
#'   Returns NULL (with a warning) if no peaks pass the threshold.
#'
#' @examples
#' peaks <- peakID(my_spectrum, smooth = 5, min_prominence = 0.05)
peakID <- function(y, smooth = 0, min_prominence = 0.01) {
  
  # Optional moving-average smoothing
  if (smooth > 1) {
    kernel <- rep(1 / smooth, smooth)
    y_s    <- as.numeric(stats::filter(y, kernel, sides = 2))
    y_s[is.na(y_s)] <- y[is.na(y_s)]   # preserve edge values
  } else {
    y_s <- y
  }
  
  # Local maxima: positions where second difference changes sign from + to -
  peak_idx <- which(diff(sign(diff(y_s))) == -2) + 1
  
  if (length(peak_idx) == 0) {
    warning("peakID: no peaks detected.")
    return(NULL)
  }
  
  # Normalised prominence: height above the global minimum
  prom      <- (y_s[peak_idx] - min(y_s)) / (max(y_s) - min(y_s))
  keep_idx  <- peak_idx[prom >= min_prominence]
  
  if (length(keep_idx) == 0) {
    warning("peakID: no peaks exceed prominence threshold (", min_prominence, ").")
    return(NULL)
  }
  
  data.frame(
    index          = keep_idx,
    value          = y[keep_idx],
    smoothed_value = y_s[keep_idx],
    prominence     = prom[prom >= min_prominence]
  )
}


#' Refine a peak position by fitting a local Gaussian
#'
#' Takes a coarsely detected peak index, extracts a window of points around
#' it, and fits a Gaussian to recover a sub-sample peak position. Supports
#' log-x fitting for DLS-style size distributions.
#'
#' @param x          numeric  Full x-axis vector.
#' @param y          numeric  Full y-axis vector.
#' @param peak_idx   integer  Index of the detected peak in x/y.
#' @param window     integer  Points on each side of the peak to use
#'                    (default: 3).
#' @param show_plot  logical  TRUE plots the data and fitted Gaussian
#'                    (default: TRUE).
#' @param log_x      logical  TRUE fits in log(x) space, returning exp(mu)
#'                    as the peak position (default: FALSE).
#'
#' @return numeric  Refined x-coordinate of the peak maximum. Returns
#'   x[peak_idx] if the Gaussian fit fails.
peakValue <- function(x, y, peak_idx, window = 3, show_plot = TRUE,
                      log_x = FALSE) {
  
  n   <- length(x)
  idx <- seq(max(1, peak_idx - window), min(n, peak_idx + window))
  x_w <- x[idx]
  y_w <- y[idx]
  
  if (log_x) {
    message("peakValue: fitting in log(x) space")
    x_fit  <- log(x_w)
    mu0    <- log(x[peak_idx])
    x_lbl  <- "x (log scale)"
  } else {
    x_fit  <- x_w
    mu0    <- x[peak_idx]
    x_lbl  <- "x"
  }
  
  y0    <- min(y_w)
  A0    <- max(y_w) - y0
  sig0  <- (max(x_fit) - min(x_fit)) / 4
  
  fit <- try(
    stats::nls(
      y_w ~ y0 + A * exp(-(x_fit - mu)^2 / (2 * sigma^2)),
      start   = list(y0 = y0, A = A0, mu = mu0, sigma = sig0),
      control = stats::nls.control(warnOnly = TRUE)
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    warning("peakValue: Gaussian fit failed at index ", peak_idx,
            " — returning raw peak position.")
    return(x[peak_idx])
  }
  
  p          <- stats::coef(fit)
  refined_x  <- if (log_x) exp(p["mu"]) else p["mu"]
  
  if (isTRUE(show_plot)) {
    x_dense <- seq(min(x_fit), max(x_fit), length.out = 200)
    y_dense <- p["y0"] + p["A"] * exp(-(x_dense - p["mu"])^2 / (2 * p["sigma"]^2))
    
    graphics::plot(x_w, y_w, pch = 19, col = "blue",
                   main = paste0("Peak at index ", peak_idx),
                   xlab = x_lbl, ylab = "y",
                   log  = if (log_x) "x" else "")
    graphics::lines(if (log_x) exp(x_dense) else x_dense, y_dense,
                    col = "red", lwd = 2)
    graphics::abline(v = refined_x, col = "darkgreen", lty = 2)
    graphics::legend("topright",
                     legend = c("Data", "Gaussian fit", "Refined peak"),
                     col    = c("blue", "red", "darkgreen"),
                     lty    = c(NA, 1, 2),
                     pch    = c(19, NA, NA))
  }
  
  return(as.numeric(refined_x))
}


#' Find local maxima using second-difference sign change (legacy helper)
#'
#' Credit: Stasia Grinberg. For valley detection, pass `-y`.
#'
#' @param x  numeric  Signal vector.
#' @param m  integer  Sensitivity — larger values detect fewer, broader peaks
#'                     (default: 3).
#'
#' @return integer vector of peak indices.
find_peaks <- function(x, m = 3) {
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks   <- sapply(which(shape < 0), function(i) {
    lo <- max(1, i - m + 1)
    hi <- min(length(x), i + m + 1)
    if (all(x[c(lo:i, (i + 2):hi)] <= x[i + 1])) i + 1 else numeric(0)
  })
  unlist(pks)
}


#' Compute the area under each peak in an xy signal
#'
#' Locates peaks with find_peaks, finds each peak's left and right minima,
#' and integrates using the trapezoidal rule. Optionally operates in log10(x)
#' space for size distributions.
#'
#' @param x         numeric  x-axis values (must be sorted or will be sorted).
#' @param y         numeric  y-axis values.
#' @param log_x     logical  TRUE log10-transforms x before integration
#'                   (default: FALSE).
#' @param ...       Arguments forwarded to find_peaks (e.g. m for sensitivity).
#'
#' @return data.frame with columns: peak_index, peak_position, area, perc_area.
#'   Returns NULL if no peaks are detected.
peakArea <- function(x, y, log_x = FALSE, ...) {
  
  ord <- order(x)
  x   <- x[ord]
  y   <- y[ord]
  
  peak_indices <- find_peaks(y, ...)
  if (is.null(peak_indices) || length(peak_indices) == 0) {
    message("peakArea: no peaks detected.")
    return(NULL)
  }
  
  if (log_x) x <- log10(x[x > 0])
  
  result <- lapply(seq_along(peak_indices), function(i) {
    idx   <- peak_indices[i]
    left  <- if (idx > 1)          which.min(y[seq_len(idx)])              else 1L
    right <- if (idx < length(y))  which.min(y[idx:length(y)]) + idx - 1L else length(y)
    
    area <- pracma::trapz(x[left:right], y[left:right])
    data.frame(
      peak_index    = idx,
      peak_position = if (log_x) 10^x[idx] else x[idx],
      area          = area
    )
  })
  
  out             <- do.call(rbind, result)
  out$perc_area   <- 100 * out$area / sum(out$area)
  return(out)
}


# ── Section 6: Statistical Tools ─────────────────────────────────────────────

#' Rapid assessment of one or more linear model fits
#'
#' Prints R², adjusted R², RMSE, MAE, F-statistic, and ANOVA p-value for
#' each model, plus Cook's distance and residual plots.
#'
#' @param models_list  lm object or list of lm objects.
modelAssessment <- function(models_list) {
  
  if (!is.list(models_list)) models_list <- list(models_list)
  
  results <- data.frame(
    Model         = character(0),
    R_squared     = numeric(0),
    Adj_R_squared = numeric(0),
    RMSE          = numeric(0),
    MAE           = numeric(0),
    F_statistic   = numeric(0),
    ANOVA_p_value = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(models_list)) {
    model   <- models_list[[i]]
    sm      <- summary(model)
    resids  <- stats::resid(model)
    
    results[i, ] <- list(
      paste("Model", i),
      round(sm$r.squared, 3),
      round(sm$adj.r.squared, 3),
      round(sqrt(mean(resids^2)), 3),
      round(mean(abs(resids)), 3),
      round(sm$fstatistic[1], 3),
      round(stats::anova(model)$`Pr(>F)`[1], 3)
    )
    
    graphics::plot(stats::cooks.distance(model), pch = 20,
                   main = paste("Cook's Distance — Model", i))
    graphics::plot(model$residuals ~ model$fitted.values,
                   xlab = "Fitted Values", ylab = "Residuals",
                   main = paste("Residual Plot — Model", i))
  }
  
  message("Good fit: R² high | RMSE & MAE low | F-stat significant (p < 0.05)")
  print(results)
  invisible(results)
}

#' Fit dose-response models to xy data with automatic model selection
#'
#' Tries a library of drc models and selects the best by AIC. Returns a
#' smooth predicted curve suitable for overlaying on a plotr() output.
#'
#' @param input_df   data.frame  Requires x and y columns. Optional id column
#'                    for per-group fitting.
#' @param model      drc model  Starting/fallback model if auto_fit is FALSE
#'                    (default: LL.2()).
#' @param auto_fit   logical  TRUE (default) tries all candidate models and
#'                    picks the best by AIC.
#' @param force_zero logical  TRUE (default) shifts predicted y so the minimum
#'                    is at zero.
#' @param split_by_id logical  TRUE fits a separate model per entry in the
#'                    `id` column (default: FALSE).
#' @param xlimits    numeric/NULL  c(xmin, xmax) for prediction range. NULL
#'                    uses c(0, max(x)).
#'
#' @return data.frame with columns: x (dense grid), y (predicted), lower
#'         (95% CI lower bound), upper (95% CI upper bound), id.
fitr <- function(input_df,
                 model       = drc::LL.2(),
                 auto_fit    = TRUE,
                 force_zero  = TRUE,
                 split_by_id = FALSE,
                 xlimits     = NULL) {
  
  candidates <- list(
    MM.2  = drc::MM.2(),  MM.3  = drc::MM.3(),
    LL.2  = drc::LL.2(),  LL.3  = drc::LL.3(),  LL.3u = drc::LL.3u(),
    LL.4  = drc::LL.4(),  LL.5  = drc::LL.5(),
    LL2.2 = drc::LL2.2(), LL2.3 = drc::LL2.3(), LL2.3u = drc::LL2.3u(),
    LL2.4 = drc::LL2.4(), LL2.5 = drc::LL2.5(),
    W1.2  = drc::W1.2(),  W1.3  = drc::W1.3(),  W1.3u = drc::W1.3u(),
    W1.4  = drc::W1.4()
  )
  
  fit_one <- function(df_sub) {
    best_fit  <- NULL
    best_name <- "LL.2"
    best_aic  <- Inf
    chosen    <- model
    
    if (isTRUE(auto_fit)) {
      for (nm in names(candidates)) {
        attempt <- tryCatch(
          drc::drm(y ~ x, data = df_sub, fct = candidates[[nm]]),
          error = function(e) NULL
        )
        if (!is.null(attempt)) {
          a <- AIC(attempt)
          if (a < best_aic) {
            best_aic  <- a
            best_fit  <- attempt
            best_name <- nm
            chosen    <- candidates[[nm]]
          }
        }
      }
    } else {
      best_fit <- tryCatch(
        drc::drm(y ~ x, data = df_sub, fct = chosen),
        error = function(e) NULL
      )
      if (!is.null(best_fit)) best_aic <- AIC(best_fit)
    }
    
    if (is.null(best_fit)) {
      warning("fitr: no model converged for id '",
              if ("id" %in% names(df_sub)) df_sub$id[1] else "unknown", "'")
      return(NULL)
    }
    
    message("fitr: best model = ", best_name, "  AIC = ", round(best_aic, 2))
    print(stats::coef(best_fit))
    
    x_range <- if (is.null(xlimits)) {
      seq(0, max(df_sub$x, na.rm = TRUE), length.out = 1000)
    } else {
      seq(xlimits[1], xlimits[2], length.out = 1000)
    }
    
    pred     <- data.frame(x = x_range)
    pred_mat <- stats::predict(best_fit, newdata = pred,
                               interval = "confidence", level = 0.95)
    pred$y     <- pred_mat[, "Prediction"]
    pred$lower <- pred_mat[, "Lower"]
    pred$upper <- pred_mat[, "Upper"]
    
    if (isTRUE(force_zero)) {
      baseline   <- pred$y[which.min(pred$x)]
      pred$y     <- pred$y     - baseline
      pred$lower <- pred$lower - baseline
      pred$upper <- pred$upper - baseline
    }
    
    pred$id <- if ("id" %in% names(df_sub)) df_sub$id[1] else NA
    pred
  }
  
  df_list <- if (isTRUE(split_by_id) && "id" %in% names(input_df)) {
    split(input_df, input_df$id)
  } else {
    list(input_df)
  }
  
  out <- do.call(rbind, Filter(Negate(is.null), lapply(df_list, fit_one)))
  names(out) <- c("x", "y", "lower", "upper", "id")
  return(out)
}

# ── Section 7: Efflux / Flux Data ────────────────────────────────────────────

#' Prepare cuvette efflux data for plotting
#'
#' Removes pre-injection data, resets time to zero, normalises, and then
#' averages replicates identified by the `sample` column. Handles the ".1",
#' ".2" suffixes that some fluorimeters append to repeated sample names.
#'
#' @param input_df    data.frame  Output of cuvette_framer / cuvette_importer.
#'                     Requires columns: x (time), y (signal), sample, id.
#' @param inject_list numeric vector  Injection time(s) in the same units as x.
#'                     If fewer times than samples are given, the list is
#'                     recycled.
#' @param normalise   character  Normalisation applied after time-zeroing:
#'                     "yes"/"y" → 0-1 × 100, "zero"/"0" → subtract first
#'                     value, "min" → subtract minimum (default: "yes").
#' @param time_round  integer  Decimal places to round x to before aggregating
#'                     replicates — reduce if traces are noisy (default: 1).
#'
#' @return data.frame with columns: id (sample name), x (time), y (mean),
#'   stdev.
effluxDataPrep <- function(input_df,
                           inject_list,
                           normalise  = "yes",
                           time_round = 1) {
  
  samples   <- unique(input_df$sample)
  n_samp    <- length(samples)
  n_inj     <- length(inject_list)
  
  message("effluxDataPrep: processing ", n_samp, " sample(s): ",
          paste(samples, collapse = ", "))
  message("  Injection times: ", paste(inject_list, collapse = ", "))
  
  # Recycle injection times if fewer than sample count
  if (n_samp > n_inj) {
    inject_list <- rep(inject_list, length.out = n_samp)
    message("  Recycled injection times to match ", n_samp, " samples.")
  }
  
  # Unique split key: sample + file id to handle same sample name across files
  input_df$split_key <- paste0(input_df$sample, input_df$id)
  
  working_df <- do.call(rbind, mapply(
    function(sub_df, t_inject) {
      sub_df  <- subset(sub_df, sub_df$x > t_inject)
      sub_df$x <- sub_df$x - sub_df$x[1]
      sub_df$x <- round(sub_df$x, time_round)
      
      nrm <- tolower(normalise)
      if (nrm %in% c("y", "yes")) {
        sub_df$y <- norm(sub_df$y) * 100
      } else if (nrm %in% c("0", "zero")) {
        sub_df$y <- sub_df$y - sub_df$y[1]
      } else if (nrm == "min") {
        sub_df$y <- sub_df$y - min(sub_df$y, na.rm = TRUE)
      }
      
      sub_df
    },
    split(input_df, input_df$split_key),
    inject_list,
    SIMPLIFY = FALSE
  ))
  
  # Strip fluorimeter-added numeric suffixes (e.g. "Sample.1" → "Sample")
  cleaned_names <- unique(sub("\\.[0-9]+$", "", unique(working_df$sample)))
  message("effluxDataPrep: cleaned sample names: ",
          paste(cleaned_names, collapse = ", "))
  
  working_df$sample <- sub("\\.[0-9]+$", "", working_df$sample)
  
  # Average replicates per sample × time using modal average + SD
  out_df <- stats::aggregate(y ~ sample + x, data = working_df,
                             function(vals) c(avg = modalAvg(vals), sd = sd(vals)))
  out_df <- do.call(data.frame, out_df)
  names(out_df) <- c("id", "x", "y", "stdev")
  
  message("effluxDataPrep: done — ", nrow(out_df), " rows returned")
  return(out_df)
}


#' Background-corrected maximum efflux values for pairwise sample/background
#'
#' Extracts maximum y values up to `t_max` for four groups (bkg1, sample1,
#' bkg2, sample2) and returns background-subtracted values with propagated SD.
#'
#' @param input_df  data.frame  Output of effluxDataPrep. Requires: x, y,
#'                   stdev, id.
#' @param labels    character vector  Exactly four group names in order:
#'                   c(bkg1, sample1, bkg2, sample2).
#' @param t_max     numeric  Upper time limit for maximum extraction.
#'
#' @return data.frame with columns: id (sample1, sample2), y (corrected max),
#'   stdev (propagated).
maxFlux <- function(input_df, labels, t_max) {
  
  stopifnot(length(labels) == 4)
  
  sub_df <- subset(input_df, input_df$x < t_max)
  
  get_max_row <- function(group_id) {
    g <- subset(sub_df, sub_df$id == group_id)
    g[which.max(g$y), , drop = FALSE]
  }
  
  bkg1  <- get_max_row(labels[1])
  dat1  <- get_max_row(labels[2])
  bkg2  <- get_max_row(labels[3])
  dat2  <- get_max_row(labels[4])
  
  rbind(
    data.frame(id    = labels[2],
               y     = dat1$y - bkg1$y,
               stdev = sqrt(dat1$stdev^2 + bkg1$stdev^2)),
    data.frame(id    = labels[4],
               y     = dat2$y - bkg2$y,
               stdev = sqrt(dat2$stdev^2 + bkg2$stdev^2))
  )
}


#' Estimate initial efflux rate from the linear phase of each trace
#'
#' Fits a linear model over several short time windows and returns their
#' mean slope as the initial rate estimate.
#'
#' @param input_df     data.frame  Single-sample trace with columns: x, y, id.
#' @param rate_windows numeric vector  Upper time limits (from t = 0) defining
#'                      the early linear windows (default: c(1, 1.5, 2.5)).
#'
#' @return data.frame with columns: y (mean slope), stdev (SD across windows),
#'   id. Returns NULL if any window has no data.
initRate <- function(input_df, rate_windows = c(1, 1.5, 2.5)) {
  
  slopes <- sapply(rate_windows, function(t_end) {
    sub_df <- subset(input_df, x < t_end)
    
    if (nrow(sub_df) == 0) {
      message("initRate: no data for x < ", t_end, " in id '",
              if ("id" %in% names(input_df)) unique(input_df$id) else "?", "'")
      return(NULL)
    }
    if (nrow(sub_df) == 1) {
      message("initRate: only 1 point for x < ", t_end, " — returning NA slope")
      return(NA_real_)
    }
    
    stats::coef(stats::lm(y ~ x, data = sub_df))[2]
  })
  
  # If any window produced NULL (no data), skip the whole sample
  if (any(vapply(slopes, is.null, logical(1)))) return(NULL)
  
  data.frame(
    y     = mean(unlist(slopes), na.rm = TRUE),
    stdev = sd(unlist(slopes),   na.rm = TRUE),
    id    = unique(input_df$id)
  )
}


# ── Section 8: Spectroscopy & Quantification ─────────────────────────────────

#' Quantify polymer concentration from a UV absorbance calibration curve
#'
#' Loads a stored calibration curve, fits a linear model, and applies it
#' to the input data at the specified quantification wavelength.
#'
#' Currently only DIBMA is supported.
#'
#' @param input_df  data.frame  Requires x (wavelength) and y (absorbance)
#'                   columns. An `id` column enables per-sample quantification.
#' @param quant_wl  numeric  Wavelength to quantify at (default: 215 nm).
#' @param poltype   character  Polymer type: currently "DIBMA" only.
#' @param dil_facs  numeric/vector  Dilution factor(s). If a single value is
#'                   given it is applied to all samples.
#' @param show_plot logical  TRUE prints a calibration curve plot with
#'                   confidence interval (default: FALSE).
#'
#' @return data.frame: input_df rows at `quant_wl` with added columns:
#'   uM (estimated concentration) and corrected_uM (× dil_facs).
polQuant <- function(input_df,
                     quant_wl  = 215,
                     poltype   = "DIBMA",
                     dil_facs  = 1,
                     show_plot = FALSE) {
  
  # ── Load calibration data ────────────────────────────────────────────────────
  if (poltype == "DIBMA") {
    cal_dir   <- "~/Onedrive/Lab_Data/UV-VIS/SCANS/20210916_DIBMA_C_Curve/"
    cal_files <- c("10um", "7p5", "5um1", "3p5", "2p5um", "1p25um")
    cal_concs <- c(10, 7.5, 5, 3.5, 2.5, 1.25)
    
    cal_df      <- data_framr(cal_files, cal_dir, normalise = "zero")
    cal_df$id   <- cal_concs[match(cal_df$id, cal_files)]
    cal_df      <- subset(cal_df, cal_df$x == quant_wl)
    cal_df$x    <- cal_df$id   # x = concentration for model
  } else {
    stop("polQuant: only 'DIBMA' is currently supported.")
  }
  
  fit       <- stats::lm(y ~ x, data = cal_df)
  intercept <- stats::coef(fit)[1]
  slope     <- stats::coef(fit)[2]
  message("polQuant: calibration fit  intercept = ", round(intercept, 4),
          "  slope = ", round(slope, 4))
  
  if (isTRUE(show_plot)) {
    fitted_ci     <- as.data.frame(stats::predict(fit, interval = "confidence"))
    cal_plot_df   <- cbind(cal_df, fitted_ci)
    plt <- plotr(cal_plot_df, xlab = "Concentration (µM)",
                 ylab = paste0("A", quant_wl), marker = "point") +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lwr, ymax = upr),
        fill = "black", alpha = 0.2, colour = NA
      ) +
      ggplot2::geom_line(ggplot2::aes(y = fit), size = 0.5, colour = "black")
    print(plt)
  }
  
  out      <- subset(input_df, input_df$x == quant_wl)
  out$uM   <- (out$y - intercept) / slope
  try(out$corrected_uM <- out$uM * dil_facs)
  return(out)
}


#' Quantify DNA concentration from A260 absorbance
#'
#' Applies the standard Beer-Lambert conversion for double- and single-
#' stranded DNA, with optional subtraction of DIBMA absorbance at 260 nm.
#' Returns concentration in nM.
#'
#' @param input_df          data.frame  Must have x (wavelength) and y (A260)
#'                            columns. x must contain 260.
#' @param dil_fac           numeric  Dilution factor applied before measurement.
#' @param mol_weight        numeric  DNA molecular weight in g/mol (obtain from
#'                            IDT OligoAnalyzer or calculate from sequence).
#' @param perc_ss           numeric/vector  Percentage single-stranded content
#'                            (0–100). Scalar or one value per sample row.
#' @param dibma_conc        numeric/vector  DIBMA concentration in the sample
#'                            in µM (used to subtract the polymer A260
#'                            contribution). 0 = no DIBMA (default).
#'
#' @return input_df with added columns: undiluted, pol_a260, to_subtract,
#'   corr, ds_ug, ss_ug, ugml, ugul, dna_mol_weight (= mol_weight), M, nM.
DNAQuant <- function(input_df,
                     dil_fac    = 1,
                     mol_weight,
                     perc_ss    = 0,
                     dibma_conc = 0) {
  
  DS_FACTOR <- 50   # µg/mL per A260 unit for dsDNA
  SS_FACTOR <- 33   # µg/mL per A260 unit for ssDNA
  
  # Load DIBMA reference spectrum for A260 subtraction
  dib_path  <- "/media/jonah/Lickitung/SYNC/8-LAB_DATA/UV-VIS/SCANS/20220407_New_ND_Quant/"
  dib_raw   <- xy_data_importr("dibma_1in20.csv", dib_path, ext = "",
                               normalise = "zero")
  dib_uM    <- as.numeric(polQuant(dib_raw, dil_facs = 20)[5])
  
  # Restrict to A260
  if (any(input_df$x != 260))
    input_df <- subset(input_df, input_df$x == 260)
  
  input_df$undiluted <- input_df$y * dil_fac
  
  # DIBMA A260 contribution scaled to sample concentration
  dib_a260           <- dib_raw$y[dib_raw$x == 260] * 20   # 1:20 dilution reference
  input_df$pol_a260  <- dib_a260
  
  if (length(dibma_conc) < nrow(input_df)) {
    message("DNAQuant: recycling dibma_conc to match ", nrow(input_df), " sample(s).")
    dibma_conc <- rep(dibma_conc, length.out = nrow(input_df))
  }
  
  input_df$to_subtract <- (input_df$pol_a260 / dib_uM) * dibma_conc
  input_df$corr        <- input_df$undiluted - input_df$to_subtract
  
  input_df$ds_ug  <- input_df$corr * DS_FACTOR
  input_df$ss_ug  <- input_df$corr * SS_FACTOR
  input_df$perc_ss <- perc_ss
  
  # Weighted mean of ds and ss contributions per row
  input_df$ugml <- apply(
    input_df[, c("ds_ug", "ss_ug", "perc_ss")], 1,
    function(r) stats::weighted.mean(c(r["ds_ug"], r["ss_ug"]),
                                     c(100 - r["perc_ss"], r["perc_ss"]))
  )
  
  input_df$ugul           <- input_df$ugml / 1000
  input_df$dna_mol_weight <- mol_weight
  input_df$M              <- input_df$ugul / mol_weight
  input_df$nM             <- input_df$M * 1e9
  
  return(input_df)
}


#' Quantify fluorescence using a stored calibration curve
#'
#' Reads calibration curve CSVs from a named sub-folder, smooths each
#' concentration series with a loess fit, identifies the emission maximum,
#' and uses linear interpolation (approx) to convert sample signals to µM.
#' Handles single values, single spectra, and multi-id spectrum data frames.
#'
#' @param input_data      numeric, or data.frame with x/y (and optionally id).
#' @param cal_dir         character  Sub-folder name within the calibration
#'                          curve root directory (default: "6HB0C_Cy5/").
#' @param show_plot       logical  TRUE prints the calibration curve with CI
#'                          (default: TRUE).
#' @param subset_to_em    logical  TRUE subsets the returned data to the
#'                          emission maximum wavelength only (default: TRUE).
#'
#' @return data.frame with a `uM` column containing estimated concentrations.
Quantr <- function(input_data,
                   cal_dir      = "6HB0C_Cy5/",
                   show_plot    = TRUE,
                   subset_to_em = TRUE) {
  
  CAL_ROOT <- "~/Onedrive/Lab_Data/CAL_CURVES/"
  
  # ── Discover available calibration curves ────────────────────────────────────
  avail <- list.dirs(CAL_ROOT, recursive = FALSE, full.names = FALSE)
  message("Quantr: available calibration curves: ", paste(avail, collapse = ", "))
  message("Quantr: using '", cal_dir, "'")
  
  # ── Load and smooth calibration data ─────────────────────────────────────────
  cal_path  <- file.path(CAL_ROOT, cal_dir)
  cal_files <- list.files(cal_path, pattern = "\\.csv$")
  if (length(cal_files) == 0)
    stop("Quantr: no CSV files found in '", cal_path, "'.")
  
  cal_df      <- data_framr(cal_files, endSlash(cal_path))
  cal_df$id   <- gsub("\\.", "", cal_df$id)
  cal_df$id   <- gsub("p", ".", cal_df$id)
  cal_df$id   <- as.numeric(gsub("[a-zA-Z]", "", cal_df$id))
  
  # Smooth each concentration curve with loess
  cal_df <- do.call(rbind, lapply(unique(cal_df$id), function(cid) {
    sub_df    <- cal_df[cal_df$id == cid, ]
    y_smooth  <- stats::predict(stats::loess(y ~ x, data = sub_df, span = 0.3),
                                newdata = data.frame(x = sub_df$x))
    data.frame(id = cid, x = sub_df$x, y = y_smooth)
  }))
  
  em_max   <- cal_df$x[which.max(cal_df$y)]
  curve_df <- subset(cal_df, cal_df$x == em_max)
  curve_df$x <- curve_df$id   # x = concentration for interpolation
  
  # Fit linear model (for plot only; quantification uses approx)
  cal_fit   <- stats::lm(y ~ x, data = curve_df)
  intercept <- stats::coef(cal_fit)[1]
  slope     <- stats::coef(cal_fit)[2]
  
  # ── Optional calibration plot ─────────────────────────────────────────────────
  if (isTRUE(show_plot)) {
    fitted_ci    <- as.data.frame(stats::predict(cal_fit, interval = "confidence"))
    plot_curve   <- cbind(curve_df, fitted_ci)
    
    x_sorted  <- sort(plot_curve$x[plot_curve$x > 0])
    log_scale <- (max(diff(x_sorted)) / min(diff(x_sorted))) > 10
    
    plt <- plotr(plot_curve, marker = "point",
                 xlab = "Concentration (µM)", ylab = "Fluorescence (AU)") +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lwr, ymax = upr),
        fill = "grey70", alpha = 0.4, colour = NA
      ) +
      ggplot2::geom_line(ggplot2::aes(y = fit), size = 0.5, colour = "black") +
      ggplot2::ggtitle(paste0("Calibration @ ", em_max, " nm"))
    
    if (log_scale) plt <- plt + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
    print(plt)
  }
  
  # ── Quantify input_data ───────────────────────────────────────────────────────
  quantify_signal <- function(sig) {
    stats::approx(x = curve_df$y, y = curve_df$id, xout = sig)$y
  }
  
  if (is.numeric(input_data)) {
    return(data.frame(original = input_data,
                      uM       = quantify_signal(input_data)))
  }
  
  if (!is.data.frame(input_data))
    stop("Quantr: input_data must be numeric or a data.frame.")
  
  if (!"x" %in% names(input_data) || !"y" %in% names(input_data))
    stop("Quantr: data.frame must have x and y columns.")
  
  smooth_and_quantify <- function(df_sub) {
    if (nrow(df_sub) < 5)
      warning("Quantr: id '", unique(df_sub$id),
              "' has < 5 rows — loess may be unreliable.")
    
    df_sub$y_smooth <- tryCatch(
      stats::predict(stats::loess(y ~ x, data = df_sub, span = 0.3)),
      error = function(e) rep(NA_real_, nrow(df_sub))
    )
    
    sig       <- df_sub$y_smooth[df_sub$x == em_max]
    df_sub$uM <- if (all(is.na(df_sub$y_smooth))) 0 else quantify_signal(sig)
    df_sub
  }
  
  if ("id" %in% names(input_data)) {
    result <- do.call(rbind, lapply(split(input_data, input_data$id),
                                    smooth_and_quantify))
    row.names(result) <- NULL
  } else {
    result <- smooth_and_quantify(input_data)
  }
  
  if (isTRUE(subset_to_em)) {
    result   <- subset(result, result$x == em_max)
    result$x <- NULL
  }
  
  message("Quantr: done — ", nrow(result), " sample(s) quantified")
  return(result)
}


# ── Section 9: Plate Reader Pipeline ─────────────────────────────────────────

#' Import and normalise plate reader kinetic + endpoint fluorescence data
#'
#' Reads a kinetic run file and a detergent (endpoint) file, extracts start,
#' kinetic, and end fluorescence ratios (wl1 / wl2), normalises each well
#' to [0, 100], and merges with a plate map.
#'
#' @param run_file          character  Path to the kinetic run xlsx file.
#' @param detergent_file    character  Path to the detergent endpoint xlsx file.
#' @param plate_map         data.frame  Must have columns: Well, id.
#' @param input_cols        integer vector  Column indices of wells to include
#'                            (default: 1:12).
#' @param start_wl1_rows    integer vector  Row range for start WL1 in run file.
#' @param start_wl2_rows    integer vector  Row range for start WL2 in run file.
#' @param start_cols        integer vector  Column range for start data.
#' @param end_wl1_rows      integer vector  Row range for end WL1 in detergent file.
#' @param end_wl2_rows      integer vector  Row range for end WL2 in detergent file.
#' @param end_cols          integer vector  Column range for end data.
#' @param kinetic_wl1_rows  integer vector  Row range for kinetic WL1 in run file.
#' @param kinetic_wl2_rows  integer vector  Row range for kinetic WL2 in run file.
#' @param kinetic_cols      integer vector  Column range for kinetic data.
#' @param time_col          integer  Column index of the time vector in run file.
#' @param compute_ratio     logical  TRUE (default) computes wl1/wl2 ratio.
#'                            FALSE uses wl1 alone.
#' @param include_start     logical  TRUE prepends the start-scan ratio as a
#'                            synthetic time point (default: TRUE).
#' @param t_min_origin      logical  TRUE resets each well's time to start at
#'                            the minimum fluorescence point, useful for
#'                            vesicle assays with variable lag (default: FALSE).
#'
#' @return data.frame with columns: x (time), raw (ratio), y (normalised 0–100),
#'   id (from plate_map), Well.
platR <- function(run_file,
                  detergent_file,
                  plate_map,
                  input_cols       = 1:12,
                  start_wl1_rows   = 52:59,   start_wl2_rows   = 81:88,
                  start_cols       = 2:13,
                  end_wl1_rows     = 47:54,   end_wl2_rows     = 76:83,
                  end_cols         = 2:13,
                  kinetic_wl1_rows = 137:267, kinetic_wl2_rows = 270:406,
                  kinetic_cols     = 2:99,
                  time_col         = 2,
                  compute_ratio    = TRUE,
                  include_start    = TRUE,
                  t_min_origin     = FALSE) {
  
  run_df <- openxlsx::read.xlsx(run_file)
  det_df <- openxlsx::read.xlsx(detergent_file)
  
  # ── Helper: extract a numeric matrix from an xlsx region, drop empty rows/cols
  extract_matrix <- function(df, rows, cols) {
    m <- apply(as.matrix(df[rows, cols]), 2, as.numeric)
    m <- m[rowSums(is.na(m)) < ncol(m), colSums(is.na(m)) < nrow(m), drop = FALSE]
    m
  }
  
  # ── End (detergent) data ────────────────────────────────────────────────────
  end_wl1  <- extract_matrix(det_df, end_wl1_rows, end_cols)
  end_ratio <- if (isTRUE(compute_ratio)) {
    end_wl2 <- extract_matrix(det_df, end_wl2_rows, end_cols)
    end_wl1 / end_wl2
  } else {
    end_wl1
  }
  end_ratio  <- stats::na.omit(end_ratio)
  
  # ── Well labels ─────────────────────────────────────────────────────────────
  well_grid  <- as.vector(outer(LETTERS[1:nrow(end_ratio)], input_cols,
                                function(r, c) paste0(r, c)))
  well_nums  <- as.integer(sub("[A-Z]", "", well_grid))
  wells_valid <- well_grid[well_nums %in% input_cols]
  wells_valid <- wells_valid[
    wells_valid %in% as.character(run_df[kinetic_wl1_rows[1] - 1, kinetic_cols])
  ]
  
  end_df <- data.frame(Well     = wells_valid,
                       EndRatio = as.vector(t(end_ratio)))
  
  # ── Start data (optional) ────────────────────────────────────────────────────
  if (isTRUE(include_start)) {
    start_wl1 <- extract_matrix(run_df, start_wl1_rows, start_cols)
    start_wl1 <- start_wl1[, colSums(is.na(start_wl1)) == 0, drop = FALSE]
    
    start_ratio <- if (isTRUE(compute_ratio)) {
      start_wl2 <- extract_matrix(run_df, start_wl2_rows, start_cols)
      start_wl2 <- start_wl2[, colSums(is.na(start_wl2)) == 0, drop = FALSE]
      start_wl1 / start_wl2
    } else {
      start_wl1
    }
    start_df <- data.frame(Well       = wells_valid,
                           StartRatio = as.vector(t(start_ratio)))
  }
  
  # ── Kinetic data ────────────────────────────────────────────────────────────
  col_names <- as.character(run_df[kinetic_wl1_rows[1] - 1, kinetic_cols])
  
  read_kinetic <- function(rows) {
    m           <- data.frame(lapply(run_df[rows, kinetic_cols], as.numeric))
    colnames(m) <- col_names
    m[, colSums(!is.na(m)) > 0, drop = FALSE]
  }
  
  tdat_wl1 <- read_kinetic(kinetic_wl1_rows)
  
  tdat <- if (isTRUE(compute_ratio)) {
    tdat_wl2 <- read_kinetic(kinetic_wl2_rows)
    n        <- min(nrow(tdat_wl1), nrow(tdat_wl2))
    as.matrix(tdat_wl1[seq_len(n), ] / tdat_wl2[seq_len(n), ])
  } else {
    as.matrix(tdat_wl1)
  }
  
  # Drop time and temperature columns (first two)
  tdat <- tdat[, seq(3, ncol(tdat)), drop = FALSE]
  time <- as.numeric(run_df[kinetic_wl1_rows, time_col])[seq_len(nrow(tdat))]
  
  # ── Convert to long format ───────────────────────────────────────────────────
  long <- data.frame(
    Time  = rep(time, times = ncol(tdat)),
    Well  = rep(colnames(tdat), each = nrow(tdat)),
    Value = as.vector(tdat)
  )
  long <- subset(long, long$Well %in% wells_valid)
  
  # ── Normalise per well and append start/end rows ─────────────────────────────
  out <- do.call(rbind, by(long, long$Well, function(well_df) {
    t_step   <- well_df$Time[2] - well_df$Time[1]
    cur_well <- unique(well_df$Well)
    
    if (isTRUE(t_min_origin))
      well_df <- subset(well_df, well_df$Time >= well_df$Time[which.min(well_df$Value)])
    
    end_row <- data.frame(
      Time  = well_df$Time[nrow(well_df)] + t_step,
      Value = end_df$EndRatio[end_df$Well == cur_well],
      Well  = cur_well
    )
    
    tmp <- if (isTRUE(include_start)) {
      start_row <- data.frame(
        Time  = well_df$Time[1] - t_step,
        Value = start_df$StartRatio[start_df$Well == cur_well],
        Well  = cur_well
      )
      rbind(start_row, well_df, end_row)
    } else {
      rbind(well_df, end_row)
    }
    
    tmp$y           <- norm(tmp$Value) * 100
    tmp$y[is.na(tmp$y)] <- 100
    tmp$x           <- tmp$Time
    tmp
  }))
  
  # ── Merge plate map ──────────────────────────────────────────────────────────
  out <- merge(out, plate_map, by = "Well", all.x = TRUE, sort = FALSE)
  out <- data.frame(x = out$Time, raw = out$Value, y = out$y,
                    id = out$id, Well = out$Well)
  
  message("platR: done — ", nrow(out), " rows, ", length(unique(out$Well)),
          " well(s)")
  return(out)
}