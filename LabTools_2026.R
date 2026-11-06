
# General-purpose utilities for data import, processing, plotting, and
# quantification used across the JLC/Howorka lab.
#




.jlc_pkgs <- c(
  "ggplot2", "grid", "gridExtra", "cowplot", "zoo",
  "readxl", "scales", "openxlsx", "ggthemes", "ggside",
  "patchwork", "drc"
)
.not_installed <- .jlc_pkgs[!(.jlc_pkgs %in% installed.packages()[, "Package"])]
if (length(.not_installed)) install.packages(.not_installed)
invisible(lapply(.jlc_pkgs, require, character.only = TRUE))



#' Remove panel labels from facet plots
#'
#' Usage: `my_plot + strip_blank`
lab_rmv <- ggplot2::theme(
  strip.background   = ggplot2::element_blank(),
  strip.text.x       = ggplot2::element_blank(),
  strip.text.y       = ggplot2::element_blank(),
  strip.background.x = ggplot2::element_blank(),
  strip.background.y = ggplot2::element_blank()
)

#' Removes gridlines, backgrounds, borders, and axis decorations from
#' ggside marginal panels.
#' Usage: `my_plot + side_panel_clean`
side_panel_clean <- ggplot2::theme(
  ggside.panel.grid       = ggplot2::element_blank(),
  ggside.panel.background = ggplot2::element_blank(),
  ggside.panel.border     = ggplot2::element_blank(),
  ggside.axis.text        = ggplot2::element_blank(),
  ggside.axis.ticks       = ggplot2::element_blank()
)

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


#' Rolling mean and SD across multiple y columns of a data frame
#'
#' Computes a centred rolling mean and SD for each y column.
#'
#' The smoothing window can be specified in physical time using
#' `smooth_width`, in ns. The function automatically determines the
#' timestep from the spacing between rows and converts the requested
#' smoothing width into an appropriate number of data points.
#'
#' The x column is assumed to contain time in ns. If the sampling
#' timestep cannot or should not be inferred from x, it can be supplied
#' manually using `timestep_ps`, in ps.
#'
#' For example, if the trajectory has a timestep of 0.1 ns,
#' `smooth_width = 5` corresponds to approximately 51 data points.
#' If the trajectory is sampled every 1 ps, the same 5 ns smoothing
#' window corresponds to approximately 5001 data points.
#'
#' If `smooth_width = NULL`, the legacy `block` argument is used instead.
#'
#' Any non-numeric column other than x_col is auto-detected as a grouping
#' column (e.g. "id", "region", or both). The data is split by the unique
#' combination of those columns, and the rolling mean is applied within
#' each group.
#'
#' @param input_df data.frame Input containing an x column, one or more
#'   y columns, and optionally one or more grouping columns.
#' @param block integer Legacy rolling window size in data points.
#'   Only used when smooth_width = NULL. Default: 5.
#' @param x_col character or integer Name or index of the x column.
#' @param id character Label used to create an `id` column when no grouping
#'   column already produces one.
#' @param group_cols character vector Names of columns to group by.
#'   If NULL (default), auto-detected as every non-numeric column other
#'   than x_col.
#' @param timestep_ps numeric Optional manually supplied timestep in ps.
#'   If NULL (default), the timestep is automatically calculated from
#'   the spacing between x values. The x column is assumed to be in ns.
#' @param smooth_width numeric Width of smoothing window in ns.
#'   Default: 1 ns. Set to NULL to use the legacy `block` argument.
#'
#' @return data.frame with columns: x, y (rolling mean), stdev
#'   (rolling SD), variable (source y column name), any grouping columns,
#'   and id.
run_avg <- function(input_df,
                    block = 5,
                    x_col = 1,
                    id = "id_missing",
                    group_cols = NULL,
                    timestep_ps = NULL,
                    smooth_width = 1) {

  
  if (is.numeric(x_col)) {
    x_pos <- x_col
    x_name <- names(input_df)[x_pos]
  } else {
    x_name <- x_col
    if (!(x_name %in% names(input_df))) {
      stop(
        "run_avg: x_col '",
        x_name,
        "' not found in input_df."
      )
    }
  }

  
  nm <- names(input_df)
  if (is.null(nm)) {
    nm <- rep("", ncol(input_df))
  }
  missing_nm <- is.na(nm) | nm == ""
  if (any(missing_nm)) {
    nm[missing_nm] <- paste0(
      "y",
      seq_len(sum(missing_nm))
    )
    names(input_df) <- make.unique(nm)
    
    message(
      "run_avg: ",
      sum(missing_nm),
      " column(s) had missing names; auto-assigned ",
      paste(nm[missing_nm], collapse = ", "),
      "."
    )
  }
  

  if (is.null(group_cols)) {
    
    group_cols <- names(input_df)[
      names(input_df) != x_name &
        !sapply(input_df, is.numeric)
    ]
  } else {
    
    missing_groups <- setdiff(
      group_cols,
      names(input_df)
    )
    if (length(missing_groups) > 0) {
      stop(
        "run_avg: group_cols not found in input_df: ",
        paste(missing_groups, collapse = ", ")
      )
    }
  }
  
  .compute <- function(df) {

    y_cols <- setdiff(
      names(df),
      c(x_name, group_cols)
    )
    y_cols <- y_cols[
      sapply(
        df[y_cols],
        is.numeric
      )
    ]
    if (length(y_cols) == 0) {
      stop(
        "run_avg: no numeric y columns found besides x_col '",
        x_name,
        "'."
      )
    }

    x_vals <- df[[x_name]]

    if (!is.null(smooth_width)) {
      dx <- median(
        diff(x_vals),
        na.rm = TRUE
      )
      if (!is.finite(dx) || dx <= 0) {
        stop(
          "run_avg: could not determine a valid timestep from x_col."
        )
      }
      
      if (is.null(timestep_ps)) {
        
        timestep_ns <- dx
        message(
          "run_avg: timestep automatically detected as ",
          signif(timestep_ns, 6),
          " ns."
        )
      } else {
        
        if (!is.numeric(timestep_ps) ||
            length(timestep_ps) != 1 ||
            !is.finite(timestep_ps) ||
            timestep_ps <= 0) {
          stop(
            "run_avg: timestep_ps must be a single positive numeric value."
          )
        }
        timestep_ns <- timestep_ps / 1000
        message(
          "run_avg: timestep manually supplied as ",
          signif(timestep_ps, 6),
          " ps (",
          signif(timestep_ns, 6),
          " ns)."
        )
      }
      block <- round(
        smooth_width / timestep_ns
      )
      if (block %% 2 == 0) {
        block <- block + 1
      }
      block <- max(
        block,
        3
      )
     message(
        "run_avg: smoothing = ",
        smooth_width,
        " ns; block = ",
        block,
        " points."
      )
    }
    
    result_list <- lapply(
      y_cols,
      function(y_name) {
        
        y_vals <- df[[y_name]]
        rm_vals <- rollmean(
          y_vals,
          k = block,
          fill = NA,
          align = "center"
        )
        sd_vals <- rollapply(
          y_vals,
          width = block,
          FUN = sd,
          fill = NA,
          align = "center"
        )
        
        out <- data.frame(
          x        = x_vals,
          y        = rm_vals,
          stdev    = sd_vals,
          variable = y_name,
          stringsAsFactors = FALSE
        )
        
        out
      }
    )
    
    out_df <- do.call(
      rbind,
      result_list
    )
    rownames(out_df) <- NULL
    out_df
  }
  
 if (length(group_cols) > 0) {
    
    grp_key <- do.call(
      interaction,
      c(
        as.list(input_df[group_cols]),
        list(
          drop = TRUE,
          sep = "_"
        )
      )
    )
  
    results <- lapply(
      split(input_df, grp_key),
      function(sub_df) {
        
        if (nrow(sub_df) == 0) {
          return(NULL)
        }
        
        out <- .compute(sub_df)
        
        
        # Carry grouping variables through
        for (gc in group_cols) {
          out[[gc]] <- sub_df[[gc]][1]
        }
        
        out
      }
    )
    out_df <- do.call(
      rbind,
      results
    )
    rownames(out_df) <- NULL
  } else {
    
    out_df <- .compute(
      input_df
    )
  }
  
  
  # ==========================================================
  # Backwards compatibility
  # ==========================================================
  
  if (!("id" %in% names(out_df))) {
    out_df$id <- id
  }
  
  out_df
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


#' Read a 96-well plate fluorescence file (xlsx) into a data frame
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


#' General-purpose xy scatter / line plot
#'
#' Builds a ggplot2 line or point plot from a data frame. Supports error
#' ribbons or bars, manual colours, faceting, log axes, and axis limits.
#'
#' @param input_df    data.frame containing columns named by xvar and yvar.
#' @param xvar        Character. X-axis column (default: "x").
#' @param yvar        Character. Y-axis column (default: "y").
#' @param xlab        Character/NULL. X-axis label; defaults to xvar.
#' @param ylab        Character/NULL. Y-axis label; defaults to yvar.
#' @param colour_by   Character/NULL. Column used for colour/fill grouping.
#'                    NULL plots everything in a single colour.
#' @param leg_loc     Character. Legend position (default: "right").
#' @param error       Logical or character. TRUE/"full" enables full SD;
#'                    "half" enables half SD; FALSE/"none" disables errors.
#'                    Also accepts "y"/"yes" and "n"/"no".
#' @param error_type  Character. "ribbon" (default) or "errorbar".
#' @param cols        Character/NULL. Manual colour values.
#' @param facet       Character/NULL. Column name to facet by (rows ~ facet).
#' @param marker      Character. "line" (default) or "point".
#' @param marker_size Numeric. Geom size (default: 1).
#' @param n_ticks     Integer. Approximate number of x-axis breaks.
#' @param log         Character. Log-scale axis: "x", "y", or "none".
#' @param sci_notation Logical. TRUE enables scientific axis notation.
#' @param leg_title   Character/NULL. Legend title; NULL removes the title.
#' @param plot_ylim   Numeric/NULL. Y-axis limits c(min, max).
#' @param plot_xlim   Numeric/NULL. X-axis limits c(min, max).
#' @param txt_size    Numeric. Base font size (default: 12).
#'
#' @return ggplot2 object.
plotr <- function(input_df,
                  xvar         = "x",
                  yvar         = "y",
                  xlab         = NULL,
                  ylab         = NULL,
                  colour_by    = "id",
                  leg_loc      = "right",
                  error        = FALSE,
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
  
  # ── Basic validation ──────────────────────────────────────────────────────
  
  if (!is.data.frame(input_df)) {
    stop("input_df must be a data.frame.")
  }
  
  if (!all(c(xvar, yvar) %in% names(input_df))) {
    stop("input_df must contain columns named '", xvar, "' and '", yvar, "'.")
  }
  
  if (!is.null(colour_by) && !colour_by %in% names(input_df)) {
    stop("colour_by column '", colour_by, "' not found in input_df.")
  }
  
  if (!is.null(facet) && !facet %in% names(input_df)) {
    stop("facet column '", facet, "' not found in input_df.")
  }
  
  if (is.null(xlab)) xlab <- xvar
  if (is.null(ylab)) ylab <- yvar
  
  marker <- tolower(marker)
  error_type <- tolower(error_type)
  log <- tolower(log)
  
  if (!marker %in% c("line", "point")) {
    stop("marker must be 'line' or 'point'.")
  }
  
  if (!error_type %in% c("ribbon", "errorbar")) {
    stop("error_type must be 'ribbon' or 'errorbar'.")
  }
  
  if (!log %in% c("x", "y", "none")) {
    stop("log must be 'x', 'y', or 'none'.")
  }
  
  if (!is.logical(sci_notation) ||
      length(sci_notation) != 1L ||
      is.na(sci_notation)) {
    stop("sci_notation must be TRUE or FALSE.")
  }
  
  # ── Interpret error setting ───────────────────────────────────────────────
  
  if (is.logical(error) && length(error) == 1L && !is.na(error)) {
    
    err_mode <- if (error) "full" else "none"
    
  } else if (is.character(error) && length(error) == 1L &&
             !is.na(error)) {
    
    err_mode <- tolower(error)
    
    # Backwards-compatible aliases
    if (err_mode %in% c("y", "yes", "true")) {
      err_mode <- "full"
    }
    
    if (err_mode %in% c("n", "no", "false")) {
      err_mode <- "none"
    }
    
  } else {
    stop("error must be TRUE/FALSE or a valid character option.")
  }
  
  if (!err_mode %in% c("full", "half", "none")) {
    stop("error must be TRUE, FALSE, 'full', 'half', or 'none'.")
  }
  
  # ── Error data ────────────────────────────────────────────────────────────
  
  # Work on a local copy; don't modify the original data frame.
  plot_df <- input_df
  
  if (err_mode != "none") {
    
    if (!"stdev" %in% names(plot_df)) {
      stop("An 'stdev' column is required when error is enabled.")
    }
    
    err_factor <- if (err_mode == "half") 0.5 else 1
    
    plot_df$ymax <- plot_df[[yvar]] + err_factor * plot_df$stdev
    plot_df$ymin <- plot_df[[yvar]] - err_factor * plot_df$stdev
  }
  

  # ── Base aesthetics ──────────────────────────────────────────────────────
  
  base_aes <- if (is.null(colour_by)) {
    
    ggplot2::aes(
      x = .data[[xvar]],
      y = .data[[yvar]]
    )
    
  } else {
    
    ggplot2::aes(
      x      = .data[[xvar]],
      y      = .data[[yvar]],
      colour = .data[[colour_by]],
      fill   = .data[[colour_by]]
    )
  }
  # ── Base plot ─────────────────────────────────────────────────────────────
  
  plt <- ggplot2::ggplot(plot_df, base_aes) +
    ggplot2::theme_bw(base_size = txt_size) +
    ggplot2::theme(
      text            = ggplot2::element_text(size = txt_size),
      legend.position = leg_loc
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    My_Theme(txt_size = txt_size)
  
  # ── Error display ─────────────────────────────────────────────────────────
  
  if (err_mode != "none") {
    
    err_aes <- ggplot2::aes_string(
      ymin = "ymin",
      ymax = "ymax"
    )
    
    if (error_type == "ribbon") {
      
      plt <- plt +
        ggplot2::geom_ribbon(
          mapping = err_aes,
          alpha   = 0.25,
          colour  = NA
        )
      
    } else {
      
      plt <- plt +
        ggplot2::geom_errorbar(
          mapping = err_aes,
          alpha   = 0.5,
          colour  = "black",
          width   = marker_size / 5
        )
    }
  }
  
  # ── Marker ────────────────────────────────────────────────────────────────
  
  if (marker == "line") {
    
    plt <- plt +
      ggplot2::geom_path(size = marker_size)
    
  } else {
    
    plt <- plt +
      ggplot2::geom_point(size = marker_size)
  }
  
  # ── Colours ───────────────────────────────────────────────────────────────
  
  if (!is.null(cols) && !is.null(colour_by)) {
    
    plt <- plt +
      ggplot2::scale_colour_manual(
        name   = leg_title,
        values = cols,
        drop   = TRUE
      ) +
      ggplot2::scale_fill_manual(
        name   = leg_title,
        values = cols,
        drop   = TRUE
      )
  }
  
  # ── Legend ────────────────────────────────────────────────────────────────
  
  if (!is.null(colour_by)) {
    
    plt <- plt +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          override.aes = list(size = 10, shape = 20)
        ),
        fill = "none"
      )
    
    if (is.null(leg_title)) {
      plt <- plt +
        ggplot2::theme(legend.title = ggplot2::element_blank())
    }
  }
  
  # ── Facet ─────────────────────────────────────────────────────────────────
  
  if (!is.null(facet)) {
    plt <- plt +
      ggplot2::facet_grid(stats::reformulate(facet, "."))
  }
  
  # ── Axis scales ───────────────────────────────────────────────────────────
  
  make_scale <- function(axis, is_log, lims, use_sci) {
    
    scale_fn <- if (is_log) {
      
      if (axis == "x") {
        ggplot2::scale_x_log10
      } else {
        ggplot2::scale_y_log10
      }
      
    } else {
      
      if (axis == "x") {
        ggplot2::scale_x_continuous
      } else {
        ggplot2::scale_y_continuous
      }
    }
    
    args <- list()
    
    if (!is_log && !is.null(lims)) {
      
      sb <- stefBreaks(lims)
      
      args$breaks <- sb[[1]]
      args$labels <- sb[[2]]
      
    } else if (!use_sci) {
      
      args$labels <- function(v) {
        format(v, scientific = FALSE, trim = TRUE)
      }
    }
    
    do.call(scale_fn, args)
  }
  
  plt <- plt +
    make_scale("x", log == "x", plot_xlim, sci_notation) +
    make_scale("y", log == "y", plot_ylim, sci_notation)
  
  # ── Axis limits ───────────────────────────────────────────────────────────
  
  if (!is.null(plot_xlim) || !is.null(plot_ylim)) {
    
    plt <- plt +
      ggplot2::coord_cartesian(
        xlim = plot_xlim,
        ylim = plot_ylim
      )
  }
  
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
    input_df$ymin <- input_df$y - input_df$stdev
  } else if (err_mode %in% c("half", "0.5")) {
    input_df$ymax <- input_df$y + 0.5 * input_df$stdev
    input_df$ymin <- input_df$y - 0.5 * input_df$stdev
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
#' 

#' Detect local maxima in a numeric vector
#'
#' Identifies peaks by locating sign changes in the first difference.
#' Optional smoothing reduces noise sensitivity, and a prominence threshold
#' filters out weak peaks.
#'
#' @param y               numeric Input signal vector.
#' @param x               numeric Optional x-coordinates. If supplied, data
#'                          are sorted by x before peak detection.
#' @param smooth          integer Moving-average window width. 0 = no smoothing.
#' @param min_prominence  numeric Minimum peak height relative to the global
#'                          signal range, on a [0, 1] scale.
#'
#' @return Data frame with columns: index, x, value, smoothed_value, prominence.
#'   Returns NULL if no peaks pass the threshold.
# install.packages("pracma") # if needed

peakID <- function(y, x = NULL, smooth = 0,
                   min_prominence = 0.01) {
  
  if (!is.null(x)) {
    if (length(x) != length(y)) {
      stop("x and y must have the same length.")
    }
    
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
  }
  
  if (smooth > 1) {
    kernel <- rep(1 / smooth, smooth)
    y_s <- as.numeric(stats::filter(y, kernel, sides = 2))
    y_s[is.na(y_s)] <- y[is.na(y_s)]
  } else {
    y_s <- y
  }
  
  # Detect local maxima
  peak_idx <- which(diff(sign(diff(y_s))) == -2) + 1
  
  if (length(peak_idx) == 0) {
    warning("peakID: no peaks detected.")
    return(NULL)
  }
  
  # Calculate local prominence for each peak
  prom <- vapply(peak_idx, function(i) {
    left <- i - 1
    while (left > 1 && y_s[left - 1] < y_s[i]) {
      left <- left - 1
    }
    
    right <- i + 1
    while (right < length(y_s) && y_s[right + 1] < y_s[i]) {
      right <- right + 1
    }
    
    min(y_s[i] - min(y_s[left:i]),
        y_s[i] - min(y_s[i:right]))
  }, numeric(1))
  
  # Normalise prominence to global signal range
  signal_range <- diff(range(y_s))
  
  if (signal_range == 0) return(NULL)
  
  prom_norm <- prom / signal_range
  keep <- prom_norm >= min_prominence
  
  if (!any(keep)) return(NULL)
  
  result <- data.frame(
    index = peak_idx[keep],
    value = y[peak_idx[keep]],
    smoothed_value = y_s[peak_idx[keep]],
    prominence = prom_norm[keep]
  )
  
  if (!is.null(x)) {
    result$x <- x[peak_idx[keep]]
    result <- result[, c("index", "x", "value",
                         "smoothed_value", "prominence")]
  }
  
  result
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


#' Compute the area under each peak in an xy signal
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


#' Prepare cuvette efflux data for plotting
#'
#' Removes pre-injection data, resets time to zero, normalises each
#' replicate, aligns replicates to a common time grid, and averages.
#'
#' @param input_df    Data frame containing x, y, sample, and id.
#' @param inject_list Injection times, recycled across replicates if needed.
#' @param normalise   "yes"/"y", "zero"/"0", "min", or "none".
#' @param time_round  Decimal places for the common time grid.
#'
#' @return Data frame with columns id, x, y, stdev.

effluxDataPrep <- function(input_df,
                           inject_list,
                           normalise = "yes",
                           time_round = 1) {
  
  # Injection time mapping: source file ID -> injection time
  source_ids <- unique(as.character(input_df$id))
  
  if (length(inject_list) != length(source_ids)) {
    stop(
      "Supply one injection time per source file ID. IDs are: ",
      paste(source_ids, collapse = ", ")
    )
  }
  
  injection_map <- setNames(inject_list, source_ids)
  
  # Remove only fluorimeter replicate suffixes
  input_df$condition <- sub(
    "\\.[0-9]+$", "", as.character(input_df$sample)
  )
  
  # Keep original recording name to distinguish replicates
  input_df$recording <- as.character(input_df$sample)
  
  groups <- split(input_df, input_df$recording)
  
  prepared <- lapply(groups, function(d) {
    
    t_inject <- injection_map[as.character(d$id[1])]
    
    d <- d[order(d$x), ]
    d <- d[d$x > t_inject, ]
    
    if (nrow(d) < 2) return(NULL)
    
    # Reset time independently for each recording
    d$x <- d$x - d$x[1]
    
    # Normalise each recording independently
    nrm <- tolower(normalise)
    
    if (nrm %in% c("yes", "y")) {
      yrange <- diff(range(d$y, na.rm = TRUE))
      if (!is.finite(yrange) || yrange == 0) return(NULL)
      
      d$y <- (d$y - min(d$y, na.rm = TRUE)) / yrange * 100
      
    } else if (nrm %in% c("zero", "0")) {
      d$y <- d$y - d$y[1]
      
    } else if (nrm == "min") {
      d$y <- d$y - min(d$y, na.rm = TRUE)
    }
    
    d <- d[is.finite(d$x) & is.finite(d$y), ]
    
    # Handle duplicate times within a recording
    d <- aggregate(y ~ x, data = d, FUN = mean)
    
    d$recording <- as.character(groups[[1]]$recording[1])
    d
  })
  
  # Rebuild recording metadata robustly
  prepared <- Filter(Negate(is.null), prepared)
  
  if (!length(prepared)) stop("No valid recordings after injection.")
  
  # Assign condition from the recording name
  for (i in seq_along(prepared)) {
    prepared[[i]]$recording <- names(prepared)[i]
    prepared[[i]]$condition <- sub(
      "\\.[0-9]+$", "", names(prepared)[i]
    )
  }
  
  # Common grid across recordings, in the same time units as x
  dt <- 1 / (10^time_round)
  max_time <- max(vapply(prepared, function(d) max(d$x), numeric(1)))
  grid <- seq(0, floor(max_time / dt) * dt, by = dt)
  
  # Interpolate within each recording's measured range
  aligned <- do.call(rbind, lapply(prepared, function(d) {
    
    grid_use <- grid[
      grid >= min(d$x) & grid <= max(d$x)
    ]
    
    data.frame(
      condition = d$condition[1],
      recording = d$recording[1],
      x = grid_use,
      y = approx(d$x, d$y, xout = grid_use, ties = mean)$y
    )
  }))
  
  # Average across recordings in each condition
  out <- aggregate(
    y ~ condition + x,
    data = aligned,
    FUN = function(v) c(
      mean = mean(v),
      sd = if (length(v) > 1) sd(v) else NA_real_,
      n = length(v)
    )
  )
  
  out <- do.call(data.frame, out)
  names(out) <- c("id", "x", "y", "stdev", "n")
  
  out
}

#' Background-corrected maximum efflux values for pairwise sample/background
#'
#' Extracts maximum y values up to t_max` for four groups (bkg1, sample1,
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
  dib_path  <- "~/Onedrive/Lab_Data/UV-VIS/SCANS/20220407_New_ND_Quant/"
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
