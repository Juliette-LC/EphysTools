# =============================================================================
# Nanopore Electrophysiology Analysis Toolkit
# =============================================================================
# Functions for importing, processing, and visualising current recordings
# from Oxford Nanopore (FAST5/HDF5) and OrbitMini (.abf) devices.
#
# Conventions used throughout:
#   - Data frames passed as the first argument are always named `input_df`
#   - Boolean-like parameters use TRUE/FALSE or NULL (not "y"/"n" strings)
#   - File paths are constructed with file.path() throughout
#   - message() is used for progress reporting; warning() for recoverable issues
# =============================================================================

# -------------- Dependencies ----------------------------#
## install packages if required
pkgs <- c("rhdf5","rhdf5filters","Rhdf5lib", "cowplot", "ggpubr", "scales", "signal","readABF",
          "data.table", "zoo", "docstring", "readABF", "asdetect","fitdistrplus", "scattermore", 
          "svglite", "parallel", "MASS", "future.apply" )

not_installed <- pkgs[!(pkgs %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)

# ### will need bioconductor to install rhdf5, rhsdffilters
#  if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  BiocManager::install(version = "3.19")
# 
# BiocManager::install("rhdf5")## load in packages
lapply(pkgs, require, character.only=TRUE)

## Below should not be necessary as should be self-reliant now 
# needed for myTheme in plotting functions
#source("~/Onedrive/Scripts/R/jlc_lab_tools.R")
source("~/Onedrive/Scripts/R/LabTools_2026.R")




# ── Section 1: Data Import ────────────────────────────────────────────────────

#' Read OrbitMini electrophysiology recordings (.abf files)
#'
#' Reads one or more ABF files from a folder, optionally corrects DC offset,
#' crops to a time window, and returns a long-format data frame. A QC plot of
#' all channels is printed unless suppressed.
#'
#' @param direc           character  Absolute path of the data directory.
#' @param folder          character  Experiment sub-folder inside `direc`.
#' @param file_path       character/NULL  Single filename, character vector of
#'                          filenames, or NULL to read every .abf in `folder`.
#' @param id              character/NULL  Label appended as an `id` column.
#' @param downsample      numeric/NULL  Chunk factor for downsampling the
#'                          returned data. NULL skips downsampling.
#' @param offset_fix      character  "robust" → offsetCorrect(); "hard" →
#'                          subtract the first sample (zero-normalise);
#'                          anything else → skip.
#' @param skip_plot       logical  TRUE suppresses the QC plot.
#' @param plot_sec_axis   logical  TRUE adds a secondary voltage axis on the
#'                          QC plot.
#' @param plot_downsample numeric  Chunk factor for the QC plot only; does not
#'                          affect the returned data.
#' @param write_csv       logical  TRUE writes a CSV alongside each .abf file.
#' @param time_subset     numeric vector/NULL  c(t_start, t_end) crops the
#'                          trace to this window.
#' @param convert_to_pA   logical  TRUE multiplies current × 1000 (nA → pA).
#'
#' @return data.frame with columns: time (s), current (pA or nA), voltage (mV),
#'   cm (channel label), id (experiment label).
readOrbit <- function(direc,
                      folder,
                      file_path       = NULL,
                      id              = NULL,
                      downsample      = NULL,
                      offset_fix      = "none",
                      skip_plot       = FALSE,
                      plot_sec_axis   = FALSE,
                      plot_downsample = 20007,
                      write_csv       = FALSE,
                      time_subset     = NULL,
                      convert_to_pA   = TRUE) {
  
  data_path <- file.path(direc, folder)
  
  # ── Internal helper: load and pre-process one ABF file ──────────────────────
  read_one_abf <- function(abf_file) {
    full_path   <- file.path(data_path, abf_file)
    file_stem   <- sub("\\..*$", "", abf_file)
    csv_out     <- file.path(data_path, paste0(file_stem, "-analysed.csv"))
    
    message("  Reading: ", full_path)
    
    dat        <- readABF(full_path)
    dat        <- as.data.frame(dat, sweep = 1)
    names(dat) <- c("time", "current", "voltage")
    
    if (isTRUE(write_csv)) {
      message("  Writing CSV: ", csv_out)
      fwrite(dat, csv_out, row.names = FALSE)
    }
    
    if (!is.null(downsample) && is.numeric(downsample))
      dat <- chunker(dat, chunk = downsample)
    
    if (offset_fix == "robust")
      dat <- offsetCorrect(dat)
    if (offset_fix == "hard")
      dat$current <- dat$current - dat$current[1]
    
    dat        <- na.omit(dat)
    dat$id     <- file_stem
    dat
  }
  
  # ── Resolve which files to read ─────────────────────────────────────────────
  if (!is.null(file_path) && length(file_path) == 1) {
    message("readOrbit: reading single file — ", file_path)
    return(read_one_abf(file_path))
  }
  
  dat_files <- if (length(file_path) > 1) {
    message("readOrbit: reading ", length(file_path), " specified files")
    file_path
  } else {
    found <- grep("\\.abf$", list.files(data_path), value = TRUE)
    message("readOrbit: found ", length(found), " .abf files in ", data_path)
    found
  }
  
  df <- rbindlist(
    lapply(dat_files, read_one_abf),
    use.names = TRUE, fill = TRUE
  )
  
  # ── Stitch multi-part recordings split by the device (> 10 min) ─────────────
  # The device appends a 3-digit numeric suffix when it splits a recording.
  # Detected by comparing unique id count vs unique channel (cm) count.
  df$cm <- sub(".*_(.*)_.*", "\\1", df$id)
  
  if (length(unique(df$id)) > length(unique(df$cm))) {
    message("readOrbit: detected split recordings — stitching time axes")
    df$part <- substr(df$id, nchar(df$id) - 2, nchar(df$id))
    
    df <- rbindlist(lapply(split(df, df$cm), function(sdf) {
      sdf       <- sdf[order(sdf$part, sdf$time), ]
      max_times <- tapply(sdf$time, sdf$part, max)
      offsets   <- c(0, cumsum(max_times[-length(max_times)]))
      names(offsets) <- names(max_times)
      
      for (part_id in names(offsets))
        sdf$time[sdf$part == part_id] <-
        sdf$time[sdf$part == part_id] + offsets[part_id]
      
      sdf[, c("time", "current", "voltage", "id")]
    }), use.names = TRUE)
  }
  
  # ── Optional time crop ───────────────────────────────────────────────────────
  if (!is.null(time_subset)) {
    message("readOrbit: cropping to t = [", time_subset[1], ", ", time_subset[2], "] s")
    df <- subset(df, df$time > time_subset[1] & df$time < time_subset[2])
  }
  
  # ── QC plot ──────────────────────────────────────────────────────────────────
  if (!isTRUE(skip_plot)) {
    message("readOrbit: generating QC plot")
    tmp <- as.data.frame(df)
    if (is.numeric(plot_downsample))
      tmp <- rbindlist(lapply(split(tmp, tmp$id), chunker, chunk = plot_downsample))
    print(
      currPlot(tmp, rel_time = FALSE, sec_axis = plot_sec_axis, facet = "l") +
        ggplot2::ggtitle(folder)
    )
  }
  
  # ── Convert nA → pA ──────────────────────────────────────────────────────────
  if (isTRUE(convert_to_pA)) {
    message("readOrbit: converting current nA → pA")
    df$current <- df$current * 1000
  }
  
  rownames(df) <- NULL
  message("readOrbit: done — ", nrow(df), " rows, ", length(unique(df$id)), " channel(s)")
  return(as.data.frame(df))
}


#' Import raw nanopore signal from a FAST5/HDF5 file
#'
#' Scales the raw ADC signal for a single channel to physical units (pA / mV),
#' optionally caching the result as a CSV for faster re-reads.
#'
#' @param direc       character  Base working directory.
#' @param folder      character  Sub-folder containing the FAST5 file.
#' @param file_name   character  FAST5 filename.
#' @param channel     integer    Channel number to extract.
#' @param bias_vol    numeric    Voltage scale factor (default: −5).
#' @param save_csv    logical    TRUE writes a cached CSV alongside the FAST5.
#' @param re_extract  logical    TRUE forces re-extraction even if a CSV exists.
#' @param ...         Reserved for future arguments.
#'
#' @return data.table with columns: current (pA), time (s), voltage (mV),
#'   channel (int).
rawImport <- function(direc,
                      folder,
                      file_name,
                      channel,
                      bias_vol    = -5,
                      save_csv    = FALSE,
                      re_extract  = FALSE,
                      ...) {
  
  channel   <- as.integer(channel)
  full_path <- file.path(direc, folder, file_name)
  
  # Derive a short run ID (last 5 chars before extension) for cache naming
  run_id   <- substr(substr(file_name, nchar(file_name) - 10, nchar(file_name)), 1, 5)
  cache_dir <- file.path(direc, folder, paste0(run_id, "_CSV"))
  csv_path  <- file.path(cache_dir, paste0("run_", run_id, "-channel_", channel, ".csv"))
  
  # ── Return cached CSV if available ──────────────────────────────────────────
  if (!isTRUE(re_extract) && file.exists(csv_path)) {
    message("  rawImport: reading cached CSV for channel ", channel)
    dt <- data.table::fread(csv_path)
    if ("channel" %in% names(dt))
      dt[, channel := as.integer(channel)]
    return(dt)
  }
  
  # ── Extract and scale from FAST5 ────────────────────────────────────────────
  message("  rawImport: extracting raw signal for channel ", channel)
  rhdf5::h5closeAll()
  
  raw  <- rhdf5::h5read(full_path, paste0("/Raw/Channel_", channel, "/Signal"))
  meta <- rhdf5::h5readAttributes(
    full_path, paste0("/IntermediateData/Channel_", channel, "/Meta/")
  )
  
  scale   <- meta$range / meta$digitisation
  current <- (raw + meta$offset) * scale
  time    <- seq(0, length(current) - 1) / meta$sample_rate
  
  voltage_meta <- rhdf5::h5read(full_path, "/Device/MetaData")
  voltage      <- voltage_meta$bias_voltage * bias_vol
  
  dt <- data.table::data.table(
    current = as.numeric(current),
    time    = as.numeric(time),
    voltage = as.numeric(voltage),
    channel = as.integer(channel)
  )
  
  if (isTRUE(save_csv)) {
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    data.table::fwrite(dt, csv_path)
    message("  rawImport: saved CSV → ", csv_path)
  }
  
  rhdf5::h5closeAll()
  return(dt)
}


#' Extract event-level data for a single channel from a FAST5 file
#'
#' Reads the pre-computed event table (mean current per event) from the
#' IntermediateData group. Used internally by runReporter.
#'
#' @param channel   integer  Channel number.
#' @param data_path character  Full path to the FAST5 file.
#'
#' @return data.frame with columns: time (s), current (pA), channel.
channelExtract <- function(channel, data_path) {
  
  rhdf5::h5closeAll()
  
  events <- data.frame(
    rhdf5::h5read(data_path, paste0("/IntermediateData/Channel_", channel, "/Events"))
  )
  meta <- data.frame(
    rhdf5::h5readAttributes(data_path, paste0("/IntermediateData/Channel_", channel, "/Meta"))
  )
  
  sample_rate <- meta$sample_rate
  
  df          <- data.frame(
    time    = events$start / sample_rate,
    current = events$mean,
    channel = channel
  )
  
  return(df)
}

# Keep the old name as an alias so existing scripts don't break
channelExtrct <- channelExtract


# ── Section 2: Warm-up & Mux Handling ────────────────────────────────────────

#' Detect chip warm-up period and mux-change timepoints
#'
#' Reads the raw signal from one channel to locate voltage transitions,
#' identifies the repeating voltage programme used for mux cycling, and
#' returns the time at which the warm-up ends plus the times of each mux
#' change.
#'
#' @param full_path           character  Full path to the FAST5 file.
#' @param return_mux_times    logical/NULL  Non-NULL returns mux change times
#'                              as the second list element.
#' @param channel             integer  Channel to use for voltage detection
#'                              (default: 1).
#' @param bias_vol            numeric  Voltage scale factor (default: −5).
#' @param return_voltage_steps logical/NULL  Non-NULL adds a third list element:
#'                              a data frame of (time, voltage) at each step.
#' @param mux_count           integer  Number of muxes on the device (default: 4
#'                              for MinION). Reduce if the run crashed early.
#' @param voltage_program     numeric/NULL  Override automatic programme
#'                              detection by supplying the voltage sequence.
#' @param ...                 Additional args forwarded to h5read / h5readAttributes.
#'
#' @return list:
#'   [[1]] numeric  Warm-up end time (s) — crop data to `time > [[1]]`.
#'   [[2]] numeric vector  Mux-change times (s).  (if return_mux_times != NULL)
#'   [[3]] data.frame (time, voltage) at each step. (if return_voltage_steps != NULL)
warmUpDetector <- function(full_path,
                           return_mux_times    = NULL,
                           channel             = 1,
                           bias_vol            = -5,
                           return_voltage_steps = NULL,
                           mux_count           = 4,
                           voltage_program     = NULL,
                           ...) {
  
  dots <- list(...)
  
  # Helper: call a function forwarding only args it accepts
  safe_call <- function(fun, args) {
    do.call(fun, c(args, dots[names(dots) %in% names(formals(fun))]))
  }
  
  # ── Load and scale raw signal ────────────────────────────────────────────────
  message("warmUpDetector: reading channel ", channel, " from ", basename(full_path))
  raw  <- data.frame(safe_call(rhdf5::h5read,
                               list(full_path, paste0("/Raw/Channel_", channel, "/Signal"))))
  meta <- data.frame(safe_call(rhdf5::h5readAttributes,
                               list(full_path, paste0("/IntermediateData/Channel_", channel, "/Meta/"))))
  
  scale      <- meta$range / meta$digitisation
  scaled_raw <- data.frame(current = (raw + meta$offset) * scale)
  sample_idx <- seq(0, nrow(scaled_raw) - 1)
  
  device_meta <- data.frame(safe_call(rhdf5::h5read, list(full_path, "/Device/MetaData")))
  
  raw_df <- data.frame(
    current = scaled_raw$current,
    time    = sample_idx / meta$sample_rate,
    voltage = device_meta$bias_voltage * bias_vol
  )
  
  # ── Locate voltage transitions ───────────────────────────────────────────────
  dt <- data.table::data.table(voltage = raw_df$voltage, index = seq_along(raw_df$voltage))
  dt[, is_first := (voltage != data.table::shift(voltage, fill = voltage[1] - 1)),
     by = data.table::rleid(voltage)]
  v_changes <- dt[is_first == TRUE, .(voltage, first_index = index)]
  
  if (nrow(v_changes) < 2 && is.null(voltage_program)) {
    warning("warmUpDetector: no voltage changes detected — cannot determine warm-up period.")
    return(list(0, numeric(0)))
  }
  
  # ── Determine (or accept) the voltage programme ──────────────────────────────
  v_seq <- v_changes$voltage
  
  if (is.null(voltage_program)) {
    message("warmUpDetector: auto-detecting voltage programme (mux_count = ", mux_count, ")")
    voltage_program <- find_pattern(v_seq, target_repeats = mux_count)
    message(
      "  Detected programme: ", paste(voltage_program, collapse = " "),
      "\n  Full sequence:      ", paste(v_seq, collapse = " "),
      "\n  If incorrect, set `voltage_program` manually or adjust `mux_count`."
    )
  }
  
  vp_len <- length(voltage_program)
  
  # ── Assign mux labels and locate the warm-up zero step ──────────────────────
  v_changes[, mux := (seq_len(.N) - 1) %/% vp_len + 1]
  
  zero_rows      <- which(v_changes$voltage == 0)
  zero_rows      <- zero_rows[-1]   # first zero is the warm-up; skip it
  zero_lengths   <- unlist(v_changes$first_index[zero_rows] -
                             unlist(v_changes$first_index[zero_rows - 1]))
  zero_step_len  <- mean(stats::na.omit(zero_lengths))
  
  init_v_change      <- v_changes[2, 2]
  warm_up_start_idx  <- init_v_change - zero_step_len
  
  # Prepend a mux-0 entry for the warm-up period, then adjust mux-1 start
  v_changes <- rbind(
    data.table::data.table(voltage = 0, first_index = 0, mux = 0),
    v_changes
  )
  v_changes$first_index[v_changes$mux == 1][1] <- warm_up_start_idx
  
  init_rmv        <- raw_df$time[unlist(v_changes[2, 2])]
  message("warmUpDetector: warm-up end time = ", round(init_rmv, 2), " s")
  
  # ── Mux change times ─────────────────────────────────────────────────────────
  mux_change_idx   <- v_changes$first_index[c(TRUE, diff(v_changes$mux) != 0)]
  mux_change_times <- raw_df$time[unlist(mux_change_idx)]
  
  if (!is.null(return_voltage_steps)) {
    voltage_steps_df <- raw_df[unlist(v_changes$first_index), c("time", "voltage")]
    return(list(init_rmv, mux_change_times, voltage_steps_df))
  }
  
  return(list(init_rmv, mux_change_times))
}


#' Split a channel-mux string vector into a two-column data frame
#'
#' Channel-mux identifiers encode the channel as all but the last character
#' and the mux as the final character (e.g. "1141" → channel 114, mux 1).
#'
#' @param cm_list  character vector  Channel-mux identifiers.
#'
#' @return data.frame with columns: channel (character), mux (character).
cmSplit <- function(cm_list) {
  data.frame(
    channel = vapply(cm_list, function(x) substr(x, 1, nchar(x) - 1), character(1)),
    mux     = vapply(cm_list, function(x) substr(x, nchar(x), nchar(x)), character(1)),
    stringsAsFactors = FALSE
  )
}


#' Assign mux labels and remove warm-up data from a current trace
#'
#' @param input_df  data.frame  Requires columns: current, time.
#' @param channel   integer     Channel number (written into the `channel` column).
#' @param init_rmv  list        Output of warmUpDetector:
#'                                [[1]] warm-up end time,
#'                                [[2]] mux-change times.
#'
#' @return data.frame with added columns: channel, mux (integer), cm (character).
muxSplitr <- function(input_df, channel, init_rmv) {
  
  input_df          <- subset(input_df, input_df$time > init_rmv[[1]])
  input_df$channel  <- channel
  input_df$mux      <- findInterval(input_df$time, stats::na.omit(init_rmv[[2]]))
  input_df$cm       <- paste0(input_df$channel, input_df$mux)
  
  return(input_df)
}


# ── Section 3: Run-Level QC ───────────────────────────────────────────────────

#' Generate per-channel current trace plots for a run report
#'
#' Loops over channels, extracts event-level data, removes warm-up, and
#' saves a multi-page PDF of current vs. time plots.
#'
#' @param direc           character  Base working directory.
#' @param folder          character  Sub-folder containing the FAST5 file.
#' @param file_name_list  character vector  FAST5 filenames to process.
#' @param channels        integer vector  Channels to include (default: 1–512).
#' @param plot_ylim       numeric vector/NULL  Y-axis limits, e.g. c(−200, 200).
#' @param out_file_name   character/NULL  Output PDF name stem. Default:
#'                          "run-report-<runID>.pdf".
#' @param high_res        logical  TRUE uses geom_line (slow); FALSE uses
#'                          geom_scattermore (fast, default).
#' @param ...             Forwarded to warmUpDetector.
runReporter <- function(direc,
                        folder,
                        file_name_list,
                        channels      = seq_len(512),
                        plot_ylim     = NULL,
                        out_file_name = NULL,
                        high_res      = FALSE,
                        ...) {
  
  lapply(file_name_list, function(file_name) {
    full_path <- file.path(direc, folder, file_name)
    run_id    <- substr(substr(file_name, nchar(file_name) - 10, nchar(file_name)), 1, 5)
    
    message("runReporter: processing ", file_name, " (run ", run_id, ")")
    init_rmv <- warmUpDetector(full_path, ...)
    
    plt_list <- lapply(channels, function(ch) {
      cdata <- channelExtract(ch, full_path)
      cdata <- subset(cdata, cdata$time > init_rmv[[1]])
      
      plt <- ggplot2::ggplot(cdata, ggplot2::aes(x = time, y = current)) +
        ggplot2::ylab("Current (pA)") +
        ggplot2::xlab("Time (s)") +
        My_Theme() +
        ggplot2::ggtitle(as.character(ch))
      
      plt <- if (isTRUE(high_res)) plt + ggplot2::geom_line()
      else                  plt + geom_scattermore()
      
      if (!is.null(plot_ylim))
        plt <- plt + ggplot2::ylim(plot_ylim)
      
      plt
    })
    
    out_file <- file.path(
      direc, folder,
      paste0(if (is.null(out_file_name)) "run-report" else out_file_name,
             "-", run_id, ".pdf")
    )
    
    message("runReporter: saving ", length(plt_list), " plots → ", out_file)
    arranged <- gridExtra::marrangeGrob(grobs = plt_list, nrow = 2, ncol = 2)
    tryCatch(
      ggplot2::ggsave(out_file, arranged, width = 29.7, height = 21, units = "cm"),
      error = function(e) warning("runReporter: ggsave failed — ", e$message)
    )
    
    rhdf5::h5closeAll()
  })
  
  invisible(NULL)
}


# ── Section 4: Signal Processing ─────────────────────────────────────────────

#' Correct DC offset in a current trace
#'
#' Estimates the open-pore baseline from the beginning of the trace and
#' subtracts it. A secondary check using a modal average over the mid-trace
#' region catches cases where the initial estimate fails.
#'
#' @param input_df    data.table  Requires columns: current, time, voltage.
#' @param method      character  "robust" (default) uses openPoreCalc on the
#'                     initial segment; "min" subtracts the minimum value.
#' @param ...         Forwarded to openPoreCalc.
#'
#' @return data.table with the current column offset-corrected in place.
offsetCorrect <- function(input_df, method = c("robust", "min"), ...) {
  
  method <- match.arg(method)
  data.table::setDT(input_df)
  input_df[, rel_time := time - time[1]]
  
  init_cur <- initSelector(input_df)
  
  if (method == "robust") {
    offset <- tryCatch(
      as.numeric(openPoreCalc(init_cur, ...)),
      error = function(e) NA_real_
    )
    
    if (!is.na(offset) && abs(offset) >= 0.005) {
      input_df[, current := current - offset]
      message("  offsetCorrect: primary offset = ", round(offset, 4))
    }
    
    # Secondary check over the 20–65 % time window
    offset_check <- modalAvg(
      input_df[rel_time > max(rel_time) * 0.2 &
                 rel_time < max(rel_time) * 0.65, current]
    )
    if (!is.na(offset_check) && abs(offset_check) > 0.05) {
      input_df[, current := current - offset_check]
      message("  offsetCorrect: secondary offset = ", round(offset_check, 4))
    }
    
  } else {
    offset <- min(init_cur, na.rm = TRUE)
    input_df[, current := current - offset]
    message("  offsetCorrect: min offset = ", round(offset, 4))
  }
  
  return(input_df)
}


#' Find open-pore current for each voltage step
#'
#' Calls openPoreCalc per voltage step and joins the result back as an `opc`
#' column. An adjustment factor nudges the estimate away from any partial
#' blockade events remaining after selection.
#'
#' @param input_df        data.table  Requires columns: voltage, current.
#' @param sd_adjust       numeric  Shifts each OPC estimate by this amount
#'                          (positive = towards zero; default: 1). Set to 0
#'                          to disable.
#' @param ...             Forwarded to openPoreCalc.
#'
#' @return data.table with an added `opc` column, or NULL if no valid OPC
#'   values could be calculated.
opcFinder <- function(input_df, sd_adjust = 1, ...) {
  
  if (!data.table::is.data.table(input_df))
    data.table::setDT(input_df)
  
  input_df[, v_step := data.table::rleid(voltage)]
  
  opc_tbl <- input_df[, .(opc = openPoreCalc(current, ...)), by = v_step]
  opc_tbl <- opc_tbl[!is.na(opc)]
  
  if (nrow(opc_tbl) == 0) {
    warning("opcFinder: no valid OPC values — returning NULL.")
    return(NULL)
  }
  
  if (!is.null(sd_adjust) && sd_adjust != 0) {
    opc_tbl[, opc := data.table::fifelse(opc > 0,
                                         opc - sd_adjust,
                                         opc + sd_adjust)]
  }
  
  message("  opcFinder: OPC estimated for ", nrow(opc_tbl), " voltage step(s)")
  opc_tbl[input_df, on = "v_step"]
}


#' Estimate the open-pore current from a current vector
#'
#' Trims transient edges, focuses on the upper quantile range (where the open
#' pore signal lives), and fits a robust linear intercept. Falls back to the
#' median if the robust fit fails.
#'
#' @param current_values  numeric  Current samples for one voltage step.
#' @param ...             Reserved; ignored.
#'
#' @return numeric scalar: estimated open-pore current.
openPoreCalc <- function(current_values, ...) {
  
  current_values <- current_values[is.finite(current_values)]
  
  if (length(current_values) < 10)
    return(stats::median(current_values, na.rm = TRUE))
  
  n       <- length(current_values)
  trimmed <- current_values[seq.int(n * 0.25, n * 0.9)]
  
  if (length(trimmed) < 10)
    return(stats::median(current_values, na.rm = TRUE))
  
  lower    <- stats::quantile(trimmed, 0.60, names = FALSE)
  upper    <- stats::quantile(trimmed, 0.95, names = FALSE)
  filtered <- trimmed[trimmed > lower & trimmed <= upper]
  
  if (length(filtered) < 10)
    return(stats::median(trimmed, na.rm = TRUE))
  
  result <- tryCatch(
    suppressMessages(stats::coef(MASS::rlm(filtered ~ 1, method = "MM"))[1]),
    error = function(e) stats::median(filtered, na.rm = TRUE)
  )
  
  as.numeric(result)
}


#' Apply a Butterworth low-pass filter to a current trace
#'
#' Designs a Butterworth filter via `signal::butter()` and applies it once.
#' The sampling rate is resolved from the argument, the HDF5 file metadata,
#' or estimated from the time column (in that order).
#'
#' @param input_df       data.table  Requires columns: current, time.
#' @param sample_rate    numeric/NULL  Sampling rate in Hz. NULL → auto-detect.
#' @param data_path      character/NULL  Path to HDF5 file for metadata lookup.
#' @param filter_order   integer  Filter order (default: 4).
#' @param filter_cutoff  numeric  Cutoff frequency in Hz (default: 10 000).
#' @param filter_type    character  One of "low", "high", "stop", "pass".
#' @param ...            Reserved; ignored.
#'
#' @return data.table with the current column replaced by the filtered signal.
#'   Leading/trailing 1 % edge artefacts are trimmed and time is reset to zero.
besselFilt <- function(input_df,
                       sample_rate   = NULL,
                       data_path     = NULL,
                       filter_order  = 4,
                       filter_cutoff = 10000,
                       filter_type   = c("low", "high", "stop", "pass"),
                       ...) {
  
  filter_type <- match.arg(filter_type)
  
  # ── Resolve sampling rate ────────────────────────────────────────────────────
  if (is.null(sample_rate) && !is.null(data_path) && file.exists(data_path)) {
    sample_rate <- tryCatch({
      meta <- rhdf5::h5readAttributes(data_path, "/IntermediateData/Channel_1/Meta/")
      meta$sample_rate
    }, error = function(e) {
      warning("besselFilt: could not read sample rate from file — ", e$message)
      NULL
    })
  }
  
  if (is.null(sample_rate) && "time" %in% names(input_df)) {
    diffs       <- diff(head(input_df$time, 100))
    sample_rate <- 1 / mean(diffs)
    message("  besselFilt: estimated sample rate = ", round(sample_rate), " Hz")
  }
  
  if (is.null(sample_rate)) {
    warning("besselFilt: sample rate unknown, defaulting to 4000 Hz")
    sample_rate <- 4000
  }
  
  if (!data.table::is.data.table(input_df))
    data.table::setDT(input_df)
  
  nyquist <- sample_rate / 2
  message("  besselFilt: sample rate = ", sample_rate,
          " Hz, Nyquist = ", nyquist, " Hz, cutoff = ", filter_cutoff, " Hz")
  
  if (filter_cutoff >= nyquist) {
    old <- filter_cutoff
    filter_cutoff <- nyquist * 0.9
    warning("besselFilt: cutoff (", old, " Hz) ≥ Nyquist; adjusted to ",
            round(filter_cutoff, 1), " Hz")
  }
  
  norm_cutoff <- filter_cutoff / nyquist
  if (norm_cutoff <= 0 || norm_cutoff >= 1)
    stop("besselFilt: invalid normalised cutoff: ", norm_cutoff)
  
  bw_filter <- signal::butter(filter_order, norm_cutoff, type = filter_type)
  input_df[, current := as.numeric(signal::filter(bw_filter, current))]
  
  # Trim 1 % edge artefacts
  n     <- nrow(input_df)
  start <- max(1L, as.integer(n * 0.01))
  end   <- min(n,  as.integer(n * 0.99))
  if (end > start) input_df <- input_df[start:end]
  
  input_df[, time := time - time[1]]
  return(input_df)
}


#' Extract the initial baseline segment of a current trace
#'
#' Selects data from the first contiguous voltage step, skips the first second
#' (voltage transition artefact), caps at `window_s` seconds, and removes
#' outliers via IQR filtering.
#'
#' @param input_df  data.frame/data.table  Requires: current, time, voltage.
#' @param window_s  numeric  Seconds of data to use (default: 10).
#' @param ...       Reserved; ignored.
#'
#' @return Numeric vector of baseline current samples, or numeric(0) if empty.
initSelector <- function(input_df, window_s = 10, ...) {
  
  input_df$rel_time  <- input_df$time - input_df$time[1]
  input_df$v_counter <- cumsum(c(TRUE,
                                 diff(as.numeric(as.character(input_df$voltage))) != 0))
  
  seg <- input_df$current[
    input_df$v_counter == 1 &
      input_df$rel_time > 1 &
      input_df$rel_time < window_s
  ]
  
  if (length(seg) == 0) {
    warning("initSelector: no data in initial segment — returning empty vector.")
    return(numeric(0))
  }
  
  q25 <- stats::quantile(seg, 0.25, na.rm = TRUE)
  q75 <- stats::quantile(seg, 0.75, na.rm = TRUE)
  iqr <- q75 - q25
  seg[seg > (q25 - 1.5 * iqr) & seg < (q75 + 1.5 * iqr)]
}


#' Detect the repeating pattern in a voltage sequence
#'
#' Finds the shortest sub-sequence that, when tiled, reproduces the full
#' vector. Used by warmUpDetector to identify the mux-cycling programme.
#'
#' @param x              numeric  Full voltage sequence.
#' @param target_repeats integer  Expected number of repetitions (default: 4).
#'
#' @return numeric vector: the repeating pattern.
find_pattern <- function(x, target_repeats = 4) {
  
  n              <- length(x)
  ideal_len      <- n / target_repeats
  lengths_sorted <- seq_len(n %/% 2)[order(abs(seq_len(n %/% 2) - ideal_len))]
  
  for (plen in lengths_sorted) {
    pattern  <- x[seq_len(plen)]
    repeated <- rep(pattern, length.out = n)
    if (all(x == repeated[seq_len(n)]))
      return(pattern)
  }
  
  return(x)   # no pattern found — return full sequence unchanged
}


#' Simple row-decimation downsampler
#'
#' Keeps every `by`-th row. Fast but does not anti-alias; suitable for
#' visualisation only.
#'
#' @param input_df  data.frame  Input data.
#' @param by        integer  Keep 1 in every `by` rows (default: 100).
#'
#' @return data.frame with ~1/by as many rows.
downSamp <- function(input_df, by = 100) {
  input_df[seq(1, nrow(input_df), by = by), , drop = FALSE]
}


# ── Section 5: Channel / Mux Read Pipelines ──────────────────────────────────

#' Mid-resolution channel/mux reader with QC plots
#'
#' Extracts raw current traces for a list of channel-mux combinations from a
#' FAST5 file, applies optional offset correction, generates per-cm plots,
#' and returns either full trace data or per-voltage insertion metrics.
#'
#' @param direc         character  Root directory.
#' @param folder        character  Sub-folder containing the FAST5 file.
#' @param file_name     character  FAST5 filename.
#' @param cm_list       character vector  Channel-mux strings, e.g. c("11","22").
#' @param out_file_name character/NULL  PDF output name stem (default: "cm_plots").
#' @param parallel      logical  TRUE uses mclapply across available cores.
#' @param offset_fix    logical  TRUE applies offsetCorrect().
#' @param all_data      logical  TRUE returns full trace; FALSE returns metrics.
#' @param plot_fast     logical  TRUE uses faster downsampled plot rendering.
#' @param ...           Forwarded to warmUpDetector, offsetCorrect, currPlot,
#'                       and insertionAnalyser.
#'
#' @return data.table of per-voltage insertion metrics (all_data = FALSE) or
#'   full trace data (all_data = TRUE).
cmReadr <- function(direc,
                    folder,
                    file_name,
                    cm_list,
                    out_file_name = NULL,
                    parallel      = FALSE,
                    offset_fix    = TRUE,
                    all_data      = FALSE,
                    plot_fast     = FALSE,
                    ...) {
  
  dots      <- list(...)
  safe_call <- function(fun, args)
    do.call(fun, c(args, dots[names(dots) %in% names(formals(fun))]))
  
  full_path <- file.path(direc, folder, file_name)
  run_id    <- substr(substr(file_name, nchar(file_name) - 10, nchar(file_name)), 1, 5)
  
  message("cmReadr: detecting warm-up and voltage steps for run ", run_id)
  init_rmv <- safe_call(
    warmUpDetector,
    list(full_path, return_mux_times = "y", return_voltage_steps = "y")
  )
  
  voltage_df        <- init_rmv[[3]]
  voltage_df        <- voltage_df[order(voltage_df$time), ]
  voltage_df        <- voltage_df[!is.na(voltage_df$time), ]
  
  channels <- unique(substr(cm_list, 1, nchar(cm_list) - 1))
  message("cmReadr: processing ", length(channels), " channel(s)")
  
  process_channel <- function(i) {
    tryCatch({
      ch <- channels[i]
      message("  Working on channel ", ch)
      cdata <- channelExtract(ch, full_path)
      cdata <- muxSplitr(cdata, ch, init_rmv)
      cdata <- subset(cdata, cdata$cm %in% cm_list)
      if (!nrow(cdata)) return(NULL)
      
      lapply(split(cdata, cdata$cm), function(cm_df) {
        if (nrow(cm_df) < 2) {
          message("  Skipping cm ", unique(cm_df$cm), " — < 2 data points")
          return(NULL)
        }
        cm_df$rel_time <- cm_df$time - cm_df$time[1]
        cm_df$voltage  <- voltage_df$voltage[findInterval(cm_df$time, voltage_df$time)]
        
        if (isTRUE(offset_fix))
          cm_df <- safe_call(offsetCorrect, list(cm_df))
        
        plt <- safe_call(currPlot, list(cm_df, rel_time = FALSE,
                                        plot_fast = plot_fast)) +
          ggplot2::facet_grid(cm ~ .) + My_Theme() +
          ggplot2::ylab("Current (pA)") + ggplot2::xlab("Time (s)")
        
        insdf <- if (!isTRUE(all_data)) {
          safe_call(insertionAnalyser, list(cm_df))
        } else {
          cm_df
        }
        list(plt = plt, insdf = insdf)
      })
    }, error = function(e) {
      message("  Error on channel ", channels[i], ": ", e$message)
      NULL
    })
  }
  
  result <- if (isTRUE(parallel)) {
    message("cmReadr: running in parallel mode")
    parallel::mclapply(seq_along(channels), process_channel,
                       mc.cores = parallel::detectCores())
  } else {
    message("cmReadr: running in single-core mode")
    lapply(seq_along(channels), process_channel)
  }
  
  flat    <- Filter(Negate(is.null), unlist(result, recursive = FALSE))
  plt_list <- Filter(
    function(p) inherits(p, "gg"),
    lapply(flat, `[[`, "plt")
  )
  combined <- rbindlist(
    Filter(Negate(is.null), lapply(flat, `[[`, "insdf")),
    use.names = TRUE, fill = TRUE
  )
  
  out_file <- file.path(direc, folder,
                        paste0(if (is.null(out_file_name)) "cm_plots" else out_file_name, "-", run_id, ".pdf")
  )
  
  if (length(plt_list) > 0) {
    arranged <- gridExtra::marrangeGrob(grobs = plt_list, nrow = 2, ncol = 2)
    message("cmReadr: saving ", length(plt_list), " plot(s) → ", out_file)
    tryCatch(
      ggplot2::ggsave(out_file, arranged, width = 29.7, height = 21, units = "cm"),
      error = function(e) message("cmReadr: ggsave failed — ", e$message)
    )
  } else {
    message("cmReadr: no valid plots — PDF not written.")
  }
  
  rhdf5::h5closeAll()
  return(combined)
}


#' High-throughput channel/mux reader for IV and event analysis
#'
#' The main per-channel processing pipeline. Reads raw signal, assigns mux
#' labels, applies offset correction, generates PNG plots, and returns either
#' per-voltage insertion metrics (for IV curves) or the full current trace.
#' Supports resuming from a previously saved insertiondf CSV.
#'
#' @param direc           character  Base directory.
#' @param folder          character  Sub-folder containing the FAST5 file.
#' @param file_name       character  FAST5 filename.
#' @param cm_list         character/list  Channel-mux strings to process.
#' @param all_data        logical  TRUE returns full trace; FALSE (default)
#'                          returns insertionAnalyser metrics.
#' @param re_analyse      logical  TRUE re-analyses even if a CSV already exists.
#' @param parallel        logical  TRUE uses mclapply on half available cores.
#' @param offset_fix      logical  TRUE applies offsetCorrect().
#' @param downsample      numeric/NULL  Decimation factor.
#' @param save_plots      logical  TRUE saves per-cm PNG plots.
#' @param init_rmv        list/NULL  Pre-computed warmUpDetector output. Computed
#'                          internally when NULL (pass it in from eventAna to
#'                          avoid redundant recomputation).
#' @param save_csv        logical  TRUE (default) writes insertiondf and
#'                          per-channel raw CSVs.
#' @param ...             Forwarded by name to: offsetCorrect (method),
#'                          besselFilt (filter_order, filter_cutoff,
#'                          filter_type, sample_rate), opcFinder (sd_adjust),
#'                          eventDetect.
#'
#' @return data.table of insertion metrics or full trace data.
runReadr <- function(direc,
                     folder,
                     file_name,
                     cm_list,
                     all_data    = FALSE,
                     re_analyse  = FALSE,
                     parallel    = FALSE,
                     offset_fix  = FALSE,
                     downsample  = NULL,
                     save_plots  = TRUE,
                     init_rmv    = NULL,
                     save_csv    = TRUE,
                     ...) {
  
  dots              <- list(...)
  offset_args       <- dots[names(dots) %in% "method"]
  
  full_path         <- file.path(direc, folder, file_name)
  dir_path          <- file.path(direc, folder)
  run_id            <- substr(substr(file_name, nchar(file_name) - 10,
                                     nchar(file_name)), 1, 5)
  plot_dir          <- file.path(dir_path, paste0(run_id, "_PLOTS"))
  
  # ── Parse cm_list into a data frame of channel + mux ────────────────────────
  input_channels      <- cmSplit(as.character(unlist(cm_list)))
  input_channels$cm   <- paste0(input_channels$channel, input_channels$mux)
  cm_list_chr         <- as.character(unlist(cm_list))
  
  # ── Resume: load previously saved insertion metrics ──────────────────────────
  csv_files <- list.files(dir_path)
  csv_files <- csv_files[grepl(run_id, csv_files) & grepl("insertiondf", csv_files)]
  
  if (isTRUE(all_data)) re_analyse <- TRUE
  
  if (!isTRUE(re_analyse) && length(csv_files) > 0) {
    message("runReadr: loading previously saved insertion metrics")
    if (length(csv_files) > 1)
      stop("runReadr: multiple insertiondf CSVs found for run ", run_id, " — ambiguous.")
    out_df <- try(data.table::fread(file.path(dir_path, csv_files[1]), fill = TRUE))
  }
  
  if (exists("out_df")) {
    done_cms       <- as.character(unique(out_df$cm))
    input_channels <- subset(input_channels, !input_channels$cm %in% done_cms)
    message("runReadr: skipping ", length(done_cms), " previously analysed cm(s)")
  }
  
  channels <- unique(stats::na.omit(input_channels[, 1]))
  
  # ── Warm-up detection ────────────────────────────────────────────────────────
  if (is.null(init_rmv)) {
    message("runReadr: running warmUpDetector for run ", run_id)
    init_rmv <- tryCatch(
      warmUpDetector(full_path, return_mux_times = "y"),
      error = function(e) {
        warning("runReadr: warmUpDetector failed — ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(init_rmv) || length(init_rmv[[2]]) == 0) {
      warning("runReadr: no mux change times detected — skipping ", file_name)
      return(NULL)
    }
  }
  
  # ── Per-channel processing ───────────────────────────────────────────────────
  if (length(channels) == 0) {
    message("runReadr: no new channels to process — returning cached data")
  } else {
    message("runReadr: processing channels: ", paste(channels, collapse = ", "))
    
    if (isTRUE(save_plots) && !dir.exists(plot_dir))
      dir.create(plot_dir, recursive = TRUE)
    
    apply_fun <- if (isTRUE(parallel)) {
      cores <- max(1L, parallel::detectCores(logical = FALSE) %/% 2L)
      message("runReadr: parallel mode — ", cores, " cores")
      function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
    } else {
      message("runReadr: single-core mode")
      lapply
    }
    
    outdat <- rbindlist(apply_fun(channels, function(working_chan) {
      message("  Reading channel ", working_chan)
      df <- rawImport(direc, folder, file_name, working_chan,
                      save_csv = save_csv, ...)
      
      if (data.table::is.data.table(df)) df <- as.data.frame(df)
      if ("V1" %in% colnames(df))
        df <- df[, !colnames(df) %in% "V1", drop = FALSE]
      
      keep <- intersect(c("current", "time", "voltage", "channel"), colnames(df))
      df   <- df[, keep, drop = FALSE]
      
      if (!is.null(downsample)) {
        if (!is.numeric(downsample)) {
          warning("runReadr: downsample must be numeric or NULL — ignoring")
        } else {
          df <- downSamp(df, downsample)
          if (data.table::is.data.table(df)) df <- as.data.frame(df)
        }
      }
      
      df <- muxSplitr(df, working_chan, init_rmv)
      df <- df[df$cm %in% cm_list_chr, , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
      
      cm_result <- rbindlist(lapply(split(df, df$cm), function(cm_df) {
        data.table::setDT(cm_df)
        data.table::set(cm_df, j = "rel_time", value = cm_df$time - cm_df$time[1])
        cm_df <- cm_df[rel_time > (max(rel_time) / 100)]
        data.table::set(cm_df, j = "rel_time", value = cm_df$rel_time - cm_df$rel_time[1])
        data.table::set(cm_df, j = "current",  value = as.numeric(cm_df$current))
        
        if (isTRUE(offset_fix))
          cm_df <- do.call(offsetCorrect, c(list(cm_df), offset_args))
        
        if (isTRUE(save_plots)) {
          ds_df <- downSamp(cm_df)
          plt   <- currPlot(ds_df, rel_time = FALSE, plot_fast = FALSE) +
            My_Theme() +
            ggplot2::ylab("Current (pA)") +
            ggplot2::xlab("Time (s)")
          png_path <- file.path(plot_dir,
                                paste0("channel-mux_", unique(ds_df$cm), ".png"))
          grDevices::png(png_path, width = 10, height = 8, units = "in", res = 150)
          print(plt)
          grDevices::dev.off()
          rm(ds_df, plt)
        }
        
        if (!isTRUE(all_data)) {
          metrics <- insertionAnalyser(cm_df)
          rm(cm_df)
          metrics
        } else {
          cm_df
        }
      }), use.names = TRUE, fill = TRUE)
      
      rm(df)
      gc(verbose = FALSE)
      cm_result
      
    }), use.names = TRUE, fill = TRUE)
  }
  
  # ── Save insertion metrics ────────────────────────────────────────────────────
  if (!isTRUE(all_data) && exists("outdat") && isTRUE(save_csv)) {
    csv_out <- file.path(dir_path, paste0(run_id, "_insertiondf.csv"))
    message("runReadr: writing insertion metrics → ", csv_out)
    data.table::fwrite(data.table::data.table(outdat), csv_out, append = TRUE)
  }
  
  if (exists("out_df") && exists("outdat")) {
    out_df <- rbindlist(list(out_df, outdat), use.names = TRUE, fill = TRUE)
  } else if (exists("outdat")) {
    out_df <- outdat
  }
  
  data.table::set(out_df, j = "run_id", value = run_id)
  out_df <- out_df[cm %in% cm_list_chr]
  message("runReadr: done — ", nrow(out_df), " rows returned")
  return(out_df)
}


# ── Section 6: Event Detection ────────────────────────────────────────────────

#' Full event-detection pipeline for a list of FAST5 files
#'
#' For each file, detects the warm-up period, loops over channels and mux
#' combinations, optionally filters the trace, and calls eventDetect to locate
#' translocation events. Results are written incrementally to CSV files.
#'
#' @param direc               character  Working directory.
#' @param folder              character  Sub-folder within `direc` (default: "").
#' @param data_file_list      list  FAST5 filenames to process.
#' @param cm_list_list        list  List of cm vectors, one per file.
#' @param csv_name_list       character  Output CSV name stems, one per file.
#' @param use_filter          logical  Apply Butterworth filter? (default: TRUE).
#' @param sample_rate         numeric/NULL  Sampling rate in Hz. NULL = auto-detect.
#' @param parallel            logical  Use parallel processing? (default: FALSE).
#' @param offset_fix          logical  Apply offset correction? (default: TRUE).
#' @param resume              logical  Skip previously processed cms? (default: TRUE).
#' @param ...                 Forwarded to besselFilt (filter_cutoff, filter_order,
#'                              filter_type), offsetCorrect (method),
#'                              opcFinder (sd_adjust), eventDetect
#'                              (event_dwell_threshold).
#'
#' @return Invisible NULL. Results are written to CSV files.
eventAna <- function(direc,
                     folder            = "",
                     data_file_list,
                     cm_list_list,
                     csv_name_list,
                     use_filter        = TRUE,
                     sample_rate       = 4000,
                     parallel          = FALSE,
                     offset_fix        = TRUE,
                     resume            = TRUE,
                     ...) {
  
  dots <- list(...)
  
  filter_args <- list(
    filter_order  = dots$filter_order  %||% 2,
    filter_cutoff = dots$filter_cutoff %||% 1500,
    filter_type   = dots$filter_type   %||% "low"
  )
  offset_args <- dots[names(dots) %in% "method"]
  opc_args    <- dots[names(dots) %in% "sd_adjust"]
  event_args  <- dots[names(dots) %in% "event_dwell_threshold"]
  
  base_dir <- if (nchar(folder) > 0) file.path(direc, folder) else direc
  
  for (i in seq_along(data_file_list)) {
    
    data_file <- data_file_list[[i]]
    cm_list   <- as.character(cm_list_list[[i]])
    csv_out   <- file.path(base_dir, paste0(csv_name_list[[i]], ".csv"))
    full_path <- file.path(base_dir, data_file)
    run_id    <- substr(data_file, max(1, nchar(data_file) - 35),
                        nchar(data_file) - 7)
    
    message("\n=== File ", i, "/", length(data_file_list), ": ",
            basename(data_file), " (run: ", run_id, ") ===")
    
    # ── Resume: skip already-processed cms ──────────────────────────────────
    if (isTRUE(resume) && file.exists(csv_out)) {
      done_cms <- unique(data.table::fread(csv_out, select = "cm")$cm)
      cm_list  <- setdiff(cm_list, done_cms)
      if (length(done_cms) > 0)
        message("Skipping ", length(done_cms), " previously processed cm(s)")
    }
    if (length(cm_list) == 0) { message("All cms processed — skipping file"); next }
    
    # ── Warm-up detection ────────────────────────────────────────────────────
    message("Detecting voltage programme...")
    init_rmv <- tryCatch(
      warmUpDetector(full_path),
      error = function(e) { warning("warmUpDetector failed: ", e$message); NULL }
    )
    if (is.null(init_rmv)) next
    
    # ── Resolve sample rate ──────────────────────────────────────────────────
    if (is.null(sample_rate)) {
      sample_rate <- tryCatch({
        meta <- rhdf5::h5readAttributes(full_path, "/IntermediateData/Channel_1/Meta/")
        meta$sample_rate
      }, error = function(e) { message("Using default 4000 Hz"); 4000 })
    } else if (sample_rate < 100) {
      sample_rate <- sample_rate * 1000   # kHz → Hz
    }
    message("Sample rate: ", sample_rate, " Hz")
    
    cm_split <- cmSplit(cm_list)
    channels <- as.integer(unique(cm_split$channel))
    message("Processing ", length(channels), " channel(s)")
    
    for (ch_idx in seq_along(channels)) {
      ch <- channels[ch_idx]
      message(sprintf("  Channel %d (%d/%d)", ch, ch_idx, length(channels)))
      
      raw_data <- tryCatch(
        rawImport(direc = base_dir, folder = "", file_name = data_file,
                  channel = ch, save_csv = FALSE, re_extract = TRUE),
        error = function(e) {
          warning("  Failed to extract channel ", ch, ": ", e$message)
          NULL
        }
      )
      if (is.null(raw_data)) next
      
      raw_data          <- muxSplitr(raw_data, ch, init_rmv)
      raw_data          <- raw_data[raw_data$cm %in% cm_list, ]
      if (nrow(raw_data) == 0) next
      raw_data$run_id   <- run_id
      data.table::setDT(raw_data)
      
      if (isTRUE(offset_fix))
        raw_data <- do.call(offsetCorrect, c(list(raw_data), offset_args))
      
      raw_data <- do.call(opcFinder, c(list(raw_data), opc_args))
      if (is.null(raw_data) || nrow(raw_data) == 0) next
      
      for (cm_id in unique(raw_data$cm)) {
        cm_data <- raw_data[cm == cm_id]
        
        if (isTRUE(use_filter)) {
          cm_data <- do.call(
            besselFilt,
            c(list(input_df = cm_data, sample_rate = sample_rate,
                   data_path = full_path), filter_args)
          )
        }
        
        do.call(eventDetect, c(list(input_df = cm_data, csv_out_path = csv_out),
                               event_args))
        rm(cm_data)
      }
      rm(raw_data); gc(verbose = FALSE)
    }
  }
  
  message("\n=== eventAna: complete ===")
  invisible(NULL)
}


#' Detect translocation events in a filtered current trace
#'
#' Identifies current blockade events as contiguous segments where the
#' absolute current falls below the open-pore current estimate. Events
#' shorter than `min_dwell` are discarded. Results are appended to a CSV.
#'
#' @param input_df    data.table  Requires: current, time, voltage, opc, cm.
#'                      Optionally: run_id.
#' @param csv_out_path character  Path to the output CSV (appended to).
#' @param min_dwell   numeric  Minimum event duration in seconds (default: 0.001).
#' @param ...         Reserved; ignored.
#'
#' @return Invisible NULL. Events are written to `csv_out_path`.
eventDetect <- function(input_df,
                        csv_out_path,
                        min_dwell = 0.001,
                        ...) {
  
  if (!data.table::is.data.table(input_df))
    data.table::setDT(input_df)
  
  if (!"run_id" %in% names(input_df)) {
    input_df[, run_id := "unknown_run"]
    warning("eventDetect: 'run_id' column missing — using placeholder")
  }
  
  # Baseline at 0 mV
  zero_curr <- if (any(input_df$voltage == 0)) {
    tryCatch(
      openPoreCalc(input_df[voltage == 0, current]),
      error = function(e) { message("  Could not estimate 0 mV baseline; using 0"); 0 }
    )
  } else {
    message("  No 0 mV data; using 0 as baseline")
    0
  }
  
  voltages <- unique(input_df$voltage)
  message(sprintf("  eventDetect: scanning %d voltage(s)", length(voltages)))
  
  for (v in voltages) {
    dt <- input_df[voltage == v]
    if (nrow(dt) < 3) next
    
    dt <- dt[abs(current) < abs(opc)]
    if (nrow(dt) < 2) next
    
    dt[, time := time - time[1]]
    dt[, time_diff := c(0, diff(time))]
    
    expected_step <- stats::median(dt$time_diff[dt$time_diff > 0])
    if (is.na(expected_step) || expected_step == 0)
      expected_step <- mean(diff(head(dt$time, 100)))
    
    dt[, event_id := cumsum(time_diff > (expected_step * 1.5))]
    
    event_sizes   <- dt[, .N, by = event_id]
    valid_ids     <- event_sizes[N >= 2, event_id]
    dt            <- dt[event_id %in% valid_ids]
    if (nrow(dt) == 0) next
    
    events_list   <- vector("list", length(unique(dt$event_id)))
    event_counter <- 0L
    
    for (eid in unique(dt$event_id)) {
      ev       <- dt[event_id == eid]
      if (nrow(ev) < 2) next
      duration <- max(ev$time) - min(ev$time)
      if (duration < min_dwell) next
      
      avg_curr <- stats::median(ev$current, na.rm = TRUE)
      sd_curr  <- stats::sd(ev$current, na.rm = TRUE)
      if (is.na(sd_curr)) sd_curr <- 0
      
      amplitude <- tryCatch(
        100 * norm(c(zero_curr, avg_curr, ev$opc[1]), type = "2"),
        error = function(e) 100 * abs(avg_curr / ev$opc[1])
      )
      
      event_counter <- event_counter + 1L
      events_list[[event_counter]] <- data.table::data.table(
        cm          = ev$cm[1],
        avg_current = avg_curr,
        sd_current  = sd_curr,
        duration    = duration,
        amplitude   = amplitude,
        opc         = ev$opc[1],
        run_id      = ev$run_id[1],
        voltage     = v
      )
    }
    
    if (event_counter > 0) {
      events_df <- rbindlist(events_list[seq_len(event_counter)])
      data.table::fwrite(
        events_df, csv_out_path,
        append    = file.exists(csv_out_path),
        col.names = !file.exists(csv_out_path)
      )
      message(sprintf("    %+4d mV: %d event(s)", v, event_counter))
    }
  }
  
  invisible(NULL)
}


# ── Section 7: Metric Extraction ─────────────────────────────────────────────

#' Summarise a current trace into per-voltage statistics
#'
#' Returns one row per voltage step with current averages, noise (SD), linear
#' drift (slope), curvature, and conductance.
#'
#' @param input_df       data.frame  Requires: time, current, voltage, cm.
#' @param avg_type       character  Central tendency: "mean", "median", or
#'                         "mode" (default). Mode uses modalAvg().
#' @param voltage_subset numeric/NULL  Restrict to these voltages only.
#'
#' @return data.frame with one row per voltage step and columns:
#'   per_voltage_current, per_voltage_current_sd, voltage,
#'   per_voltage_slope, per_voltage_curve, cm,
#'   per_voltage_cond, per_chan_curve, per_chan_cond, per_chan_cond_stdev.
insertionAnalyser <- function(input_df,
                              avg_type       = "mode",
                              voltage_subset = NULL) {
  
  input_df <- stats::na.omit(input_df)
  
  if (!is.null(voltage_subset))
    input_df <- subset(input_df, input_df$voltage %in% voltage_subset)
  
  vsteps <- unique(input_df$voltage)
  
  outdf <- rbindlist(lapply(vsteps, function(v) {
    step_df <- subset(input_df, input_df$voltage == v)
    # Drop the first 5 % of each voltage step to skip transition artefacts
    step_df <- step_df[ceiling(nrow(step_df) / 20):nrow(step_df), ]
    
    mavg <- switch(
      avg_type,
      mean   = mean(step_df$current),
      median = stats::median(step_df$current),
      mode   = modalAvg(step_df$current),
      {
        warning("insertionAnalyser: avg_type must be 'mean', 'median', or 'mode' — using mode.")
        modalAvg(step_df$current)
      }
    )
    
    vsd <- stats::sd(step_df$current)
    
    # Fast linear and quadratic fits using .lm.fit() (avoids formula overhead)
    t   <- step_df$time
    I   <- step_df$current
    t_c <- t - mean(t)
    t2  <- t_c^2 - mean(t_c^2)
    
    slope <- .lm.fit(cbind(1, t_c), I)$coefficients[2]
    curve <- .lm.fit(cbind(1, t_c, t2), I)$coefficients[3]
    
    data.frame(
      per_voltage_current    = mavg,
      per_voltage_current_sd = vsd,
      voltage                = v,
      per_voltage_slope      = slope,
      per_voltage_curve      = curve
    )
  }), use.names = TRUE)
  
  outdf$cm <- unique(input_df$cm)
  
  outdf$per_voltage_cond <- outdf$per_voltage_current / outdf$voltage
  outdf$per_voltage_cond[outdf$voltage == 0] <- 0
  
  outdf$per_chan_curve      <- mean(
    abs(outdf$per_voltage_curve[outdf$voltage != 0]) /
      abs(outdf$voltage[outdf$voltage != 0])
  )
  outdf$per_chan_cond       <- mean(outdf$per_voltage_cond)
  outdf$per_chan_cond_stdev <- stats::sd(outdf$per_voltage_cond)
  
  return(outdf)
}


#' Membrane QC from IV data using per-channel current statistics
#'
#' Classifies each channel-mux as "Good", "Intact", "Disrupted", or "Burst"
#' based on current SD and maximum current thresholds. Optionally saves a
#' summary CSV and returns data in a format compatible with barplotr.
#'
#' @param direc               character  Working directory.
#' @param folder              character  Sub-folder.
#' @param file_name           character  FAST5 filename.
#' @param re_analyse          logical  TRUE re-analyses even if a CSV exists.
#' @param channels            integer vector  Channels to analyse (default: 1–512).
#' @param sd_thresholds       numeric vector  c(good_sd, intact_sd). Channels
#'                              with SD < good_sd are "Good"; good_sd–intact_sd
#'                              are "Intact"; above intact_sd are "Disrupted".
#' @param burst_threshold     numeric  Channels whose max |current| exceeds this
#'                              are classified "Burst" (default: 770 pA).
#' @param output_for_plot     logical  TRUE returns a frequency table compatible
#'                              with barplotr (default: TRUE).
#' @param remove_warmup       logical  FALSE disables warm-up removal (for old
#'                              Minnow data).
#' @param id                  character/NULL  Optional label appended to output.
#' @param ...                 Forwarded to warmUpDetector.
#'
#' @return data.frame of membrane quality classifications, either raw per-cm or
#'   aggregated counts (when output_for_plot = TRUE).
memQCfromIV <- function(direc,
                        folder,
                        file_name,
                        re_analyse       = FALSE,
                        channels         = seq_len(512),
                        sd_thresholds    = c(25, 70),
                        burst_threshold  = 770,
                        output_for_plot  = TRUE,
                        remove_warmup    = TRUE,
                        id               = NULL,
                        ...) {
  
  dir_path  <- file.path(direc, folder)
  full_path <- file.path(dir_path, file_name)
  run_id    <- substr(substr(file_name, nchar(file_name) - 10, nchar(file_name)), 1, 5)
  
  csv_files <- list.files(dir_path)
  csv_files <- csv_files[grepl(run_id, csv_files) & grepl("_memQCfromIV", csv_files)]
  
  if (!isTRUE(re_analyse) && length(csv_files) > 0) {
    message("memQCfromIV: loading previously saved QC data")
    if (length(csv_files) > 1)
      stop("memQCfromIV: multiple memQC files found for run ", run_id)
    out_df <- try(data.table::fread(file.path(dir_path, csv_files[1]), fill = TRUE))
    
  } else {
    message("memQCfromIV: analysing ", length(channels), " channels")
    
    if (isTRUE(remove_warmup))
      init_rmv <- warmUpDetector(full_path, channel = 1, return_mux_times = "y", ...)
    
    out_df <- do.call(rbind, lapply(channels, function(ch) {
      cdata <- channelExtract(ch, full_path)
      
      cdata <- if (isTRUE(remove_warmup)) {
        muxSplitr(cdata, channel = ch, init_rmv = init_rmv)
      } else {
        muxSplitOld(cdata)
      }
      
      do.call(rbind, lapply(split(cdata, cdata$cm), function(cm_df) {
        data.frame(
          curr_avg = mean(cm_df$current),
          curr_sd  = stats::sd(cm_df$current),
          curr_max = max(abs(cm_df$current)),
          cm       = unique(cm_df$cm)
        )
      }))
    }))
    
    out_df$quality <- ifelse(
      out_df$curr_max > burst_threshold, "Burst",
      ifelse(out_df$curr_sd < sd_thresholds[1], "Good",
             ifelse(out_df$curr_sd < sd_thresholds[2], "Intact", "Disrupted"))
    )
    
    message("memQCfromIV: saving → ",
            file.path(dir_path, paste0(run_id, "_memQCfromIV.csv")))
    data.table::fwrite(out_df,
                       file.path(dir_path, paste0(run_id, "_memQCfromIV.csv")))
  }
  
  if (!is.null(id)) out_df$exp <- id
  
  if (isTRUE(output_for_plot)) {
    freq_df <- as.data.frame(table(out_df$quality))
    return(data.frame(y = freq_df$Var1, id = freq_df$Freq))
  }
  
  return(out_df)
}


#' Read ONT membrane QC CSVs from a folder tree
#'
#' Recursively finds all *_memQC*.csv files (excluding previously generated
#' _memQCfromIV files), aggregates membrane quality counts, and returns a
#' long-format data frame suitable for barplotr.
#'
#' @param folder  character  Root folder to search recursively.
#' @param id      character/NULL  Optional label to overwrite the `id` column.
#'
#' @return data.frame with columns: x (quality label), y (count), id (run path),
#'   run (run path copy for compatibility).
readONTMemQC <- function(folder, id = NULL) {
  
  csv_files <- list.files(folder, pattern = "\\.csv$",
                          recursive = TRUE, full.names = TRUE)
  csv_files <- csv_files[!grepl("_memQCfromIV", csv_files)]
  message("readONTMemQC: found ", length(csv_files), " QC CSV(s)")
  
  out_df <- do.call(rbind, lapply(csv_files, function(f) {
    mqc   <- read.csv(f)
    cnts  <- as.data.frame(table(mqc$membrane_quality))
    names(cnts) <- c("x", "y")
    cnts$id <- gsub("\\.csv$", "", f)
    cnts
  }))
  
  if (!is.null(id)) out_df$id <- id
  out_df$run <- out_df$id
  out_df$id  <- out_df$x
  return(out_df)
}


# ── Section 8: Conductance / Size Conversion ─────────────────────────────────

#' Convert between nanopore conductance and estimated pore diameter
#'
#' Implements the Howorka 2020 Nat. Protoc. formula:
#'   G = k × (π d² / 4L  +  π d)
#' Pass G to solve for d, or d to solve for G. Pass exactly one of KCl or NaCl.
#'
#' @param G         numeric/NULL  Conductance in nS.
#' @param KCl       numeric/NULL  KCl concentration in M.
#' @param NaCl      numeric/NULL  NaCl concentration in M.
#' @param d         numeric/NULL  Pore diameter in nm.
#' @param pore_length  numeric  Pore length in nm (default: 3.74 nm, lipid bilayer).
#' @param G_error   numeric/NULL  Conductance measurement error for error
#'                   propagation. NULL skips error estimation.
#'
#' @return Invisibly: list(estimate, error). Also prints a summary line.
cond2size <- function(G          = NULL,
                      KCl        = NULL,
                      NaCl       = NULL,
                      d          = NULL,
                      pore_length = 3.74,
                      G_error    = NULL) {
  
  # Molar conductivities (S m² mol⁻¹); ×100 converts to S cm⁻¹
  k_m  <- 0.0735; na_m <- 0.05011; cl_m <- 0.07635
  
  if (!is.null(KCl) && !is.null(NaCl))
    stop("cond2size: pass only one of KCl or NaCl.")
  
  k <- if (!is.null(KCl))  (k_m  + cl_m) * KCl  * 100
  else if (!is.null(NaCl)) (na_m + cl_m) * NaCl * 100
  else stop("cond2size: supply a concentration via KCl or NaCl.")
  
  if (!is.null(G) && is.null(d)) {
    d <- (sqrt(pi) * sqrt(pi * G^2 + 16 * G * k * pore_length) + pi * G) /
      (2 * pi * k)
    d <- round(d, 4)
    
    if (!is.null(G_error)) {
      dG_dd   <- k * (2 * pi * d * (4 * pore_length - pi * d)) /
        (4 * pore_length + pi * d)^2
      d_error <- G_error / abs(dG_dd)
      cat("Pore diameter:", d, "nm ±", round(d_error, 4), "nm",
          "(conductance:", G, "nS)\n")
      return(invisible(list(d, d_error)))
    }
    cat("Pore diameter:", d, "nm (conductance:", G, "nS)\n")
    return(invisible(list(d, NULL)))
    
  } else if (!is.null(d) && is.null(G)) {
    G <- k * (pi * d^2) / (4 * pore_length + pi * d)
    G <- round(G, 4)
    cat("Approximate conductance:", G, "nS (pore:", d, "nm ×", pore_length, "nm)\n")
    return(invisible(list(G, NULL)))
    
  } else {
    stop("cond2size: pass exactly one of G or d.")
  }
}


# ── Section 9: Spectral Analysis ─────────────────────────────────────────────

#' Compute a zero-padded power spectral density
#'
#' Internal helper used by powerSpecCalc.
#'
#' @param current_values  numeric  Current samples (any unit).
#' @param sample_rate     numeric  Sampling rate in Hz (default: 4000).
#'
#' @return data.frame with columns: power, freq.
powerSpecDens <- function(current_values, sample_rate = 4000) {
  
  N   <- length(current_values)
  M   <- 2^(floor(log2(N)) + 2)
  xzp <- c(current_values, rep(0, M - N))
  
  power <- (1 / N) * abs(stats::fft(xzp)[seq_len(M / 2 + 1)])^2
  freq  <- seq(sample_rate / 1000, sample_rate, length.out = length(power))
  
  data.frame(power = power, freq = freq)
}


#' Compute the PSD for a single channel-mux at a given voltage
#'
#' @param input_df          data.frame  Requires: cm, voltage, current.
#' @param cm_id             character/NULL  Cm identifier to subset to.
#' @param voltage           numeric  Voltage step to analyse.
#' @param sample_rate       numeric  Sampling rate in Hz (default: 4000).
#'
#' @return data.frame with columns: power, freq, id, voltage.
powerSpecCalc <- function(input_df,
                          cm_id       = NULL,
                          voltage,
                          sample_rate = 4000) {
  
  sub_df <- if (!is.null(cm_id)) {
    subset(input_df, input_df$cm == cm_id & input_df$voltage == voltage)
  } else {
    input_df
  }
  
  if ("current" %in% colnames(sub_df)) {
    vals <- sub_df$current / 1000   # pA → nA
  } else {
    message("powerSpecCalc: 'current' column not found — using per_channel_current")
    vals <- sub_df$per_channel_current / 1000
  }
  
  psd         <- powerSpecDens(vals, sample_rate)
  psd$id      <- cm_id
  psd$voltage <- voltage
  psd
}


#' Plot power spectral density on log-log axes
#'
#' @param psd_df    data.frame  Output of powerSpecCalc. Requires: freq,
#'                   power, id, voltage.
#' @param annotate_voltage logical  TRUE annotates the plot with the analysed
#'                              voltage (default: TRUE).
#' @param cols      character/NULL  Colours per `id`; uses ggplot defaults if NULL.
#' @param txt_size  numeric  Base font size (default: 7).
#'
#' @return ggplot object.
powerSpecPlot <- function(psd_df,
                          annotate_voltage = TRUE,
                          cols             = NULL,
                          txt_size         = 7) {
  
  #psd_df$id <- as.character(psd_df$id)
  
  plt <- ggplot2::ggplot(psd_df, ggplot2::aes(x = freq, y = power, colour = id)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::scale_x_continuous(
      trans  = "log10",
      labels = ~ format(.x, scientific = FALSE)
    ) +
    ggplot2::annotation_logticks(base = 10, sides = "bl") +
    ggplot2::xlab("\nFrequency / Hz") +
    ggplot2::ylab(expression("Power / nA"^"2" ~ "Hz"^"-1")) +
    ggplot2::theme_bw() +
    My_Theme(txt_size = txt_size)
  
  if (isTRUE(annotate_voltage))
    plt <- plt + ggplot2::annotate(
      "text", y = 100, x = 50,
      label = paste0("Calculated at ", unique(psd_df$voltage), " mV"),
      size  = (txt_size / ggplot2::.pt) / 1.5
    )
  
  if (!is.null(cols))
    plt <- plt + ggplot2::scale_colour_manual(values = cols)
  
  return(plt)
}


# ── Section 10: Plotting ──────────────────────────────────────────────────────

#' Plot current vs. time trace(s)
#'
#' The primary plotting function for representative current traces. Supports
#' absolute and relative time axes, optional secondary voltage axis, per-cm
#' or per-id colouring, and fast (scatter) or high-resolution (line) rendering.
#'
#' @param input_df      data.frame  Requires: time (s), current (pA). Optional:
#'                        voltage (mV), cm, id, stdev.
#' @param rel_time      logical  TRUE (default) resets time to zero per cm.
#' @param line_size     numeric  Geom line/point size (default: 0.25).
#' @param plot_ylim     numeric/NULL  Y-axis limits, e.g. c(-200, 200).
#' @param plot_xlim     numeric/NULL  X-axis limits.
#' @param show_stdev    logical  TRUE adds a ribbon for the `stdev` column.
#' @param colour        character  Single colour for all data (default: "black").
#' @param sec_axis      logical  TRUE adds a secondary y-axis for voltage.
#' @param colour_by     character of column name to colour plots by eg "id" or "cm"
#'                        any other value = single colour.
#' @param txt_size      numeric  Base font size (default: 7).
#' @param plot_fast     logical  TRUE uses geom_scattermore (faster, default).
#' @param facet         character/NULL  "l" facets by cm using facet_wrap.
#'
#' @return ggplot object.
currPlot <- function(input_df,
                     rel_time   = TRUE,
                     line_size  = 0.25,
                     plot_ylim  = NULL,
                     plot_xlim  = NULL,
                     show_stdev = TRUE,
                     colour     = "black",
                     sec_axis   = FALSE,
                     colour_by  = "none",
                     txt_size   = 7,
                     plot_fast  = TRUE,
                     facet      = NULL) {
  
  # ── Optionally reset time origin per channel-mux ────────────────────────────
  if (isTRUE(rel_time)) {
    if ("cm" %in% names(input_df)) {
      input_df <- do.call(rbind, lapply(split(input_df, input_df$cm), function(i) {
        i$time <- i$time - i$time[1]
        i
      }))
    } else {
      input_df$time <- input_df$time - input_df$time[1]
    }
  }
  
  # ── Validate colour column ─────────────────────────────────────────────
  use_colour <- !is.null(colour_by) &&
    tolower(colour_by) != "none" &&
    colour_by %in% names(input_df)
  
  # ── Build base geom ──────────────────────────────────────────────────────────
  geom_fun <- if (isTRUE(plot_fast)) geom_scattermore else ggplot2::geom_line
  
  if (use_colour) {
    
    out <- ggplot2::ggplot() +
      geom_fun(
        data = input_df,
        ggplot2::aes(
          x = time,
          y = current,
          colour = .data[[colour_by]]
        ),
        size = line_size
      )
    
    if (length(colour) > 1) {
      out <- out + ggplot2::scale_colour_manual(values = colour)
    }
    
  } else {
    
    out <- ggplot2::ggplot() +
      geom_fun(
        data = input_df,
        ggplot2::aes(x = time, y = current),
        size = line_size,
        colour = colour
      )
  }
  
  # ── Optional SD ribbon ───────────────────────────────────────────────────────
  if ("stdev" %in% colnames(input_df) && isTRUE(show_stdev)) {
    out <- out +
      ggplot2::geom_ribbon(
        data  = input_df,
        ggplot2::aes(x = time, ymin = current - stdev, ymax = current + stdev),
        alpha  = 0.25, colour = NA, fill = colour
      ) +
      ggplot2::geom_line(data = input_df,
                         ggplot2::aes(x = time, y = current),
                         size = line_size, colour = colour)
  }
  
  # ── Y axis + optional secondary voltage axis ─────────────────────────────────
  has_voltage <- "voltage" %in% colnames(input_df) && isTRUE(sec_axis)
  
  if (!is.null(plot_ylim)) {
    brk <- stefBreaks(plot_ylim)
    out <- out + ggplot2::scale_y_continuous(
      name   = "Current (pA)",
      breaks = brk[[1]], labels = brk[[2]], limits = plot_ylim,
      sec.axis = if (has_voltage)
        ggplot2::sec_axis(~ ., name = "Voltage (mV)",
                          breaks = brk[[1]], labels = brk[[2]])
      else ggplot2::waiver()
    )
    if (has_voltage)
      out <- out + ggplot2::geom_line(
        data   = input_df,
        ggplot2::aes(x = time, y = voltage),
        colour = "darkgrey"
      )
  } else {
    if (has_voltage) {
      out <- out +
        ggplot2::geom_line(data = input_df,
                           ggplot2::aes(x = time, y = voltage), colour = "darkgrey") +
        ggplot2::scale_y_continuous(
          name     = "Current (pA)",
          sec.axis = ggplot2::sec_axis(~ ., name = "Voltage (mV)")
        )
    } else {
      out <- out + ggplot2::ylab("Current (pA)")
    }
  }
  
  # ── X axis ───────────────────────────────────────────────────────────────────
  if (!is.null(plot_xlim)) {
    brk <- stefBreaks(plot_xlim)
    out <- out + ggplot2::scale_x_continuous(
      name = "Time (s)", breaks = brk[[1]], labels = brk[[2]], limits = plot_xlim
    )
  } else {
    out <- out + ggplot2::xlab("Time (s)")
  }
  
  # ── Optional faceting ────────────────────────────────────────────────────────
  if (!is.null(facet) && tolower(substr(facet, 1, 1)) == "l")
    out <- out + ggplot2::facet_wrap(~cm)
  
  out + My_Theme(txt_size = txt_size, base_size = txt_size)
}


#' Plot nanopore translocation events
#'
#' Supports an amplitude histogram (mode = "histogram") and a dwell-time vs.
#' amplitude scatter plot (mode = "dwells").
#'
#' @param input_df    data.frame  Output of eventDetect. Requires: amplitude,
#'                     duration, id.
#' @param mode        character  "histogram"/"h" or "dwells"/"d".
#' @param cols        character  Colour(s) per unique group (recycled if short).
#' @param n_bins      integer  Number of histogram bins (default: 15).
#' @param plot_ylim   numeric/NULL  Y-axis limits.
#' @param plot_xlim   numeric/NULL  X-axis limits.
#' @param log_axes    character  Axes to log-scale: "x", "y", "xy", or "none".
#' @param marker_size numeric  Point/tick size.
#' @param alpha_density logical  TRUE maps point transparency to count density
#'                        in dwell mode.
#' @param colour_by   character/NULL  Column name to use as the colour
#'                     aesthetic instead of "id" (e.g. "voltage").
#' @param txt_size    numeric  Base font size (default: 7).
#'
#' @return ggplot object.
eventPlotr <- function(input_df,
                       mode          = "histogram",
                       cols          = "darkgrey",
                       n_bins        = 15,
                       plot_ylim     = NULL,
                       plot_xlim     = NULL,
                       log_axes      = "none",
                       marker_size   = 1,
                       alpha_density = FALSE,
                       colour_by     = NULL,
                       txt_size      = 7) {
  
  colour_col <- if (!is.null(colour_by)) {
    if (!colour_by %in% colnames(input_df))
      stop("eventPlotr: colour_by column '", colour_by, "' not found.")
    colour_by
  } else "id"
  
  n_groups <- length(unique(input_df[[colour_col]]))
  if (n_groups > length(cols)) cols <- rep(cols, length.out = n_groups)
  
  log_x <- grepl("x", log_axes, ignore.case = TRUE)
  log_y <- grepl("y", log_axes, ignore.case = TRUE)
  
  # ── Axis-scale helpers (DRY) ─────────────────────────────────────────────────
  make_x_scale <- function(limits, log, reverse = FALSE) {
    brk   <- stefBreaks(limits, log = if (log) "y" else NULL)
    trans <- if (log && reverse) c("log10", "reverse")
    else if (log)        "log10"
    else if (reverse)    "reverse"
    else                 "identity"
    list(
      ggplot2::scale_x_continuous(trans = trans, breaks = brk[[1]],
                                  labels = brk[[2]], limits = limits),
      if (log) ggplot2::annotation_logticks(sides = "b", size = marker_size / 4)
    )
  }
  
  make_y_scale <- function(limits, log, reverse = FALSE) {
    brk   <- stefBreaks(limits, log = if (log) "y" else NULL)
    trans <- if (log && reverse) c("log10", "reverse")
    else if (log)        "log10"
    else if (reverse)    "reverse"
    else                 "identity"
    list(
      ggplot2::scale_y_continuous(trans = trans, breaks = brk[[1]],
                                  labels = brk[[2]], limits = limits),
      if (log) ggplot2::annotation_logticks(sides = "l", size = marker_size / 4)
    )
  }
  
  # ── Plot ─────────────────────────────────────────────────────────────────────
  mode_char <- substr(tolower(mode), 1, 1)
  
  if (mode_char == "d") {
    base_aes <- if (isTRUE(alpha_density)) {
      ggplot2::aes(x = duration, y = amplitude,
                   colour = .data[[colour_col]], alpha = ggplot2::after_stat(n))
    } else {
      ggplot2::aes(x = duration, y = amplitude, colour = .data[[colour_col]])
    }
    
    plt <- ggplot2::ggplot(input_df, base_aes) +
      #ggplot2::geom_count(size = marker_size) +
      #ggplot2::scale_size(range = c(marker_size, marker_size * 5))
      ggplot2::geom_point(size = marker_size) +## unsure why count was used here, reverting
      ggplot2::scale_color_manual(values = cols, name = "") +
      ggplot2::xlab("Dwell Time (s)") + ggplot2::ylab("Amplitude (%)") +
      { if (isTRUE(alpha_density)) ggplot2::guides(alpha = "none") }
    
    if (!is.null(plot_xlim)) plt <- plt + make_x_scale(plot_xlim, log_x)
    if (!is.null(plot_ylim)) plt <- plt + make_y_scale(plot_ylim, log_y, reverse = TRUE)
    
  } else if (mode_char == "h") {
    plt <- ggplot2::ggplot(input_df, ggplot2::aes(x = amplitude,
                                                  fill = .data[[colour_col]])) +
      ggplot2::geom_histogram(bins = n_bins) +
      ggplot2::scale_fill_manual(values = cols, name = "") +
      ggplot2::ylab("Frequency") + ggplot2::xlab("Amplitude (%)") +
      ggplot2::scale_x_continuous(trans = "reverse") +
      ggplot2::coord_flip()
    
    if (!is.null(plot_xlim)) plt <- plt + make_x_scale(plot_xlim, log_x, reverse = TRUE)
    if (!is.null(plot_ylim)) plt <- plt + make_y_scale(plot_ylim, log_y)
    
  } else {
    stop("eventPlotr: unknown mode '", mode, "'. Use 'histogram'/'h' or 'dwells'/'d'.")
  }
  
  plt +
    ggplot2::theme_bw(base_size = txt_size) +
    ggplot2::facet_grid(. ~ id) +
    My_Theme(txt_size = txt_size)
}

#' Plot an IV curve from per-voltage insertion metrics
#'
#' Aggregates current by voltage and condition group, overlays an optional
#' linear conductance fit, and prints conductance statistics to the console.
#'
#' @param input_df       data.frame   Requires: per_voltage_current, voltage.
#'                                    Optional: cond_group, id.
#' @param colour_by      character    Name of the column to colour points by.
#'                                    Defaults to "cond_group" if present,
#'                                    otherwise "id". Pass any column name to
#'                                    override, e.g. colour_by = "treatment".
#' @param show_fit       logical      TRUE (default) overlays a linear fit.
#' @param fit_formula    character    lm formula string (default: "y~x+0").
#' @param fit_voltage    numeric vec  Voltage range used for the fit
#'                                    (default: c(-100, 100)).
#' @param avg_type       character    "mode" (default), "mean", or "median".
#' @param plot_groups    character/NULL  Restrict to these colour_by values.
#' @param cols           character/NULL  Colours per group; ggplot defaults if NULL.
#' @param plot_xlim      numeric/NULL    X-axis limits.
#' @param plot_ylim      numeric/NULL    Y-axis limits.
#' @param error_type     character    "sem" (default) or "stdev".
#' @param marker_size    numeric      Point and line width.
#' @param offset_zero    logical      TRUE shifts current so V = 0 -> I = 0.
#' @param txt_size       numeric      Base font size (default: 7).
#' @param fit_colour     character/NULL  Fit line colour. NULL -> coloured by group.
#' @param drop_zero_v    logical      TRUE (default) drops voltage = 0 rows.
#' @param output         character    "plot" returns ggplot; "data" returns the
#'                                    filtered plot data.frame; "agg_data" returns
#'                                    a list of plot_df, rep_df, and err_df for
#'                                    debugging aggregation (default: "plot").
#'
#' @return ggplot object, data.frame, or named list depending on \code{output}.
#'
#' @examples
#' \dontrun{
#' ivPlot(df)                              # colour by cond_group (default)
#' ivPlot(df, colour_by = "id")            # colour by id column
#' ivPlot(df, colour_by = "treatment")     # colour by any other column
#' ivPlot(df, colour_by = "voltage")       # colour by voltage (discrete)
#' ivPlot(df, output = "agg_data")         # return aggregation internals for QC
#' }
ivPlot <- function(input_df,
                   colour_by    = NULL,
                   show_fit     = TRUE,
                   fit_formula  = "y~x+0",
                   fit_voltage  = c(-100, 100),
                   avg_type     = "mode",
                   plot_groups  = NULL,
                   cols         = NULL,
                   plot_xlim    = NULL,
                   plot_ylim    = NULL,
                   error_type   = "sem",
                   marker_size  = 1,
                   offset_zero  = FALSE,
                   txt_size     = 7,
                   fit_colour   = TRUE,
                   drop_zero_v  = TRUE,
                   output       = "plot") {
  
  # --- Ensure cond_group exists for backwards compatibility ------------------
  if (!"cond_group" %in% names(input_df))
    input_df$cond_group <- if ("id" %in% names(input_df)) input_df$id else "ungrouped"
  
  # --- Resolve colour column -------------------------------------------------
  # Priority: explicit colour_by arg > cond_group > id > "ungrouped"
  colour_col <- .resolve_colour_col(input_df, colour_by,
                                    default_order = c("cond_group", "id"))
  message("ivPlot: colouring by '", colour_col, "'")
  
  # Validate plot_groups against the resolved colour column
  if (is.null(plot_groups)) {
    plot_groups <- unique(input_df[[colour_col]])
  } else {
    missing_groups <- setdiff(plot_groups, unique(input_df[[colour_col]]))
    if (length(missing_groups) > 0L)
      warning("ivPlot: plot_groups values not found in '", colour_col,
              "': ", paste(missing_groups, collapse = ", "))
  }
  
  # --- Orbit device compat: rename current column if needed ------------------
  if (!"per_voltage_current" %in% names(input_df))
    input_df$per_voltage_current <- input_df$current
  
  input_df$per_voltage_current <- as.numeric(input_df$per_voltage_current)
  input_df$voltage             <- as.numeric(input_df$voltage)
  
  if (is.null(fit_voltage)) fit_voltage <- range(input_df$voltage)
  fit_voltage <- seq(fit_voltage[1L], fit_voltage[2L])
  
  # --- Aggregate current per voltage step ------------------------------------
  # The replicate unit is `cm` (one channel-mux = one pore).
  # Each cm has one current value per voltage step from insertionAnalyser.
  # We want:
  #   - plotted point = mean across cms within a group, per voltage
  #   - error         = sd (or sem) across those cms, per voltage
  #
  # colour_col groups the cms (e.g. "treatment", "cond_group").
  # If colour_col IS "cm", each cm is its own group -> no replicates -> no error.
  
  avg_fn <- if (avg_type == "mode") modalAvg else get(avg_type)
  
  # Identify the replicate column: prefer "cm", fall back to "id", else NULL
  rep_col <- if ("cm" %in% names(input_df)) "cm" else
    if ("id" %in% names(input_df) && colour_col != "id") "id" else
      NULL
  
  if (!is.null(rep_col) && colour_col != rep_col) {
    # Step 1: per-replicate average within cm x colour_col x voltage
    # (handles case where input has multiple raw rows per cm per voltage)
    rep_formula <- reformulate(c(rep_col, colour_col, "voltage"),
                               "per_voltage_current")
    rep_df <- stats::aggregate(rep_formula, input_df, avg_fn)
  } else {
    # colour_by IS the replicate column - no within-group replicates exist
    message("ivPlot: colour_by = '", colour_col, "' is the replicate column — ",
            "error bars will be NA (no replicates within each group).")
    rep_df <- input_df[, intersect(c(colour_col, "voltage", "per_voltage_current"),
                                   names(input_df))]
  }
  
  if (isTRUE(offset_zero) && any(rep_df$voltage == 0)) {
    ## zero offsetting in per colour_col manner
    zero_formula <- reformulate(colour_col, "per_voltage_current")
    
    zero_df <- stats::aggregate(
      zero_formula,
      data = rep_df[rep_df$voltage == 0, , drop = FALSE],
      FUN  = mean
    )
    
    names(zero_df)[names(zero_df) == "per_voltage_current"] <- "zero_current"
    
    # Merge offsets back onto rep_df
    rep_df <- merge(rep_df, zero_df,
                    by = colour_col,
                    all.x = TRUE)
    
    # Apply group-wise offset
    rep_df$per_voltage_current <-
      rep_df$per_voltage_current - rep_df$zero_current
    
    # Cleanup helper column
    rep_df$zero_current <- NULL
  }
  
  # Step 2: compute sd, n, and sem atomically in a single aggregate call so
  # all three statistics are guaranteed to reflect the same set of valid
  # observations per voltage x colour_col cell.
  err_formula <- reformulate(c("voltage", colour_col), "per_voltage_current")
  
  err_df <- stats::aggregate(err_formula, rep_df, function(x) {
    valid <- x[!is.na(x)]
    n     <- length(valid)
    stdev <- sd(valid)
    c(stdev = stdev, n = n, sem = stdev / sqrt(n))
  })
  
  # Unpack the matrix column produced by aggregate's multi-value FUN
  err_df <- do.call(data.frame, list(
    err_df[, c("voltage", colour_col)],
    stdev = err_df$per_voltage_current[, "stdev"],
    n     = err_df$per_voltage_current[, "n"],
    sem   = err_df$per_voltage_current[, "sem"]
  ))
  
  # Step 3: mean of per-replicate values -> one plotted point per voltage x group
  plot_formula <- reformulate(c("voltage", colour_col), "per_voltage_current")
  plot_df      <- stats::aggregate(plot_formula, rep_df, mean)
  plot_df      <- merge(plot_df, err_df, by = c("voltage", colour_col))
  
  # Diagnostic: report per-voltage n range per group, flag unbalanced designs
  message("ivPlot: n replicates per voltage x group (range per group):")
  for (grp in unique(err_df[[colour_col]])) {
    ns <- err_df$n[err_df[[colour_col]] == grp]
    message("  ", grp, ": ", min(ns), "-", max(ns),
            if (min(ns) != max(ns)) " (UNBALANCED - check for missing voltage steps)" else "")
  }
  
  # --- Filter to selected groups ---------------------------------------------
  ap_df <- plot_df[plot_df[[colour_col]] %in% plot_groups, , drop = FALSE]
  
  if (isTRUE(drop_zero_v)) {
    message("ivPlot: dropping voltage = 0 rows (note: affects offset_zero correction)")
    ap_df <- ap_df[ap_df$voltage != 0, , drop = FALSE]
  }
  
  # Return aggregation internals for QC without building the plot
  if (output == "agg_data") {
    return(list(plot_df = plot_df, rep_df = rep_df, err_df = err_df))
  }
  
  if (output == "data") return(ap_df)
  
  err_col <- if (tolower(error_type) == "stdev") "stdev" else "sem"
  
  # --- Build plot ------------------------------------------------------------
  plt <- ggplot2::ggplot(
    ap_df,
    ggplot2::aes(
      x      = voltage,
      y      = per_voltage_current,
      colour = .data[[colour_col]],
      fill   = .data[[colour_col]],
      ymin   = per_voltage_current - .data[[err_col]],
      ymax   = per_voltage_current + .data[[err_col]]
    )
  ) +
    ggplot2::geom_point(size = marker_size) +
    ggplot2::geom_errorbar(width = 2, alpha = 0.5, linewidth = marker_size) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey",
                        linetype = "dashed", linewidth = marker_size) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey",
                        linetype = "dashed", linewidth = marker_size) +
    ggplot2::xlab("Voltage (mV)") +
    ggplot2::ylab("Current (pA)") +
    My_Theme(base_size = txt_size, txt_size = txt_size) +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  
  # --- Conductance fit -------------------------------------------------------
  plot_df$cond <- plot_df$per_voltage_current / plot_df$voltage
  
  if (isTRUE(show_fit)) {
    for (grp in unique(plot_df[[colour_col]])) {
      tmp   <- stats::na.omit(plot_df[plot_df[[colour_col]] == grp, ])
      tmp   <- tmp[!is.infinite(tmp$cond), ]
      tmp$y <- tmp$per_voltage_current
      tmp$x <- tmp$voltage
      model <- try(stats::lm(
        formula = as.formula(fit_formula),
        data    = tmp[tmp$voltage %in% fit_voltage, ]
      ), silent = TRUE)
      
      message("ivPlot: conductance for '", colour_col, " = ", grp, "'")
      if (!inherits(model, "try-error")) {
        message("  Linear fit:")
        print(summary(model)$coefficients)
        message("  Modal average: ", round(modalAvg(tmp$cond), 4),
                " \u00b1 ", round(sd(tmp$cond), 4), " nS")
      } else {
        message("  Fit failed for this group.")
      }
    }
    
    fit_df <- ap_df[ap_df$voltage %in% fit_voltage, , drop = FALSE]
    
    if (isTRUE(fit_colour)) {
      plt <- plt + ggplot2::geom_smooth(
        data      = fit_df,
        ggplot2::aes(colour = .data[[colour_col]]),
        method    = "lm",
        formula   = as.formula(fit_formula),
        se        = FALSE,
        linewidth = marker_size
      )
    } else {
      plt <- plt + ggplot2::geom_smooth(
        data      = fit_df,
        method    = "lm",
        formula   = as.formula(fit_formula),
        se        = FALSE,
        linewidth = marker_size,
        colour    = fit_colour
      )
    }
  }
  
  if (!is.null(cols))      plt <- plt + ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols)
  if (!is.null(plot_ylim)) plt <- plt + ggplot2::scale_y_continuous(
    breaks = stefBreaks(plot_ylim)[[1L]],
    labels = stefBreaks(plot_ylim)[[2L]], limits = plot_ylim)
  if (!is.null(plot_xlim)) plt <- plt + ggplot2::scale_x_continuous(
    breaks = stefBreaks(plot_xlim)[[1L]],
    labels = stefBreaks(plot_xlim)[[2L]], limits = plot_xlim)
  
  return(plt)
}



# =============================================================================

#' Plot a conductance histogram
#'
#' @param input_df    data.frame   Requires: per_voltage_cond or per_chan_cond,
#'                                 voltage, cm. Optional: cond_group, id.
#' @param colour_by   character    Name of the column to fill bars by.
#'                                 Defaults to "cond_group" if present,
#'                                 otherwise "id". Pass any column name to
#'                                 override, e.g. colour_by = "salt_conc".
#' @param n_bins      integer      Number of histogram bins (default: 10).
#' @param his_voltage numeric/NULL Restrict to this voltage before plotting.
#' @param per_channel logical      TRUE uses per_chan_cond (averaged across all
#'                                 voltages). FALSE uses per_voltage_cond at
#'                                 his_voltage.
#' @param plot_groups character/NULL  Restrict to these colour_by values.
#' @param cols        character/NULL  Fill colours per group.
#' @param plot_ylim   numeric/NULL    Y-axis limits.
#' @param plot_xlim   numeric/NULL    X-axis limits.
#' @param show_n      logical      TRUE annotates plot with pore count n.
#' @param n_x         numeric/NULL X position of n label (auto if NULL).
#' @param n_y         numeric/NULL Y position of n label (auto if NULL).
#' @param txt_size    numeric      Base font size (default: 7).
#'
#' @return ggplot object.
#'
#' @examples
#' \dontrun{
#' hisPlot(df)                             # fill by cond_group (default)
#' hisPlot(df, colour_by = "id")           # fill by id
#' hisPlot(df, colour_by = "salt_conc")    # fill by any column
#' }
hisPlot <- function(input_df,
                    colour_by   = NULL,
                    n_bins      = 10,
                    his_voltage = NULL,
                    per_channel = TRUE,
                    plot_groups = NULL,
                    cols        = NULL,
                    plot_ylim   = NULL,
                    plot_xlim   = NULL,
                    show_n      = TRUE,
                    n_x         = NULL,
                    n_y         = NULL,
                    txt_size    = 7) {
  
  # --- Resolve colour column -------------------------------------------------
  colour_col <- .resolve_colour_col(input_df, colour_by,
                                    default_order = c("cond_group", "id"))
  message("hisPlot: colouring by '", colour_col, "'")
  
  if (!is.null(plot_groups))
    input_df <- input_df[input_df[[colour_col]] %in% plot_groups, , drop = FALSE]
  
  his_df <- if (!is.null(his_voltage)) {
    input_df[input_df$voltage %in% his_voltage, , drop = FALSE]
  } else {
    input_df
  }
  
  if (!isTRUE(per_channel)) {
    message("hisPlot: using per-voltage conductance — ",
            "set his_voltage to avoid inflated n.")
    his_df$per_chan_cond <- his_df$per_voltage_cond
  }
  
  his_df$per_chan_cond <- abs(as.numeric(his_df$per_chan_cond))
  his_df               <- his_df[his_df$per_chan_cond > 1e-5, , drop = FALSE]
  
  if (nrow(his_df) == 0L) {
    warning("hisPlot: no data remaining after filters.")
    return(NULL)
  }
  
  # Recycle colours if fewer supplied than groups
  n_groups <- length(unique(his_df[[colour_col]]))
  if (!is.null(cols) && length(cols) < n_groups)
    cols <- rep(cols, length.out = n_groups)
  
  # --- Build plot ------------------------------------------------------------
  plt <- ggplot2::ggplot(
    his_df,
    ggplot2::aes(x = per_chan_cond, fill = .data[[colour_col]])
  ) +
    ggplot2::geom_histogram(bins = n_bins) +
    ggplot2::ylab("Counts") +
    ggplot2::xlab("Conductance / nS") +
    ggplot2::theme_bw(base_size = txt_size) +
    My_Theme(base_size = txt_size, txt_size = txt_size) +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  
  if (!is.null(cols))      plt <- plt + ggplot2::scale_fill_manual(values = cols)
  if (!is.null(plot_ylim)) plt <- plt + ggplot2::scale_y_continuous(
    breaks = stefBreaks(plot_ylim)[[1L]],
    labels = stefBreaks(plot_ylim)[[2L]], limits = plot_ylim)
  if (!is.null(plot_xlim)) plt <- plt + ggplot2::scale_x_continuous(
    breaks = stefBreaks(plot_xlim)[[1L]],
    labels = stefBreaks(plot_xlim)[[2L]], limits = plot_xlim)
  
  if (isTRUE(show_n)) {
    n   <- length(unique(his_df$cm))
    n_x <- n_x %||% (max(his_df$per_chan_cond, na.rm = TRUE) * 0.75)
    n_y <- n_y %||% (n / 10)
    plt <- plt + ggplot2::annotate(
      "text", x = n_x, y = n_y,
      label = paste0("n = ", n),
      size  = txt_size / ggplot2::.pt
    )
  }
  
  plt
}


# =============================================================================




#' Fit probability distributions to a conductance histogram
#'
#' Adds fitted distribution curves to a \code{hisPlot} output and returns the
#' estimated peak conductance per group. Supports automatic model selection
#' by AIC, a single distribution applied to all groups, or different
#' distributions per group.
#'
#' The grouping column is resolved automatically from the plot's fill mapping
#' when \code{group_col = NULL} (the default), so you don't need to repeat
#' the \code{colour_by} argument from \code{hisPlot}. Override with an explicit
#' column name to split on a different variable entirely.
#'
#' @param his_plot     ggplot   Output of \code{hisPlot} (or any ggplot histogram
#'                              with a \code{per_chan_cond} column in its data).
#' @param distribution character or named character vector  Distribution(s) to
#'   fit. Options: "auto", "normal", "gamma", "poisson", "exponential",
#'   "lognormal", "beta", "weibull", "logistic", "chisq", "geom".
#'   Pass a named vector to assign different distributions per group, where
#'   names must match the values of the resolved grouping column.
#' @param line_size    numeric   Fitted curve line thickness (default: 0.5).
#' @param line_colours character/NULL  Colours for the fitted curves, one per
#'   group, in the same order as the groups appear in the data. NULL (default)
#'   draws all curves in black. Recycled if fewer colours than groups are
#'   supplied. Pass a named vector to assign colours explicitly by group id,
#'   e.g. \code{line_colours = c(ctrl = "blue", treated = "red")}.
#' @param group_col    character/NULL  Name of the column in the histogram data
#'   to split fits by. NULL (default) auto-detects from the plot's fill
#'   aesthetic mapping. Pass an explicit column name to override, e.g. to fit
#'   by a different grouping than what was used for colours.
#'
#' @return list with two elements:
#'   \item{plot}{ggplot with fitted distribution curves overlaid.}
#'   \item{peaks}{data.frame with columns: id, distribution, peak, error.}
#'
#' @examples
#' \dontrun{
#' p <- hisPlot(df, colour_by = "treatment")
#'
#' # Auto-detects "treatment" from the plot mapping — no need to repeat it
#' fit <- hisFitr(p)
#'
#' # Same distribution for all groups
#' fit <- hisFitr(p, distribution = "gamma")
#'
#' # Different distribution per group
#' fit <- hisFitr(p, distribution = c(ctrl = "normal", treated = "gamma"))
#'
#' # Override to split by a different column than the fill aesthetic
#' fit <- hisFitr(p, group_col = "salt_conc")
#' }
#' 
# =============================================================================
# Internal helpers — not exported
# =============================================================================

#' Detect facet variable names from a built ggplot object.
#' Returns character(0) for non-faceted plots.
.get_facet_vars <- function(gg_build) {
  facet <- gg_build$plot$facet
  if (inherits(facet, "FacetWrap")) {
    names(facet$params$facets)
  } else if (inherits(facet, "FacetGrid")) {
    c(names(facet$params$rows), names(facet$params$cols))
  } else {
    character(0L)
  }
}

#' Map facet variable values in a data slice to the PANEL integer used by
#' ggplot_build, then return the max count for that panel.
#' Falls back to the global max if the panel cannot be identified.
.panel_hist_peak <- function(gg_build, fdata, facet_vars) {
  built_df <- gg_build$data[[1L]]
  
  if (length(facet_vars) == 0L || !"PANEL" %in% names(built_df)) {
    return(max(built_df$count, na.rm = TRUE))
  }
  
  layout <- gg_build$layout$layout
  
  key <- unique(fdata[, facet_vars, drop = FALSE])
  if (nrow(key) > 1L) key <- key[1L, , drop = FALSE]
  
  match_rows <- Reduce(`&`, lapply(facet_vars, function(v) {
    layout[[v]] == key[[v]]
  }))
  
  if (!any(match_rows, na.rm = TRUE)) {
    return(max(built_df$count, na.rm = TRUE))
  }
  
  panel_id   <- layout$PANEL[match_rows][1L]
  panel_rows <- built_df[built_df$PANEL == panel_id, ]
  max(panel_rows$count, na.rm = TRUE)
}

#' Resolve a fit_range for a given group/facet combination.
#'
#' fit_range can be:
#'   NULL                         — no filtering
#'   c(lo, hi)                    — applied to all groups and facets
#'   list(g1 = c(lo, hi), ...)    — per-group; each value can itself be a
#'                                  named list keyed by facet key strings for
#'                                  per-facet control, or a plain c(lo, hi)
#'
#' facet_id is the pasted facet key string (or the group id for non-faceted
#' plots). If no per-facet entry is found, falls back to the group-level range.
.resolve_fit_range <- function(fit_range, gid, facet_id) {
  if (is.null(fit_range)) return(NULL)
  
  # Plain numeric vector — applies everywhere
  if (is.numeric(fit_range)) {
    if (length(fit_range) != 2L)
      stop("hisFitr: fit_range must be a length-2 numeric vector c(lo, hi), ",
           "a named list by group, or NULL.")
    return(fit_range)
  }
  
  if (!is.list(fit_range))
    stop("hisFitr: fit_range must be a length-2 numeric vector, a named list, or NULL.")
  
  # Named list — look up by group
  if (!gid %in% names(fit_range)) return(NULL)   # no range for this group
  
  grp_range <- fit_range[[gid]]
  
  # Group entry is itself a named list keyed by facet
  if (is.list(grp_range)) {
    if (facet_id %in% names(grp_range)) {
      r <- grp_range[[facet_id]]
    } else {
      return(NULL)   # no range for this facet within the group
    }
  } else {
    r <- grp_range   # plain c(lo, hi) for the whole group
  }
  
  if (!is.numeric(r) || length(r) != 2L)
    stop("hisFitr: fit_range entry for group '", gid, "' must be c(lo, hi).")
  r
}




# =============================================================================
# Main function
# =============================================================================

#' Fit probability distributions to a conductance histogram
#'
#' Adds fitted distribution curves to a \code{hisPlot} output and returns the
#' estimated peak conductance per group. Supports automatic model selection
#' by AIC, a single distribution applied to all groups, or different
#' distributions per group.
#'
#' @param his_plot     ggplot   Output of \code{hisPlot}.
#' @param distribution character or named character vector. Options: "auto",
#'   "normal", "gamma", "poisson", "exponential", "lognormal", "beta",
#'   "weibull", "logistic", "chisq", "geom". Pass a named vector to assign
#'   different distributions per group.
#' @param fit_range    NULL, numeric c(lo, hi), or named list. Restricts the
#'   data used for fitting to rows where \code{per_chan_cond} falls within
#'   \code{[lo, hi]}. The fitted curve is still drawn over the full axis.
#'
#'   Three levels of control:
#'   \itemize{
#'     \item \code{NULL} (default) — fit all data.
#'     \item \code{c(lo, hi)} — same range applied to every group and facet.
#'     \item Named list — per-group ranges, e.g.
#'       \code{list(ctrl = c(0, 50), treated = c(10, 80))}.
#'       Within a group entry you can go one level deeper to specify ranges
#'       per facet, using the pasted facet key as the name, e.g.
#'       \code{list(ctrl = list("drugA" = c(0, 40), "drugB" = c(5, 60)))}.
#'       Groups or facets without a matching entry are fitted on all their data.
#'   }
#' @param line_size    numeric   Line thickness (default 0.5).
#' @param line_colours character/NULL  Colours per group. NULL = all black.
#'   Recycled if fewer colours than groups. Pass a named vector to assign by
#'   group id.
#' @param group_col    character/NULL  Column to split fits by. NULL
#'   auto-detects from the plot's fill aesthetic.
#' @param min_n        integer   Minimum observations required after
#'   \code{fit_range} filtering before a fit is attempted (default 10). A
#'   descriptive error is thrown if any group/facet subset falls below this.
#'
#' @return list: \code{plot} (ggplot with curves overlaid) and \code{peaks}
#'   (data.frame with columns: id, [facet,] distribution, peak, error,
#'   [fit_range_lo, fit_range_hi]).
hisFitr <- function(his_plot,
                    distribution = "auto",
                    fit_range    = NULL,
                    line_size    = 0.5,
                    line_colours = NULL,
                    group_col    = NULL,
                    min_n        = 10L) {
  
  available <- c("auto", "normal", "gamma", "poisson", "exponential",
                 "lognormal", "beta", "weibull", "logistic", "chisq", "geom")
  message("hisFitr: available distributions: ", paste(available, collapse = ", "))
  
  # ---------------------------------------------------------------------------
  # 1. Extract data and plot metadata
  # ---------------------------------------------------------------------------
  gg_build   <- ggplot2::ggplot_build(his_plot)
  his_data   <- stats::na.omit(gg_build$plot$data)
  facet_vars <- .get_facet_vars(gg_build)
  
  if (length(facet_vars) > 0L)
    message("hisFitr: faceted plot detected — facet variable(s): ",
            paste(facet_vars, collapse = ", "))
  
  # ---------------------------------------------------------------------------
  # 2. Resolve grouping column
  # ---------------------------------------------------------------------------
  if (is.null(group_col)) {
    fill_mapping <- gg_build$plot$mapping$fill
    if (!is.null(fill_mapping)) {
      detected <- sub("^~", "", deparse(fill_mapping))
      detected <- gsub('^\\.data\\[\\["(.*)"\\]\\]$', "\\1", detected)
      if (detected %in% colnames(his_data)) {
        group_col <- detected
        message("hisFitr: auto-detected grouping column '", group_col,
                "' from plot fill mapping.")
      }
    }
    if (is.null(group_col)) {
      for (col in c("cond_group", "id")) {
        if (col %in% colnames(his_data)) { group_col <- col; break }
      }
    }
    if (is.null(group_col)) {
      message("hisFitr: no grouping column detected — fitting all data as one group.")
      his_data[["..all.."]] <- "all"
      group_col <- "..all.."
    }
  } else {
    if (!group_col %in% colnames(his_data))
      stop("hisFitr: group_col '", group_col, "' not found in histogram data.\n",
           "  Available columns: ", paste(colnames(his_data), collapse = ", "), "\n",
           "  Leave group_col = NULL to auto-detect from the plot mapping.")
  }
  
  message("hisFitr: splitting fits by column '", group_col, "'")
  
  # ---------------------------------------------------------------------------
  # 3. Colour lookup
  # ---------------------------------------------------------------------------
  x_range <- ggplot2::layer_scales(his_plot)$x$range$range
  groups  <- split(his_data, his_data[[group_col]])
  
  group_ids     <- names(groups)
  n_groups      <- length(group_ids)
  colour_lookup <- if (is.null(line_colours)) {
    stats::setNames(rep("black", n_groups), group_ids)
  } else if (!is.null(names(line_colours))) {
    missing_ids <- setdiff(group_ids, names(line_colours))
    if (length(missing_ids) > 0L)
      warning("hisFitr: no line colour specified for group(s): ",
              paste(missing_ids, collapse = ", "), " — using black.")
    defaults <- stats::setNames(rep("black", n_groups), group_ids)
    defaults[names(line_colours)] <- line_colours
    defaults
  } else {
    stats::setNames(rep(line_colours, length.out = n_groups), group_ids)
  }
  
  # ---------------------------------------------------------------------------
  # 4. AIC-based distribution selection
  # ---------------------------------------------------------------------------
  select_best_dist <- function(x) {
    candidates <- setdiff(available, "auto")
    aics <- vapply(candidates, function(d) {
      tryCatch({
        f <- MASS::fitdistr(x, d)
        2 * length(f$estimate) - 2 * f$loglik
      }, error = function(e) Inf)
    }, numeric(1L))
    
    if (all(is.infinite(aics)))
      stop("hisFitr: no distributions could be fitted.")
    
    best <- candidates[which.min(aics)]
    message("    Auto-selected: ", best,
            " (AIC = ", round(min(aics[is.finite(aics)]), 2L), ")")
    list(fit = MASS::fitdistr(x, best), dist = best)
  }
  
  # ---------------------------------------------------------------------------
  # 5. Density function builder
  # ---------------------------------------------------------------------------
  make_dens_fn <- function(dist_used, est) {
    switch(dist_used,
           normal      = function(x) stats::dnorm(x,  est["mean"],     est["sd"]),
           gamma       = function(x) stats::dgamma(x, est["shape"],    est["rate"]),
           poisson     = function(x) stats::dpois(round(x),            est["lambda"]),
           exponential = function(x) stats::dexp(x,                    est["rate"]),
           lognormal   = function(x) stats::dlnorm(x, est["meanlog"],  est["sdlog"]),
           beta        = function(x) stats::dbeta(x,  est["shape1"],   est["shape2"]),
           weibull     = function(x) stats::dweibull(x, est["shape"],  est["scale"]),
           logistic    = function(x) stats::dlogis(x, est["location"], est["scale"]),
           chisq       = function(x) stats::dchisq(x,                  est["df"]),
           geom        = function(x) stats::dgeom(round(x),            est["prob"])
    )
  }
  
  # ---------------------------------------------------------------------------
  # 6. Analytic peak helper
  # ---------------------------------------------------------------------------
  analytic_peak <- function(dist_used, est, dens_fn, x_range) {
    peak_g <- switch(dist_used,
                     normal      = est["mean"],
                     gamma       = if (est["shape"] > 1) (est["shape"] - 1) / est["rate"] else NA,
                     exponential = 0,
                     lognormal   = exp(est["meanlog"] - est["sdlog"]^2),
                     weibull     = if (est["shape"] > 1)
                       est["scale"] * ((est["shape"] - 1) / est["shape"])^(1 / est["shape"])
                     else 0,
                     chisq       = if (est["df"] >= 2) est["df"] - 2 else 0,
                     logistic    = est["location"],
                     NA_real_
    )
    
    if (is.na(peak_g)) {
      opt    <- stats::optimize(function(x) -dens_fn(x), interval = range(x_range))
      peak_g <- opt$minimum
      message("    Peak estimated numerically: ", round(peak_g, 4L))
    } else {
      message("    Peak (analytic): ", round(peak_g, 4L))
    }
    peak_g
  }
  
  # ---------------------------------------------------------------------------
  # 7. Main fitting loop — outer: group, inner: facet subset
  # ---------------------------------------------------------------------------
  all_results <- vector("list", length(groups))
  
  for (i in seq_along(groups)) {
    gdata <- groups[[i]]
    gid   <- names(groups)[i]
    
    dist_req <- if (length(distribution) > 1L) {
      if (!gid %in% names(distribution))
        stop("hisFitr: no distribution specified for group '", gid,
             "'. Provided names: ", paste(names(distribution), collapse = ", "))
      as.character(distribution[[gid]])
    } else {
      as.character(distribution)
    }
    
    message("hisFitr: group '", gid, "' — ", dist_req)
    
    # Split by facet (if any)
    if (length(facet_vars) > 0L) {
      facet_key     <- do.call(paste, c(gdata[facet_vars], sep = "__"))
      facet_subsets <- split(gdata, facet_key)
    } else {
      facet_subsets        <- list(gdata)
      names(facet_subsets) <- gid
    }
    
    sub_results <- vector("list", length(facet_subsets))
    
    for (j in seq_along(facet_subsets)) {
      fdata    <- facet_subsets[[j]]
      facet_id <- names(facet_subsets)[j]
      
      message("  facet subset: '", facet_id, "' — n = ", nrow(fdata))
      
      # -----------------------------------------------------------------------
      # Apply fit_range filter to the fitting data only
      # -----------------------------------------------------------------------
      rng <- .resolve_fit_range(fit_range, gid, facet_id)
      
      if (!is.null(rng)) {
        fdata_fit <- fdata[fdata$per_chan_cond >= rng[1L] &
                             fdata$per_chan_cond <= rng[2L], ]
        n_fit <- nrow(fdata_fit)
        message("    fit_range [", rng[1L], ", ", rng[2L], "] — ",
                n_fit, " / ", nrow(fdata), " observations retained for fitting.")
        
        if (n_fit < min_n)
          stop(
            "hisFitr: insufficient data to fit.\n",
            "  Group    : '", gid, "'\n",
            "  Facet    : '", facet_id, "'\n",
            "  fit_range: [", rng[1L], ", ", rng[2L], "]\n",
            "  In range : ", n_fit, " observation(s)\n",
            "  Total    : ", nrow(fdata), " observation(s) in this subset\n",
            "  Required : min_n = ", min_n, "\n",
            "  → Widen fit_range, or set min_n to a lower value if a fit ",
            "is still meaningful with fewer observations."
          )
      } else {
        fdata_fit <- fdata
      }
      
      # -----------------------------------------------------------------------
      # Fit distribution (on filtered data only)
      # -----------------------------------------------------------------------
      if (dist_req == "auto") {
        res       <- select_best_dist(fdata_fit$per_chan_cond)
        fit       <- res$fit
        dist_used <- res$dist
      } else {
        dist_used <- dist_req
        fit <- tryCatch(
          MASS::fitdistr(fdata_fit$per_chan_cond, dist_used),
          error = function(e)
            stop("hisFitr: fitdistr failed for group '", gid,
                 "', facet '", facet_id, "': ", e$message)
        )
      }
      
      # -----------------------------------------------------------------------
      # Build density curve over the FULL axis range (not just fit_range)
      # -----------------------------------------------------------------------
      x_vals    <- seq(min(x_range), max(x_range), length.out = 1000L)
      est       <- fit$estimate
      dens_fn   <- make_dens_fn(dist_used, est)
      raw_dens  <- dens_fn(x_vals)
      
      # Scale to per-panel histogram peak using unfiltered fdata
      panel_peak  <- .panel_hist_peak(gg_build, fdata, facet_vars)
      scaled_dens <- raw_dens * (panel_peak / max(raw_dens, na.rm = TRUE))
      
      peak_g <- analytic_peak(dist_used, est, dens_fn, x_range)
      
      # curve_df carries facet variable columns so ggplot routes correctly
      curve_df      <- data.frame(x = x_vals, y = scaled_dens, id = gid)
      curve_df$y[curve_df$y <= 0] <- NA
      
      for (fv in facet_vars) {
        curve_df[[fv]] <- fdata[[fv]][1L]
      }
      
      his_plot <- his_plot +
        ggplot2::geom_line(
          data        = curve_df,
          ggplot2::aes(x = x, y = y, group = id),
          linewidth   = line_size,
          colour      = colour_lookup[[gid]],
          inherit.aes = FALSE
        )
      
      sub_results[[j]] <- data.frame(
        id           = gid,
        facet        = facet_id,
        distribution = dist_used,
        n_hist       = nrow(fdata),
        peak         = as.numeric(peak_g),
        error        = min(fit$sd, na.rm = TRUE),
        fit_range_lo = if (!is.null(rng)) rng[1L] else NA_real_,
        fit_range_hi = if (!is.null(rng)) rng[2L] else NA_real_
      )
    }
    
    all_results[[i]] <- do.call(rbind, sub_results)
  }
  
  peaks_df <- do.call(rbind, all_results)
  
  # Drop columns that carry no information in the non-faceted / no-range case
  if (length(facet_vars) == 0L)     peaks_df$facet        <- NULL
  if (all(is.na(peaks_df$fit_range_lo))) {
    peaks_df$fit_range_lo <- NULL
    peaks_df$fit_range_hi <- NULL
  }
  
  print(peaks_df)
  list(plot = his_plot, peaks = peaks_df)
}

### legacy code left in for posterity

#' 
#' #' Plot an IV curve from per-voltage insertion metrics
#' #'
#' #' Aggregates current by voltage and condition group, overlays an optional
#' #' linear conductance fit, and prints conductance statistics to the console.
#' #'
#' #' @param input_df       data.frame  Requires: per_voltage_current, voltage.
#' #'                         Optional: cond_group, id.
#' #' @param show_fit       logical  TRUE (default) overlays a linear fit.
#' #' @param fit_formula    character  lm formula string (default: "y~x+0").
#' #' @param fit_voltage    numeric vector  Voltage range for the fit
#' #'                         (default: c(-100, 100)).
#' #' @param avg_type       character  "mode" (default), "mean", or "median".
#' #' @param plot_groups    character/NULL  Restrict to these cond_group values.
#' #' @param cols           character/NULL  Colours per group.
#' #' @param plot_xlim      numeric/NULL  X-axis limits.
#' #' @param plot_ylim      numeric/NULL  Y-axis limits.
#' #' @param error_type     character  "sem" (default) or "stdev".
#' #' @param marker_size    numeric  Point and line width.
#' #' @param offset_zero    logical  TRUE shifts current so V = 0 → I = 0.
#' #' @param txt_size       numeric  Base font size.
#' #' @param fit_colour     character/NULL  Fit line colour. NULL → coloured by group.
#' #' @param drop_zero_v    logical  TRUE (default) drops voltage = 0 rows.
#' #' @param output         character  "plot" returns ggplot; "data" returns the
#' #'                         aggregated data frame.
#' #'
#' #' @return ggplot object or data.frame depending on `output`.
#' ivPlot <- function(input_df,
#'                    show_fit     = TRUE,
#'                    fit_formula  = "y~x+0",
#'                    fit_voltage  = c(-100, 100),
#'                    avg_type     = "mode",
#'                    plot_groups  = NULL,
#'                    cols         = NULL,
#'                    plot_xlim    = NULL,
#'                    plot_ylim    = NULL,
#'                    error_type   = "sem",
#'                    marker_size  = 1,
#'                    offset_zero  = FALSE,
#'                    txt_size     = 7,
#'                    fit_colour   = "black",
#'                    drop_zero_v  = TRUE,
#'                    output       = "plot") {
#'   
#'   if (!"cond_group" %in% names(input_df))
#'     input_df$cond_group <- if ("id" %in% names(input_df)) input_df$id else "ungrouped"
#'   
#'   if (!"per_voltage_current" %in% names(input_df))
#'     input_df$per_voltage_current <- input_df$current
#'   
#'   input_df$per_voltage_current <- as.numeric(input_df$per_voltage_current)
#'   input_df$voltage             <- as.numeric(input_df$voltage)
#'   
#'   if (is.null(plot_groups)) plot_groups <- unique(input_df$cond_group)
#'   if (is.null(fit_voltage)) fit_voltage <- range(input_df$voltage)
#'   fit_voltage <- seq(fit_voltage[1], fit_voltage[2])
#'   
#'   # ── Aggregate ────────────────────────────────────────────────────────────────
#'   grp_vars <- if ("id" %in% names(input_df)) c("cond_group", "voltage", "id")
#'   else c("cond_group", "voltage")
#'   
#'   avg_fn  <- if (avg_type == "mode") modalAvg else get(avg_type)
#'   # plot_df <- stats::aggregate(
#'   #   per_voltage_current ~ .,
#'   #   data = input_df[, c(grp_vars, "per_voltage_current")],
#'   #   FUN  = avg_fn
#'   # )
#'   plot_df <- input_df[
#'     ,
#'     .(per_voltage_current = avg_fn(per_voltage_current)),
#'     by = grp_vars
#'   ]
#'   # ── Error bars ───────────────────────────────────────────────────────────────
#'   sd_df   <- stats::aggregate(per_voltage_current ~ voltage + cond_group, input_df, sd)
#'   n_df    <- stats::aggregate(per_voltage_current ~ voltage + cond_group, input_df, length)
#'   err_df  <- merge(sd_df, n_df, by = c("voltage", "cond_group"))
#'   names(err_df)[3:4] <- c("stdev", "n")
#'   err_df$sem <- err_df$stdev / sqrt(err_df$n)
#'   plot_df    <- merge(plot_df, err_df, by = c("voltage", "cond_group"))
#'   
#'   ap_df <- subset(plot_df, plot_df$cond_group %in% plot_groups)
#'   
#'   if (isTRUE(offset_zero))
#'     ap_df$per_voltage_current <- ap_df$per_voltage_current -
#'     ap_df$per_voltage_current[ap_df$voltage == 0]
#'   
#'   if (isTRUE(drop_zero_v)) {
#'     message("ivPlot: dropping voltage = 0 rows (note: affects offset correction)")
#'     ap_df <- subset(ap_df, ap_df$voltage != 0)
#'   }
#'   
#'   err_col <- if (tolower(error_type) == "stdev") "stdev" else "sem"
#'   
#'   plt <- ggplot2::ggplot(
#'     ap_df,
#'     ggplot2::aes(
#'       x = voltage, y = per_voltage_current, colour = cond_group,
#'       ymin = per_voltage_current - .data[[err_col]],
#'       ymax = per_voltage_current + .data[[err_col]]
#'     )
#'   ) +
#'     ggplot2::geom_point(size = marker_size) +
#'     ggplot2::geom_errorbar(width = 2, alpha = 0.5, linewidth = marker_size) +
#'     ggplot2::geom_hline(yintercept = 0, colour = "grey",
#'                         linetype = "dashed", linewidth = marker_size) +
#'     ggplot2::geom_vline(xintercept = 0, colour = "grey",
#'                         linetype = "dashed", linewidth = marker_size) +
#'     ggplot2::xlab("Voltage (mV)") + ggplot2::ylab("Current (pA)") +
#'     My_Theme(base_size = txt_size, txt_size = txt_size) +
#'     ggplot2::theme(legend.title = ggplot2::element_blank())
#'   
#'   # ── Conductance fit ──────────────────────────────────────────────────────────
#'   plot_df$cond <- plot_df$per_voltage_current / plot_df$voltage
#'   
#'   if (isTRUE(show_fit)) {
#'     for (grp in unique(plot_df$cond_group)) {
#'       tmp      <- stats::na.omit(subset(plot_df, plot_df$cond_group == grp))
#'       tmp      <- subset(tmp, !is.infinite(tmp$cond))
#'       tmp$y    <- tmp$per_voltage_current
#'       tmp$x    <- tmp$voltage
#'       model    <- try(stats::lm(
#'         formula = as.formula(fit_formula),
#'         data    = subset(tmp, tmp$voltage %in% fit_voltage)
#'       ))
#'       message("ivPlot: conductance for group '", grp, "'")
#'       if (!inherits(model, "try-error")) {
#'         message("  Linear fit:")
#'         print(summary(model)$coefficients)
#'         message("  Modal average: ", modalAvg(tmp$cond), " ± ", sd(tmp$cond), " nS")
#'       }
#'     }
#'     
#'     fit_df <- subset(ap_df, ap_df$voltage %in% fit_voltage)
#'     plt <- plt + ggplot2::geom_smooth(
#'       data      = fit_df,
#'       ggplot2::aes(colour = cond_group),
#'       method    = "lm",
#'       formula   = as.formula(fit_formula),
#'       se        = FALSE,
#'       linewidth = marker_size,
#'       colour    = fit_colour
#'     )
#'   }
#'   
#'   if (!is.null(cols))      plt <- plt + ggplot2::scale_colour_manual(values = cols)
#'   if (!is.null(plot_ylim)) plt <- plt + ggplot2::scale_y_continuous(
#'     breaks = stefBreaks(plot_ylim)[[1]],
#'     labels = stefBreaks(plot_ylim)[[2]], limits = plot_ylim
#'   )
#'   if (!is.null(plot_xlim)) plt <- plt + ggplot2::scale_x_continuous(
#'     breaks = stefBreaks(plot_xlim)[[1]],
#'     labels = stefBreaks(plot_xlim)[[2]], limits = plot_xlim
#'   )
#'   
#'   if (output == "plot") return(plt)
#'   if (output == "data") return(ap_df)
#'   stop("ivPlot: output must be 'plot' or 'data'.")
#' }
#' 
#' 
#' #' Plot a conductance histogram
#' #'
#' #' @param input_df    data.frame  Requires: per_voltage_cond or per_chan_cond,
#' #'                     voltage, cm. Optional: cond_group, id.
#' #' @param n_bins      integer  Number of histogram bins (default: 10).
#' #' @param his_voltage numeric/NULL  Restrict to this voltage before plotting.
#' #' @param per_channel logical  TRUE uses per_chan_cond (averaged over all
#' #'                     voltages); FALSE uses per_voltage_cond at his_voltage.
#' #' @param plot_groups character/NULL  Restrict to these cond_group values.
#' #' @param cols        character/NULL  Fill colours per group.
#' #' @param plot_ylim   numeric/NULL  Y-axis limits.
#' #' @param plot_xlim   numeric/NULL  X-axis limits.
#' #' @param show_n      logical  TRUE annotates with pore count (n).
#' #' @param n_x         numeric/NULL  X position of n label (auto if NULL).
#' #' @param n_y         numeric/NULL  Y position of n label (auto if NULL).
#' #' @param txt_size    numeric  Base font size (default: 7).
#' #'
#' #' @return ggplot object.
#' hisPlot <- function(input_df,
#'                     n_bins      = 10,
#'                     his_voltage = NULL,
#'                     per_channel = TRUE,
#'                     plot_groups = NULL,
#'                     cols        = NULL,
#'                     plot_ylim   = NULL,
#'                     plot_xlim   = NULL,
#'                     show_n      = TRUE,
#'                     n_x         = NULL,
#'                     n_y         = NULL,
#'                     txt_size    = 7) {
#'   
#'   if (!is.null(plot_groups))
#'     input_df <- subset(input_df, input_df$cond_group %in% plot_groups)
#'   
#'   his_df <- if (!is.null(his_voltage)) {
#'     subset(input_df, input_df$voltage %in% his_voltage)
#'   } else {
#'     input_df
#'   }
#'   
#'   if (!isTRUE(per_channel)) {
#'     message("hisPlot: using per-voltage conductance — set his_voltage to avoid inflated n.")
#'     his_df$per_chan_cond <- his_df$per_voltage_cond
#'   }
#'   
#'   his_df$per_chan_cond <- abs(as.numeric(his_df$per_chan_cond))
#'   his_df               <- subset(his_df, his_df$per_chan_cond > 0.00001)
#'   
#'   fill_aes <- if ("cond_group" %in% colnames(input_df)) {
#'     ggplot2::aes(x = per_chan_cond, fill = cond_group)
#'   } else if ("id" %in% colnames(input_df)) {
#'     ggplot2::aes(x = per_chan_cond, fill = id)
#'   } else {
#'     ggplot2::aes(x = per_chan_cond)
#'   }
#'   
#'   if (!is.null(cols) && length(cols) < length(unique(input_df$cond_group)))
#'     cols <- rep(cols, length.out = length(unique(input_df$cond_group)))
#'   
#'   plt <- ggplot2::ggplot(his_df, fill_aes) +
#'     ggplot2::geom_histogram(bins = n_bins) +
#'     ggplot2::ylab("Counts") + ggplot2::xlab("Conductance / nS") +
#'     ggplot2::theme_bw(base_size = txt_size) +
#'     My_Theme(base_size = txt_size, txt_size = txt_size) +
#'     ggplot2::theme(legend.title = ggplot2::element_blank())
#'   
#'   if (!is.null(cols))      plt <- plt + ggplot2::scale_fill_manual(values = cols)
#'   if (!is.null(plot_ylim)) plt <- plt + ggplot2::scale_y_continuous(
#'     breaks = stefBreaks(plot_ylim)[[1]],
#'     labels = stefBreaks(plot_ylim)[[2]], limits = plot_ylim
#'   )
#'   if (!is.null(plot_xlim)) plt <- plt + ggplot2::scale_x_continuous(
#'     breaks = stefBreaks(plot_xlim)[[1]],
#'     labels = stefBreaks(plot_xlim)[[2]], limits = plot_xlim
#'   )
#'   
#'   if (isTRUE(show_n)) {
#'     n   <- length(unique(his_df$cm))
#'     n_x <- n_x %||% (max(his_df$per_chan_cond) * 0.75)
#'     n_y <- n_y %||% (n / 10)
#'     plt <- plt + ggplot2::annotate("text", x = n_x, y = n_y,
#'                                    label = paste0("n = ", n),
#'                                    size  = txt_size / ggplot2::.pt)
#'   }
#'   
#'   return(plt)
#' }
#' 
#' 
#' #' Fit probability distributions to a conductance histogram
#' #'
#' #' Adds fitted distribution curves to a hisPlot output and returns the
#' #' estimated peak conductance per group. Supports automatic model selection
#' #' by AIC, a single distribution applied to all groups, or different
#' #' distributions per group.
#' #'
#' #' @param his_plot     ggplot  Output of hisPlot (or any ggplot histogram).
#' #' @param distribution character/named character vector  Distribution(s) to fit:
#' #'   "auto", "normal", "gamma", "poisson", "exponential", "lognormal",
#' #'   "beta", "weibull", "logistic", "chisq", or "geom". Pass a named vector
#' #'   to assign different distributions per group (names must match group IDs).
#' #' @param line_size    numeric  Fitted curve line thickness (default: 0.5).
#' #' @param group_col    character  Name of the grouping column in the histogram
#' #'   data (default: "id").
#' #'
#' #' @return list:
#' #'   $plot   ggplot with fitted curves overlaid.
#' #'   $peaks  data.frame (id, distribution, peak, error).
#' hisFitr <- function(his_plot,
#'                     distribution = "auto",
#'                     line_size    = 0.5,
#'                     group_col    = "id") {
#'   
#'   available <- c("auto", "normal", "gamma", "poisson", "exponential",
#'                  "lognormal", "beta", "weibull", "logistic", "chisq", "geom")
#'   message("hisFitr: available distributions: ", paste(available, collapse = ", "))
#'   
#'   # ── AIC-based automatic selection ───────────────────────────────────────────
#'   select_best_dist <- function(x) {
#'     candidates <- setdiff(available, "auto")
#'     aics <- vapply(candidates, function(d) {
#'       tryCatch({
#'         f <- MASS::fitdistr(x, d)
#'         2 * length(f$estimate) - 2 * f$loglik
#'       }, error = function(e) Inf)
#'     }, numeric(1))
#'     
#'     if (all(is.infinite(aics)))
#'       stop("hisFitr: no distributions could be fitted to this group.")
#'     
#'     best <- candidates[which.min(aics)]
#'     message("  Auto-selected: ", best, " (AIC = ", round(min(aics), 2), ")")
#'     list(fit = MASS::fitdistr(x, best), dist = best)
#'   }
#'   
#'   # ── Extract histogram data ───────────────────────────────────────────────────
#'   gg_build  <- ggplot2::ggplot_build(his_plot)
#'   his_data  <- stats::na.omit(gg_build$plot$data)
#'   x_range   <- ggplot2::layer_scales(his_plot)$x$range$range
#'   hist_peak <- max(gg_build$data[[1]]$count)
#'   
#'   groups <- if (group_col %in% colnames(his_data)) {
#'     split(his_data, his_data[[group_col]])
#'   } else {
#'     list(all = his_data)
#'   }
#'   
#'   results <- vector("list", length(groups))
#'   
#'   for (i in seq_along(groups)) {
#'     gdata    <- groups[[i]]
#'     gid      <- names(groups)[i]
#'     
#'     # Determine which distribution to use for this group
#'     dist_req <- if (length(distribution) > 1) {
#'       if (!gid %in% names(distribution))
#'         stop("hisFitr: no distribution specified for group '", gid, "'.")
#'       as.character(distribution[gid])
#'     } else {
#'       as.character(distribution)
#'     }
#'     
#'     if (dist_req == "auto") {
#'       res          <- select_best_dist(gdata$per_chan_cond)
#'       fit          <- res$fit
#'       dist_used    <- res$dist
#'     } else {
#'       dist_used    <- dist_req
#'       fit          <- MASS::fitdistr(gdata$per_chan_cond, dist_used)
#'     }
#'     
#'     message("hisFitr: group '", gid, "' — fitting ", dist_used)
#'     
#'     # ── Density curve ────────────────────────────────────────────────────────
#'     x_vals <- seq(min(x_range), max(x_range), length.out = 1000)
#'     est    <- fit$estimate
#'     
#'     dens_fn <- switch(dist_used,
#'                       normal      = function(x) stats::dnorm(x, est["mean"], est["sd"]),
#'                       gamma       = function(x) stats::dgamma(x, est["shape"], est["rate"]),
#'                       poisson     = function(x) stats::dpois(x, est["lambda"]),
#'                       exponential = function(x) stats::dexp(x, est["rate"]),
#'                       lognormal   = function(x) stats::dlnorm(x, est["meanlog"], est["sdlog"]),
#'                       beta        = function(x) stats::dbeta(x, est["shape1"], est["shape2"]),
#'                       weibull     = function(x) stats::dweibull(x, est["shape"], est["scale"]),
#'                       logistic    = function(x) stats::dlogis(x, est["location"], est["scale"]),
#'                       chisq       = function(x) stats::dchisq(x, est["df"]),
#'                       geom        = function(x) stats::dgeom(x, est["prob"])
#'     )
#'     
#'     raw_dens   <- dens_fn(x_vals)
#'     scaled_dens <- raw_dens * (hist_peak / max(raw_dens, na.rm = TRUE))
#'     
#'     # ── Peak conductance (analytic where possible) ────────────────────────────
#'     peak_g <- switch(dist_used,
#'                      normal      = est["mean"],
#'                      gamma       = if (est["shape"] > 1) (est["shape"] - 1) / est["rate"] else NA,
#'                      exponential = 0,
#'                      lognormal   = exp(est["meanlog"] - est["sdlog"]^2),
#'                      weibull     = if (est["shape"] > 1)
#'                        est["scale"] * ((est["shape"] - 1) / est["shape"])^(1 / est["shape"])
#'                      else 0,
#'                      chisq       = if (est["df"] >= 2) est["df"] - 2 else 0,
#'                      logistic    = est["location"],
#'                      NA_real_
#'     )
#'     
#'     if (is.na(peak_g)) {
#'       opt    <- stats::optimize(function(x) -dens_fn(x),
#'                                 interval = c(min(x_range), max(x_range)))
#'       peak_g <- opt$minimum
#'     }
#'     
#'     curve_df <- data.frame(x = x_vals, y = scaled_dens, id = gid)
#'     curve_df$y[curve_df$y <= 0] <- NA
#'     
#'     his_plot <- his_plot +
#'       ggplot2::geom_line(
#'         data         = curve_df,
#'         ggplot2::aes(x = x, y = y, group = id),
#'         size         = line_size,
#'         colour       = "black",
#'         inherit.aes  = FALSE
#'       )
#'     
#'     results[[i]] <- data.frame(
#'       id           = gid,
#'       distribution = dist_used,
#'       peak         = peak_g,
#'       error        = min(fit$sd, na.rm = TRUE)
#'     )
#'   }
#'   
#'   peaks_df <- do.call(rbind, results)
#'   print(peaks_df)
#'   
#'   return(list(plot = his_plot, peaks = peaks_df))
#' }


# ── Section 11: Clustering ────────────────────────────────────────────────────

#' Cluster event data and return the dominant clusters
#'
#' Optionally splits the data by grouping variables before clustering. Returns
#' only rows belonging to the `top_n` largest clusters per group.
#'
#' @param input_df     data.frame  Input data.
#' @param split_by     character/NULL  Column(s) to split by before clustering.
#'                       NULL clusters the whole data frame.
#' @param cluster_vars character vector  Numeric columns to cluster on.
#' @param method       character  "kmeans" (default) or "hclust".
#' @param k            integer  Number of clusters (default: 3).
#' @param top_n        integer  Return rows from the `top_n` largest clusters
#'                       (default: 2).
#' @param seed         integer/NULL  RNG seed for reproducible kmeans.
#' @param ...          Forwarded to kmeans() or hclust().
#'
#' @return data.frame containing rows from the dominant clusters.
clustR <- function(input_df,
                   split_by     = NULL,
                   cluster_vars,
                   method       = "kmeans",
                   k            = 3,
                   top_n        = 2,
                   seed         = NULL,
                   ...) {
  
  input_df <- as.data.frame(input_df)
  
  cluster_one <- function(df) {
    sel <- df[, cluster_vars, drop = FALSE]
    
    if (!all(vapply(sel, is.numeric, logical(1)))) {
      warning("clustR: all cluster_vars must be numeric — skipping this group.")
      return(NULL)
    }
    if (nrow(sel) < k) {
      warning("clustR: fewer rows (", nrow(sel), ") than clusters (", k,
              ") — skipping.")
      return(NULL)
    }
    
    tryCatch({
      labels <- if (method == "kmeans") {
        if (!is.null(seed)) set.seed(seed)
        stats::kmeans(sel, centers = k, ...)$cluster
      } else if (method == "hclust") {
        stats::cutree(stats::hclust(stats::dist(sel), ...), k = k)
      } else {
        stop("clustR: method must be 'kmeans' or 'hclust'.")
      }
      
      df$cluster <- labels
      top_ids    <- as.integer(
        names(sort(table(labels), decreasing = TRUE))[seq_len(top_n)]
      )
      df[df$cluster %in% top_ids, , drop = FALSE]
      
    }, error = function(e) {
      warning("clustR: clustering failed — ", e$message)
      NULL
    })
  }
  
  if (!is.null(split_by) && length(split_by) > 0) {
    split_list <- split(input_df, lapply(split_by, function(v) input_df[[v]]),
                        drop = TRUE)
    rbindlist(Filter(Negate(is.null), lapply(split_list, cluster_one)),
              use.names = TRUE, fill = TRUE)
  } else {
    cluster_one(input_df)
  }
}


norm = function(col){
  #' 0-1 normalises any vector of daya
  #' @param col    (scalar)    any list of ints to be internally normalised
  #col = as.numeric(as.character(col))
  #col = (col - min(col))/(max(col)-min(col))
  col = (col - min(col, na.rm = TRUE)) / (max(col, na.rm = TRUE) - min(col, na.rm = TRUE))
  return(col)
}


#' Detect a single upward step (insertion event) in a normalised current trace.
#'
#' Strategy: find the window with the largest peak-to-peak range (most likely
#' to straddle the step), then estimate stable pre- and post-step levels using
#' modal averages in flanking windows.  Multiple guard-rails prevent the
#' indexing from going out of bounds even when the step is near either edge.
#'
#' @param y        Numeric vector. Normalised, absolute-valued current trace.
#' @param w        Integer. Width (in samples) of the rolling window used to
#'                 detect the peak-to-peak range. Larger values are more robust
#'                 to noise but less precise. Default 200.
#' @param win      Integer. Number of samples on each side of the jump used to
#'                 estimate pre- and post-step modal averages. Default 100.
#' @param ins_cut  Numeric (0-100). Only the first `ins_cut`% of the trace is
#'                 searched for the step, on the assumption that the insertion
#'                 event occurs early. Default 60.
#' @param edge_guard Numeric (0-1). Minimum fraction of `win` samples that must
#'                 be available for a pre/post estimate to be attempted; if
#'                 fewer samples exist the estimate falls back to the grand
#'                 median of whatever is available.  Default 0.25.
#'
#' @return A named list:
#'   \item{jump_index}{Sample index of the detected step in the *original* `y`
#'         (before truncation to `ins_cut`%).}
#'   \item{pre}{Modal-average current level before the step (pA or normalised
#'         units matching `y`).}
#'   \item{post}{Modal-average current level after the step.}
#'   \item{step_size}{Signed step amplitude (post - pre).}
#'   \item{n_pre}{Number of samples used for the pre estimate.}
#'   \item{n_post}{Number of samples used for the post estimate.}

detect_step <- function(y,
                        w          = 200,
                        win        = 100,
                        ins_cut    = 60,
                        edge_guard = 0.25) {
  
  # ── 1. Truncate search region ──────────────────────────────────────────────
  # Only search the first ins_cut% of the trace; insertions are expected early.
  search_end <- floor(length(y) * ins_cut / 100)
  # Need at least 2*w points to compute a meaningful rolling range
  if (search_end < 2 * w) {
    warning(sprintf(
      "detect_step: search region (%d pts) < 2*w (%d). Expanding to full trace.",
      search_end, 2 * w
    ))
    search_end <- length(y)
  }
  y_search <- y[seq_len(search_end)]
  
  # ── 2. Locate the step via maximum rolling peak-to-peak range ─────────────
  dy <- zoo::rollapply(y_search,
                       width = w,
                       FUN   = function(z) diff(range(z)),
                       fill  = NA)
  
  # which.max ignores NAs, so this is safe even with NA-filled edges
  j <- which.max(abs(dy))   # index in y_search (and in y, since y_search ⊆ y)
  
  # ── 3. Widen the jump neighbourhood to handle slow / gradual steps ─────────
  # Allow the transition zone to span ±1% (left) and ±10% (right) of the
  # search-region length, so a sluggish rising edge doesn't bias the post mean.
  tlen <- length(y_search)
  jmin <- max(1L,           floor(j - tlen * 0.01))
  jmax <- min(length(y),    ceiling(j + tlen * 0.10))  # use full y for post
  
  # ── 4. Define pre / post index ranges with boundary clamping ──────────────
  pre_lo  <- max(1L, jmin - win)
  pre_hi  <- max(1L, jmin - 1L)          # guard: never < 1
  post_lo <- min(length(y), jmax + 1L)
  post_hi <- min(length(y), jmax + win)  # guard: never > length(y)
  
  n_pre  <- pre_hi  - pre_lo  + 1L
  n_post <- post_hi - post_lo + 1L
  
  min_pts <- max(1L, floor(win * edge_guard))   # minimum acceptable window size
  
  # ── 5. Estimate pre-step level ─────────────────────────────────────────────
  if (pre_hi < pre_lo || n_pre < min_pts) {
    # Step is so early that there are almost no pre-step points; fall back to
    # the minimum of the whole search region (which should be near baseline).
    warning(sprintf(
      "detect_step: only %d pre-step sample(s) available (need >= %d). ",
      n_pre, min_pts
    ))
    pre <- modalAvg(y_search)
  } else {
    pre <- modalAvg(y[pre_lo:pre_hi])
  }
  
  # ── 6. Estimate post-step level ────────────────────────────────────────────
  if (post_lo > post_hi || n_post < min_pts) {
    # Step is so late in the search window that there is little post-step data;
    # use everything after the jump index to the end of the full trace.
    warning(sprintf(
      "detect_step: only %d post-step sample(s) available (need >= %d). ",
      n_post, min_pts
    ))
    fallback_lo <- min(length(y), jmax + 1L)
    if (fallback_lo > length(y)) {
      post <- pre   # degenerate: no data at all after step
    } else {
      post <- modalAvg(y[fallback_lo:length(y)])
    }
  } else {
    post <- modalAvg(y[post_lo:post_hi])
  }
  
  # ── 7. Return ──────────────────────────────────────────────────────────────
  list(
    jump_index = j,
    pre        = pre,
    post       = post,
    step_size  = post - pre,
    n_pre      = n_pre,
    n_post     = n_post
  )
}


# =============================================================================
# Orbit Mini electrophysiology helpers
# orbR      — read and summarise multi-channel Orbit CSV exports
# ivCut     — trim a multi-sweep IV dataframe to the first complete IV sweep
# =============================================================================
#' Read and summarise Orbit Mini multi-channel CSV exports.
#'
#' Recursively finds all CSV files under \code{dirnames}, reshapes each from
#' wide (one column per channel) to long format, subsets to the channels
#' named in the filename (e.g. "Ch1", "Ch3"), optionally prints a current
#' trace plot per channel, and returns either a data.frame of per-voltage
#' insertion metrics (output of \code{insertionAnalyser}) or the raw reshaped
#' long-format trace data for all files combined.
#'
#' If \code{csv_out_path} is provided, results are written incrementally as
#' each file completes. On startup the function checks which source files are
#' already recorded in that CSV and skips them, allowing interrupted runs to
#' be resumed without reprocessing. Pass \code{NULL} or \code{FALSE} to
#' disable all CSV writing and resume checking.
#'
#' @param dirnames     character    Path to folder containing Orbit CSV files.
#'                                  Searched recursively.
#' @param id           character    Experiment identifier appended to the
#'                                  returned data.frame's \code{id} column.
#' @param csv_out_path character/NULL/logical  Path to the output CSV. Created
#'                                  if absent; appended to if present. Files
#'                                  already recorded in the CSV are skipped on
#'                                  resume. Pass \code{NULL} or \code{FALSE}
#'                                  to disable CSV writing (default: NULL).
#' @param print_plt    logical      If TRUE, prints a \code{currPlot} trace for
#'                                  each channel in each file. Default TRUE.
#' @param output       character    "metrics" (default) returns the
#'                                  \code{insertionAnalyser} output per channel;
#'                                  "raw" returns the reshaped long-format trace
#'                                  data.frame before insertion analysis.
#'
#' @return data.frame of \code{insertionAnalyser} metrics (output = "metrics")
#'         or raw long-format trace data (output = "raw"), both with an added
#'         \code{id} column. When writing to CSV, already-processed files are
#'         excluded from the return value as their data is only in the CSV.
#'
#' @details
#' Expected CSV column naming (Orbit Mini default export):
#'   time..ms.                  — time in milliseconds
#'   voltage.channel..mV.       — applied voltage in mV
#'   current.channel.<N>..pA.   — current for channel N in pA (multiplied x1000
#'                                to convert from nA if needed — check units)
#'
#' Channel selection: the function extracts channel labels (e.g. "Ch1", "Ch3")
#' from the CSV filename. Only columns whose channel number appears in the
#' filename are kept. This means filenames must encode which channels were
#' active, e.g. "experiment_Ch1_Ch3_2024.csv".
#'
#' The source filename is written to a \code{source_file} column in the output
#' CSV. This column is used to detect already-processed files on resume.
#'
#' @seealso \code{insertionAnalyser}, \code{currPlot}
#'
#' @examples
#' \dontrun{
#' # No CSV output — return data only
#' iv_df  <- orbR("~/data/orbit_run1/", id = "sample_A")
#'
#' # Write to CSV; resumes automatically if interrupted
#' iv_df  <- orbR("~/data/orbit_run1/", id = "sample_A",
#'                csv_out_path = "~/results/sample_A.csv")
#' iv_df  <- orbR("~/data/orbit_run1/", id = "sample_A",
#'                csv_out_path = "~/results/sample_A.csv", print_plt = FALSE)
#' raw_df <- orbR("~/data/orbit_run1/", id = "sample_A",
#'                csv_out_path = "~/results/sample_A_raw.csv", output = "raw")
#' }
orbR <- function(dirnames,
                 id,
                 csv_out_path = NULL,
                 print_plt    = TRUE,
                 output       = "metrics") {
  
  if (output != "metrics" && output != "raw")
    stop("orbR: output must be 'metrics' or 'raw'.")
  
  # Normalise csv_out_path: FALSE or NULL both mean "no CSV"
  write_csv <- !is.null(csv_out_path) && !identical(csv_out_path, FALSE)
  
  # --- Locate CSV files ------------------------------------------------------
  csv_files <- list.files(
    path       = dirnames,
    pattern    = "\\.csv$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  if (length(csv_files) == 0L) {
    warning("orbR: no CSV files found under '", dirnames, "'")
    return(NULL)
  }
  
  message("orbR: found ", length(csv_files), " CSV file(s) under '", dirnames, "'")
  
  # --- Check output CSV for already-processed files --------------------------
  already_done <- character(0L)
  
  if (write_csv && file.exists(csv_out_path)) {
    existing <- tryCatch(
      utils::read.csv(csv_out_path, nrows = -1L),
      error = function(e) {
        warning("orbR: could not read existing CSV at '", csv_out_path,
                "' — processing all files. Error: ", e$message)
        NULL
      }
    )
    if (!is.null(existing) && "source_file" %in% names(existing)) {
      already_done <- unique(existing$source_file)
      message("orbR: ", length(already_done), " file(s) already in '",
              csv_out_path, "' — skipping.")
    } else if (!is.null(existing)) {
      warning("orbR: existing CSV has no 'source_file' column — ",
              "cannot determine which files are already processed; ",
              "processing all files.")
    }
  }
  
  # Filter to unprocessed files only
  message("Files found: ", length(csv_files))
  message("Already done: ", length(already_done))
  message("Pending: ",
          sum(!basename(csv_files) %in% already_done))
  print(
    basename(csv_files)[
      !basename(csv_files) %in% already_done
    ]
  )
  pending <- csv_files[!basename(csv_files) %in% already_done]
  
  if (length(pending) == 0L) {
    message("orbR: all files already processed")
    if (write_csv && file.exists(csv_out_path)) {
      return(
        data.table::fread(csv_out_path)
      )
    }
    return(NULL)
  }
  
  message("orbR: ", length(pending), " file(s) to process.")
  
  # --- Process each file -----------------------------------------------------
  all_raw <- list()
  
  df <- do.call(rbind, lapply(pending, function(cur_csv) {
    
    fname <- basename(cur_csv)
    message("orbR: reading ", fname)
    
    tmp <- tryCatch(
      read.csv(cur_csv),
      error = function(e) {
        message("orbR: failed to read ", fname, ": ", e$message)
        return(NULL)
      }
    )
    if (is.null(tmp) || nrow(tmp) == 0L) return(NULL)
    
    # --- Reshape wide -> long (one row per time point per channel) -----------
    channel_cols <- grep("^current\\.channel\\.", names(tmp))
    
    if (length(channel_cols) == 0L) {
      message("orbR: no current channel columns found in ", fname, " — skipping.")
      return(NULL)
    }
    
    long_df <- do.call(rbind, lapply(channel_cols, function(col_idx) {
      ch_num <- sub(".*\\.(\\d+)$", "\\1", names(tmp)[col_idx])
      data.frame(
        time        = tmp$time..ms.,
        voltage     = tmp$voltage.channel..mV.,
        current     = tmp[[col_idx]] * 1000,        # nA -> pA (verify units)
        cm          = paste0(fname, "_Ch", ch_num),
        id          = paste0("Ch", ch_num),
        source_file = fname
      )
    }))
    
    # --- Subset to channels named in the filename ----------------------------
    channels_in_name <- unique(
      regmatches(fname, gregexpr("Ch[0-9]+", fname))[[1]]
    )
    
    if (length(channels_in_name) == 0L) {
      message("orbR: no channel labels (e.g. Ch1) found in filename '",
              fname, "' — keeping all channels.")
      channels_in_name <- unique(long_df$id)
    }
    
    long_df <- long_df[long_df$id %in% channels_in_name, , drop = FALSE]
    
    if (nrow(long_df) == 0L) {
      message("orbR: no data remaining after channel filter for ", fname)
      return(NULL)
    }
    
    # Store raw long data before insertion analysis
    all_raw[[fname]] <<- long_df
    
    # --- Per-channel plot and insertion analysis ------------------------------
    file_df <- do.call(rbind, lapply(split(long_df, long_df$id), function(ch_df) {
      
      if (print_plt) {
        print(
          currPlot(ch_df, txt_size = 12) +
            ggplot2::ggtitle(paste(id, fname, unique(ch_df$id)))
        )
      }
      
      if (output == "raw") return(ch_df)
      
      ins <- tryCatch(
        insertionAnalyser(ch_df),
        error = function(e) {
          message("orbR: insertionAnalyser failed for ", unique(ch_df$id),
                  " in ", fname, ": ", e$message)
          NULL
        }
      )
      if (!is.null(ins)) ins$source_file <- fname
      ins
    }))
    
    # --- Write this file's results to CSV immediately ----------------------
    if (write_csv && !is.null(file_df) && nrow(file_df) > 0L) {
      data.table::fwrite(
        file_df,
        csv_out_path,
        append    = file.exists(csv_out_path),
        col.names = !file.exists(csv_out_path)
      )
      message("orbR: written ", nrow(file_df), " row(s) for ", fname)
    }
    
    file_df
  }))
  
  # --- Return raw trace data if requested ------------------------------------
  if (output == "raw") {
    raw_df <- do.call(rbind, all_raw)
    if (is.null(raw_df) || nrow(raw_df) == 0L) {
      warning("orbR: no raw data returned across all files.")
      return(NULL)
    }
    raw_df$id <- id
    return(raw_df)
  }
  
  if (is.null(df) || nrow(df) == 0L) {
    warning("orbR: no metrics data returned across all files.")
    return(NULL)
  }
  
  df$id <- id
  df
}
#' Trim a multi-sweep IV dataframe to the first complete IV sweep.
#'
#' Orbit Mini (and similar) recordings often contain multiple repeated IV
#' sweeps plus a leading zero-voltage equilibration period. This function:
#'   1. Estimates the typical voltage step length from the full trace.
#'   2. Keeps that many rows of the zero-voltage period immediately before
#'      the first non-zero voltage step (preserving the baseline).
#'   3. Truncates at the first return to zero voltage after the sweep starts,
#'      giving exactly one complete sweep with its preceding zero baseline.
#'
#' @param input_df  data.frame   Must contain a \code{voltage} column.
#'                               Typically a single-channel slice from
#'                               \code{orbR} output or similar.
#'
#' @return data.frame subset containing one baseline + one IV sweep.
#'
#' @details
#' Step length is estimated as the median gap between voltage transitions.
#' The zero-voltage window kept before the sweep is this median step length,
#' so the baseline duration matches the duration of each voltage step.
#'
#' If the trace never returns to zero after the first non-zero step, the
#' function returns from the estimated start to the end of the dataframe
#' with a warning.
#'
#' @examples
#' \dontrun{
#' # Trim to first sweep before passing to ivPlot
#' single_sweep <- ivCut(full_trace_df)
#' ivPlot(single_sweep)
#' }
ivCut <- function(input_df) {
  
  if (!"voltage" %in% names(input_df))
    stop("ivCut: 'voltage' column not found in input_df.")
  
  if (nrow(input_df) < 3L) {
    warning("ivCut: dataframe too short to detect voltage steps — returning as-is.")
    return(input_df)
  }
  
  # --- Detect voltage transition indices and estimate step length ------------
  dv       <- c(0, diff(input_df$voltage))
  step_idx <- which(dv != 0)
  
  if (length(step_idx) < 2L) {
    warning("ivCut: fewer than 2 voltage transitions found — returning as-is.")
    return(input_df)
  }
  
  step_lengths <- diff(step_idx)
  target_len   <- round(median(step_lengths))   # typical rows per voltage step
  
  # --- Find sweep boundaries -------------------------------------------------
  # Start: one step-length before the first non-zero voltage row
  first_nonzero <- which(input_df$voltage != 0)[1L]
  
  if (is.na(first_nonzero)) {
    warning("ivCut: no non-zero voltage rows found — returning as-is.")
    return(input_df)
  }
  
  new_start <- max(1L, first_nonzero - target_len)
  
  # End: the first return to zero voltage AFTER the sweep has started
  return_to_zero_idx <- which(
    seq_len(nrow(input_df)) > first_nonzero &
      input_df$voltage == 0
  )
  
  if (length(return_to_zero_idx) == 0L) {
    warning("ivCut: no return-to-zero found after sweep start — ",
            "returning from estimated start to end of dataframe.")
    return(input_df[new_start:nrow(input_df), , drop = FALSE])
  }
  
  sweep_end <- return_to_zero_idx[1L] - 1L
  
  message("ivCut: keeping rows ", new_start, ":", sweep_end,
          " (", sweep_end - new_start + 1L, " rows, ",
          "step length estimate: ", target_len, " rows)")
  
  input_df[new_start:sweep_end, , drop = FALSE]
}

#' Run event detection across Orbit Mini multi-channel exports.
#'
#' Iterates over a list of Orbit data directories, reads raw trace data via
#' \code{orbR}, applies baseline correction, conductance filtering, Bessel
#' filtering, open-pore current estimation, and event detection via
#' \code{eventDetect}. Already-processed datasets are skipped on resume by
#' checking the \code{id} column of the output CSV.
#'
#' @param input_data_files        list/character  Named list or character vector of
#'                                     directory paths to process. Each element
#'                                     is passed to \code{orbR} as \code{dirnames}.
#' @param csv_out_path character       Path to the output CSV. Created if absent;
#'                                     appended to if present. Pass \code{NULL}
#'                                     or \code{FALSE} to disable writing and
#'                                     resume checking (default: NULL).
#' @param msize        numeric         Conductance filter threshold in nS.
#'                                     Rows with conductance >= \code{msize} are
#'                                     dropped before filtering and event
#'                                     detection (default: 5).
#' @param sample_rate  numeric         Sample rate passed to \code{besselFilt}
#'                                     in Hz (default: 2000).
#' @param min_dwell    numeric         Minimum event duration in seconds passed
#'                                     to \code{eventDetect} (default: 0.001).
#' @param min_rows     numeric         Minimum number of rows required for a cm
#'                                     group to be processed, both before and
#'                                     after conductance filtering (default: 100).
#' @param print_plt    logical         Passed to \code{orbR}; if TRUE prints a
#'                                     \code{currPlot} trace per channel
#'                                     (default: FALSE).
#' @param max_gb      numeric          size in gb overwhich files will be excluded as 
#'                                     otherwise this will crash                                    
#'
#' @return Invisible NULL. Results are written incrementally to
#'         \code{csv_out_path} by \code{eventDetect}.
#'
#' @details
#' Resume behaviour: on startup, if \code{csv_out_path} exists and contains an
#' \code{id} column, any dataset whose \code{id} is already present is skipped
#' entirely. The \code{id} for each dataset is set to \code{basename(input_data_files[[i]])}.
#'
#' Processing pipeline per cm group:
#'   1. Baseline correction  — subtract median current at 0 mV.
#'   2. Conductance filter   — drop rows where |I|/|V| >= \code{msize}.
#'   3. Bessel filter        — \code{besselFilt} at \code{sample_rate}.
#'   4. OPC estimation       — \code{opcFinder} adds per-voltage open-pore current.
#'   5. Event detection      — \code{eventDetect} writes events to CSV.
#'
#' @seealso \code{orbR}, \code{eventDetect}, \code{besselFilt}, \code{opcFinder}
#'
#' @examples
#' \dontrun{
#' # Minimal call — process all dirs, write to CSV
#' orbEventAna(input_data_files, csv_out_path = "~/results/events.csv")
#'
#' # Adjust conductance threshold and minimum dwell time
#' orbEventAna(input_data_files,
#'             csv_out_path = "~/results/events.csv",
#'             msize        = 3,
#'             min_dwell    = 0.0005)
#'
#' # Dry run — no CSV output, results discarded (useful for testing)
#' orbEventAna(input_data_files, csv_out_path = NULL)
#' }
orbEventAna <- function(input_data_files,
                        csv_out_path = NULL,
                        msize        = 5,
                        sample_rate  = 2000,
                        min_dwell    = 0.001,
                        min_rows     = 100,
                        print_plt    = FALSE,
                        max_gb       = 5) {
  
  input_data_files = filter_large_datasets(input_data_files, max_gb)
  
  
  write_csv <- !is.null(csv_out_path) && !identical(csv_out_path, FALSE)
  
  # --- Load previously processed IDs ----------------------------------------
  done_ids <- character(0L)
  
  if (write_csv && file.exists(csv_out_path)) {
    existing <- tryCatch(
      data.table::fread(csv_out_path),
      error = function(e) {
        warning("orbEventAna: could not read '", csv_out_path,
                "' — starting fresh. Error: ", e$message)
        NULL
      }
    )
    if (!is.null(existing) && "id" %in% names(existing)) {
      done_ids <- unique(existing$id)
      message(sprintf("orbEventAna: %d previously processed ID(s) found — skipping.",
                      length(done_ids)))
    }
  } else {
    message("orbEventAna: no existing output file — starting fresh.")
  }
  
  # --- Main loop over directories --------------------------------------------
  for (i in seq_along(input_data_files)) {
    
    message(sprintf("\n[%d/%d] Processing: %s", i, length(input_data_files), basename(input_data_files[[i]])))
    
    # --- Read raw trace data -------------------------------------------------
    tmp <- tryCatch(
      orbR(input_data_files[[i]],
           id          = basename(input_data_files[[i]]),
           csv_out_path = NULL,
           print_plt   = print_plt,
           output      = "raw"),
      error = function(e) {
        message(sprintf("  orbR failed: %s", e$message))
        NULL
      }
    )
    
    if (is.null(tmp) || nrow(tmp) == 0L) {
      message("  No data returned — skipping.")
      next
    }
    
    data.table::setDT(tmp)
    cur_id <- tmp$id[1L]
    
    # --- Skip already-processed datasets -------------------------------------
    if (cur_id %in% done_ids) {
      message("  Already processed — skipping.")
      next
    }
    
    # --- Process each cm group -----------------------------------------------
    for (cm in unique(tmp$cm)) {
      
      message(sprintf("  cm = %s", cm))
      
      cmdf <- tmp[which(tmp$cm == cm)]
      
      if (nrow(cmdf) < min_rows) {
        message(sprintf("    Fewer than %d rows — skipping.", min_rows))
        next
      }
      
      # --- Baseline correction using 0 mV segment ---------------------------
      if (any(cmdf$voltage == 0)) {
        baseline <- stats::median(cmdf$current[cmdf$voltage == 0], na.rm = TRUE)
        cmdf[, current := current - baseline]
        message(sprintf("    Baseline correction: %.3f pA subtracted.", baseline))
      } else {
        message("    No 0 mV rows found — skipping baseline correction.")
      }
      
      # --- Conductance calculation ------------------------------------------
      # Avoid Inf/NaN at 0 mV by restricting to non-zero voltage rows only
      cmdf[, cond := NA_real_]
      nz <- cmdf$voltage != 0
      cmdf[nz, cond := abs(current) / abs(voltage)]
      
      cond_mean <- tryCatch(modalAvg(cmdf$cond, na.rm = TRUE), error = function(e) NA_real_)
      cond_sd   <- stats::sd(cmdf$cond, na.rm = TRUE)
      message(sprintf("    Conductance: %.3f +/- %.3f nS", cond_mean, cond_sd))
      
      # --- Conductance filter -----------------------------------------------
      # Applied before expensive Bessel filter to reduce data volume
      n_before <- nrow(cmdf)
      cmdf     <- cmdf[is.na(cond) | cond < msize]
      message(sprintf("    Conductance filter: retained %d/%d rows (threshold: %.1f nS).",
                      nrow(cmdf), n_before, msize))
      
      if (nrow(cmdf) < min_rows) {
        message(sprintf("    Fewer than %d rows after conductance filter — skipping.", min_rows))
        next
      }
      
      # --- Bessel filter -----------------------------------------------------
      t_filter <- Sys.time()
      cmdf <- tryCatch(
        besselFilt(cmdf, sample_rate = sample_rate),
        error = function(e) {
          message(sprintf("    besselFilt failed: %s", e$message))
          NULL
        }
      )
      
      gc(FALSE)
      
      if (is.null(cmdf)) next
      
      message(sprintf("    Bessel filter: %.2f s.",
                      as.numeric(difftime(Sys.time(), t_filter, units = "secs"))))
      
      # --- Open-pore current estimation -------------------------------------
      cmdf <- tryCatch(
        opcFinder(cmdf),
        error = function(e) {
          message(sprintf("    opcFinder failed: %s", e$message))
          NULL
        }
      )
      
      if (is.null(cmdf)) next
      
      # --- Event detection --------------------------------------------------
      tryCatch(
        eventDetect(cmdf,
                    csv_out_path = if (write_csv) csv_out_path else tempfile(),
                    min_dwell    = min_dwell),
        error = function(e) message(sprintf("    eventDetect failed: %s", e$message))
      )
    }
    ## clean up
    rm(cmdf)
    gc(FALSE)
    # --- Mark dataset complete -----------------------------------------------
    done_ids <- c(done_ids, cur_id)
    
    rm(tmp)
    gc(verbose = FALSE)
  }
  
  message("\norbEventAna: all files processed.")
  invisible(NULL)
}


filter_large_datasets <- function(paths, max_gb = 5) {
  
  keep <- logical(length(paths))
  
  for (i in seq_along(paths)) {
    
    sz_gb <- sum(
      file.info(
        list.files(paths[i],
                   recursive = TRUE,
                   full.names = TRUE)
      )$size,
      na.rm = TRUE
    ) / 1024^3
    
    if (sz_gb > max_gb) {
      message(sprintf(
        "Skipping %s (%.2f GB)",
        basename(paths[i]),
        sz_gb
      ))
      keep[i] <- FALSE
    } else {
      keep[i] <- TRUE
    }
  }
  
  paths[keep]
}

