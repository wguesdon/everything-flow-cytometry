# The shared parts of every cytokit recipe.
#
# cytokit is the command surface that an agent calls. It exists because every
# analysis in this repository is hard wired to one FlowRepository accession, and
# a scientist arriving with their own FCS files has no way in. A recipe here
# reads a path that the caller supplies and writes one bundle, so the same step
# works on any panel.
#
# The recipes themselves live in scripts/cytokit/. This file holds what they all
# need: the argument parser, the bundle, and the readers that describe a panel.
#
# tests/testthat.R sources every file in this folder into one environment, and
# the second definition of a name silently replaces the first. Check that a name
# is free before you add it.

# A recipe writes here unless the caller says otherwise. The path is relative to
# the working directory, which the CLI sets to the repository root.
kCytokitOutputRoot <- "output"

#' Read the `--key value` arguments of a recipe
#'
#' Every recipe takes its arguments in one form, so the parser is shared. An
#' argument the recipe does not name is an error rather than a value that is
#' dropped in silence, because a misspelled flag would otherwise change nothing
#' and report success.
#'
#' @param arguments A character vector, normally from `commandArgs(TRUE)`.
#' @param allowed The argument names the recipe accepts, without the `--`.
#' @param required The names that must be present.
#' @param flags The names that take no value and become `TRUE` when present.
#' @return A named list of character values, plus `TRUE` for each flag given.
#' @examples
#' ParseCytokitArguments(c("--data", "a.fcs"), allowed = "data")
#' @export
ParseCytokitArguments <- function(arguments, allowed, required = character(0),
                                  flags = character(0)) {
  parsed <- list()
  index <- 1
  while (index <= length(arguments)) {
    token <- arguments[index]
    if (!startsWith(token, "--")) {
      stop("Expected an argument beginning with '--', found '", token, "'.")
    }
    name <- substring(token, 3)
    if (!name %in% c(allowed, flags)) {
      stop("Unknown argument '--", name, "'. This recipe accepts: ",
           paste(sort(c(allowed, flags)), collapse = ", "), ".")
    }
    if (name %in% flags) {
      parsed[[name]] <- TRUE
      index <- index + 1
      next
    }
    if (index + 1 > length(arguments) ||
        startsWith(arguments[index + 1], "--")) {
      stop("The argument '--", name, "' needs a value.")
    }
    parsed[[name]] <- arguments[index + 1]
    index <- index + 2
  }
  missing <- setdiff(required, names(parsed))
  if (length(missing) > 0) {
    stop("These arguments are required: ",
         paste(paste0("--", missing), collapse = ", "), ".")
  }
  parsed
}

#' The timestamp that names a bundle
#'
#' @param at A `POSIXct` time. Defaults to now.
#' @return A string in `YYYY_MM_DD_HHMMSS` form.
#' @examples
#' CytokitTimestamp(as.POSIXct("2026-08-18 16:30:45", tz = "UTC"))
#' @export
CytokitTimestamp <- function(at = Sys.time()) {
  format(at, "%Y_%m_%d_%H%M%S")
}

#' Create the folder that one recipe writes into
#'
#' A recipe writes one folder and nothing outside it, so a run can be read,
#' copied or deleted whole.
#'
#' @param recipe The recipe name, for example `"inspect"`.
#' @param label A short name for the input, used in the folder name.
#' @param out_root The folder to create the bundle inside.
#' @param timestamp The timestamp to use. Defaults to now.
#' @return The bundle path.
#' @examples
#' \dontrun{
#' OpenCytokitBundle("inspect", "study")
#' }
#' @export
OpenCytokitBundle <- function(recipe, label, out_root = kCytokitOutputRoot,
                              timestamp = CytokitTimestamp()) {
  safe <- gsub("[^A-Za-z0-9]+", "_", label)
  safe <- gsub("^_+|_+$", "", safe)
  if (!nzchar(safe)) {
    safe <- "input"
  }
  bundle <- file.path(out_root, paste0(recipe, "_", safe, "_", timestamp))
  dir.create(bundle, recursive = TRUE, showWarnings = FALSE)
  bundle
}

#' Write a table into a bundle
#'
#' @param bundle The bundle path.
#' @param frame The `data.frame` to write.
#' @param name The file name, with the `.csv` extension.
#' @return The path, invisibly.
#' @examples
#' \dontrun{
#' WriteBundleTable(bundle, panel, "panel.csv")
#' }
#' @export
WriteBundleTable <- function(bundle, frame, name) {
  path <- file.path(bundle, name)
  utils::write.csv(frame, path, row.names = FALSE)
  invisible(path)
}

#' Record what produced a bundle
#'
#' The manifest is what makes a result answerable later. It carries the recipe,
#' every argument, the checksum of every input and the digest of the image, so a
#' number in a report can be traced to the run that produced it and to the files
#' that run read.
#'
#' @param bundle The bundle path.
#' @param recipe The recipe name.
#' @param arguments The parsed argument list.
#' @param inputs The input files that were read.
#' @param command The command line to repeat the run.
#' @return The bundle path, invisibly.
#' @examples
#' \dontrun{
#' CloseCytokitBundle(bundle, "inspect", arguments, files, command)
#' }
#' @export
CloseCytokitBundle <- function(bundle, recipe, arguments, inputs = character(0),
                               command = NA_character_) {
  # An input can be a folder, for example a saved hierarchy. md5sum fails on a
  # folder, so a folder is expanded into the files it holds.
  expanded <- unlist(lapply(inputs, function(path) {
    if (dir.exists(path)) {
      list.files(path, recursive = TRUE, full.names = TRUE)
    } else {
      path
    }
  }), use.names = FALSE)
  checksums <- if (length(expanded) > 0) {
    present <- expanded[file.exists(expanded) & !dir.exists(expanded)]
    stats::setNames(as.character(tools::md5sum(present)), basename(present))
  } else {
    list()
  }

  # An argument that holds a container path says nothing about which study was
  # run. Two studies would otherwise carry the same manifest.
  arguments <- lapply(arguments, function(value) {
    if (is.character(value) && length(value) == 1) DisplayPath(value) else value
  })

  record <- list(
    recipe = recipe,
    written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    arguments = arguments,
    inputs = as.list(checksums),
    image = Sys.getenv("CYTOKIT_IMAGE", unset = NA_character_),
    image_digest = Sys.getenv("CYTOKIT_IMAGE_DIGEST", unset = NA_character_),
    r_version = paste(R.version$major, R.version$minor, sep = ".")
  )

  path <- file.path(bundle, "manifest.json")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    writeLines(jsonlite::toJSON(record, auto_unbox = TRUE, pretty = TRUE,
                                null = "null"), path)
  } else {
    # A manifest a person can still read is better than none.
    writeLines(paste0(names(unlist(record)), ": ", unlist(record)), path)
  }

  writeLines(utils::capture.output(utils::sessionInfo()),
             file.path(bundle, "session_info.txt"))

  writeLines(c(
    paste("#", recipe),
    "",
    "This folder was written by one cytokit run. To produce it again:",
    "",
    "```bash",
    if (is.na(command)) paste("cytokit", recipe, "...") else command,
    "```",
    "",
    "`manifest.json` records every argument and the checksum of every input.",
    "`session_info.txt` records every package version."
  ), file.path(bundle, "REPRODUCE.md"))

  invisible(bundle)
}

#' Show a path as the caller typed it
#'
#' A recipe runs inside the container and sees `/indata` and `/outdata`. A
#' person
#' reading the output typed a host path, and a command printed with the
#' container
#' path in it cannot be run. The CLI records the host paths in the environment
#' so
#' that anything printed can be mapped back.
#'
#' Both variables carry the folder that the CLI mounted, and not the path the
#' caller named, so that one substitution is correct whether `--data` named a
#' folder or one file.
#'
#' @param path A path as the recipe sees it.
#' @return The path as the caller typed it, when the mapping is known.
#' @examples
#' DisplayPath("/indata")
#' @export
DisplayPath <- function(path) {
  if (is.null(path) || is.na(path)) {
    return(path)
  }
  data_host <- Sys.getenv("CYTOKIT_DATA_HOST", unset = "")
  out_host <- Sys.getenv("CYTOKIT_OUT_HOST", unset = "")
  controls_host <- Sys.getenv("CYTOKIT_CONTROLS_HOST", unset = "")
  # A trailing slash on the host folder would double the separator, so it is
  # dropped before the two are joined.
  if (nzchar(data_host)) {
    path <- sub("^/indata/?", paste0(sub("/$", "", data_host), "/"), path)
  }
  if (nzchar(out_host)) {
    path <- sub("^/outdata/?", paste0(sub("/$", "", out_host), "/"), path)
  }
  if (nzchar(controls_host)) {
    path <- sub("^/incontrols/?",
                paste0(sub("/$", "", controls_host), "/"), path)
  }
  # A table the caller named mounts at /in_<flag>. The CLI records the folder
  # of each one, so a path printed from any of them maps back.
  for (flag in c("template", "definitions", "metadata", "counts",
                 "proportions", "gates", "clusters", "claims", "results")) {
    host <- Sys.getenv(paste0("CYTOKIT_", toupper(flag), "_HOST"), unset = "")
    if (nzchar(host)) {
      path <- sub(paste0("^/in_", flag, "/?"),
                  paste0(sub("/$", "", host), "/"), path)
    }
  }
  sub("(.)/$", "\\1", path)
}

#' List the FCS files at a path
#'
#' A recipe takes one file or a folder of them, because a scientist has either.
#'
#' @param path A file or a folder.
#' @param recursive Whether to descend into sub folders.
#' @return A character vector of file paths, sorted.
#' @examples
#' \dontrun{
#' FcsFilesIn("study/")
#' }
#' @export
FcsFilesIn <- function(path, recursive = FALSE) {
  if (!file.exists(path)) {
    stop("The path does not exist: ", DisplayPath(path))
  }
  if (!dir.exists(path)) {
    return(path)
  }
  files <- list.files(path, pattern = "\\.fcs$", ignore.case = TRUE,
                      full.names = TRUE, recursive = recursive)
  if (length(files) == 0) {
    stop("No FCS file was found in ", DisplayPath(path),
         if (recursive) "." else ". Add --recursive to look in sub folders.")
  }
  sort(files)
}

#' Describe one FCS file without reading its events
#'
#' The header carries the event count, the parameter count, the instrument and
#' the compensation keywords. Reading only the header keeps this fast on a file
#' of millions of events.
#'
#' @param path The FCS file.
#' @return A one row `data.frame`.
#' @examples
#' \dontrun{
#' DescribeFcsFile("study/sample.fcs")
#' }
#' @export
DescribeFcsFile <- function(path) {
  frame <- flowCore::read.FCS(path, which.lines = 1, truncate_max_range = FALSE,
                              transformation = FALSE)
  keywords <- flowCore::keyword(frame)
  Keyword <- function(name) {
    value <- keywords[[name]]
    if (is.null(value) || length(value) != 1) NA_character_ else
      as.character(value)
  }
  state <- ReadCompensationState(frame)

  data.frame(
    file = basename(path),
    events = suppressWarnings(as.integer(Keyword("$TOT"))),
    parameters = suppressWarnings(as.integer(Keyword("$PAR"))),
    cytometer = Keyword("$CYT"),
    acquired = Keyword("$DATE"),
    compensation_state = state$state,
    apply_compensation = state$apply_keyword,
    matrix_size = state$matrix_size,
    stringsAsFactors = FALSE
  )
}

#' Describe the panel of one FCS file, one row per detector
#'
#' The marker names are what an agent needs to draft a gating template or a cell
#' type definitions table, so this reports them beside the detector and says
#' which detectors are scatter and which carry a stain.
#'
#' @param path The FCS file.
#' @return A `data.frame` with one row per parameter.
#' @examples
#' \dontrun{
#' DescribeFcsPanel("study/sample.fcs")
#' }
#' @export
DescribeFcsPanel <- function(path) {
  frame <- flowCore::read.FCS(path, which.lines = 1, truncate_max_range = FALSE,
                              transformation = FALSE)
  channels <- DescribeChannels(frame)
  parameters <- flowCore::pData(flowCore::parameters(frame))

  is_scatter <- grepl("^(FSC|SSC)", channels$channel, ignore.case = TRUE)
  is_time <- grepl("^time$", channels$channel, ignore.case = TRUE)

  data.frame(
    channel = channels$channel,
    marker = channels$marker,
    kind = ifelse(is_scatter, "scatter",
                  ifelse(is_time, "time",
                         ifelse(channels$is_marker, "stain", "unnamed"))),
    range = suppressWarnings(as.numeric(parameters$range)),
    minimum = suppressWarnings(as.numeric(parameters$minRange)),
    maximum = suppressWarnings(as.numeric(parameters$maxRange)),
    stringsAsFactors = FALSE
  )
}

#' The marker names of a panel, in the form a table column needs
#'
#' A definitions table has one column per marker, and a gating template names a
#' marker in its `dims` column. Both need a name for each fluorescence detector.
#'
#' Many deposits leave the `$PnS` keyword empty, so the file names no markers at
#' all. The detector name then stands in for the marker, because a table has to
#' have a column name and `APC-A` is at least true. It is not the antibody, so a
#' scientist has to supply that mapping before a label means anything.
#'
#' @param panel The output of [DescribeFcsPanel()].
#' @param fallback Whether to use the detector name when a detector carries no
#'   marker. `FALSE` returns only the detectors the file actually named.
#' @return A character vector of names, without the scatter and time detectors.
#' @examples
#' \dontrun{
#' PanelMarkers(DescribeFcsPanel("study/sample.fcs"))
#' }
#' @export
PanelMarkers <- function(panel, fallback = TRUE) {
  named <- panel[panel$kind == "stain", , drop = FALSE]
  names <- named$marker[nzchar(named$marker)]
  if (!fallback) {
    return(unique(names))
  }
  unnamed <- panel$channel[panel$kind == "unnamed"]
  unique(c(names, unnamed))
}

#' Report whether a panel names its markers
#'
#' A file with no `$PnS` keyword forces every later step to work from detector
#' names. That has to be said once, plainly, rather than discovered when a cell
#' type label turns out to mean nothing.
#'
#' @param panel The output of [DescribeFcsPanel()].
#' @return A one row `data.frame` with the counts and a `state`.
#' @examples
#' \dontrun{
#' PanelNamingState(DescribeFcsPanel("study/sample.fcs"))
#' }
#' @export
PanelNamingState <- function(panel) {
  named <- sum(panel$kind == "stain")
  unnamed <- sum(panel$kind == "unnamed")
  state <- if (named == 0 && unnamed > 0) {
    "no marker names"
  } else if (unnamed > 0) {
    "some detectors unnamed"
  } else {
    "every detector named"
  }
  data.frame(named = named, unnamed = unnamed, state = state,
             stringsAsFactors = FALSE)
}

#' Decide whether a panel came from a mass cytometer or a fluorescence one
#'
#' The two need different work. A mass cytometry file carries no scatter and no
#' spillover, so a scatter gate and a compensation step are both wrong on it.
#' FR-FCM-Z244 in this archive is a CyTOF deposit with 66 detectors and no
#' scatter channel, and a recipe that assumes FSC and SSC produces nothing on
#' it.
#'
#' The check reads the instrument keyword and the channel names, because a
#' deposit carries one, the other, or both.
#'
#' @param panel A panel table from `DescribeFcsPanel`.
#' @param cytometer The `$CYT` keyword, or `NA` when the file carries none.
#' @return A list with `kind`, either `"mass"` or `"fluorescence"`, `reason`,
#'   and `has_scatter`.
#' @examples
#' panel <- data.frame(channel = c("Y89Di", "Pd102Di", "Event_length"),
#'                     kind = c("stain", "stain", "unnamed"))
#' AcquisitionKind(panel, "DVSSCIENCES-CYTOF-6.7.1014")$kind
#' @export
AcquisitionKind <- function(panel, cytometer = NA_character_) {
  has_scatter <- any(panel$kind == "scatter")
  instrument <- if (length(cytometer) == 0 || all(is.na(cytometer))) {
    ""
  } else {
    paste(unique(stats::na.omit(cytometer)), collapse = " ")
  }
  # A mass cytometer names a detector after the isotope it counts, for example
  # Y89Di or Pd102Di. Three of them is past coincidence.
  mass_channels <- grep("^[A-Z][a-z]?[0-9]{2,3}Di$", panel$channel,
                        value = TRUE)
  by_instrument <- grepl("CYTOF|HELIOS|DVSSCIENCES", instrument,
                         ignore.case = TRUE)
  by_channels <- length(mass_channels) >= 3

  if (by_instrument || by_channels) {
    reason <- if (by_instrument) {
      paste0("the instrument keyword reads ", instrument)
    } else {
      paste0(length(mass_channels), " detectors are named after an isotope")
    }
    return(list(kind = "mass", reason = reason, has_scatter = has_scatter))
  }
  list(kind = "fluorescence",
       reason = if (nzchar(instrument)) {
         paste0("the instrument keyword reads ", instrument)
       } else {
         "no instrument keyword, and no detector named after an isotope"
       },
       has_scatter = has_scatter)
}

#' Run an expression and collect what it reported instead of printing each line
#'
#' flowCore reports a malformed header once per read. On a deposit of 28 files
#' that is 112 copies of one line, and the report scrolls out of view. The
#' report still matters, so a recipe counts it and states it once. Hiding it
#' would leave a scientist reading a panel from a file that did not parse
#' cleanly.
#'
#' The reader writes some of these with `cat` rather than with `warning`, so
#' printed output is captured as well as the conditions.
#'
#' @param expression The expression to run.
#' @return A list with `value`, and `notes`, a table of the distinct lines with
#'   the number of times each was reported.
#' @examples
#' CollectNotes({message("late"); 1})$notes
#' @export
CollectNotes <- function(expression) {
  raised <- character(0)
  value <- NULL
  printed <- utils::capture.output(
    value <- withCallingHandlers(
      expression,
      warning = function(condition) {
        raised <<- c(raised, trimws(conditionMessage(condition)))
        invokeRestart("muffleWarning")
      },
      message = function(condition) {
        raised <<- c(raised, trimws(conditionMessage(condition)))
        invokeRestart("muffleMessage")
      }
    )
  )
  raised <- c(raised, trimws(printed))
  raised <- raised[nzchar(raised)]
  notes <- if (length(raised) == 0) {
    data.frame(note = character(0), times = integer(0),
               stringsAsFactors = FALSE)
  } else {
    counted <- table(raised)
    data.frame(note = names(counted), times = as.integer(counted),
               stringsAsFactors = FALSE, row.names = NULL)
  }
  list(value = value, notes = notes[order(-notes$times), , drop = FALSE])
}

# The keywords that say which specimen a file holds, and which run it came
# from. An FCS file carries no metadata table, so this is the only machine
# readable record of the design that ships with the data.
kIdentityKeywords <- c(
  original_file = "$FIL",
  specimen = "$SRC",
  tube = "TUBE NAME",
  experiment = "EXPERIMENT NAME",
  specimen_number = "$SMNO",
  project = "$PROJ",
  well = "$WELLID",
  started = "$BTIM",
  ended = "$ETIM"
)

#' Read the keywords that identify the specimen and the run
#'
#' A scientist arrives with a folder of FCS files and no metadata table, and the
#' treatment of each file has to come from somewhere. The file name is a guess.
#' These keywords are what the instrument recorded, so they are evidence.
#'
#' `$SRC` splits a folder by donor, and `TUBE NAME` gives the position in the
#' run. On FR-FCM-ZZCA the tube numbers read 1, 2, 13, 15 and 16, which says
#' that eleven tubes of that run are absent from the deposit, and `$SRC` splits
#' the five files into two groups. Neither fact is in a file name.
#'
#' @param path The FCS file.
#' @return A one row `data.frame` with one column per keyword that the file
#'   carries. A keyword the file leaves out is `NA`.
#' @examples
#' \dontrun{
#' DescribeFcsIdentity("study/sample.fcs")
#' }
#' @export
DescribeFcsIdentity <- function(path) {
  frame <- flowCore::read.FCS(path, which.lines = 1, truncate_max_range = FALSE,
                              transformation = FALSE)
  keywords <- flowCore::keyword(frame)
  values <- vapply(kIdentityKeywords, function(name) {
    value <- keywords[[name]]
    if (is.null(value) || length(value) != 1) NA_character_ else
      trimws(as.character(value))
  }, character(1))
  values[!nzchar(values) & !is.na(values)] <- NA_character_

  identity <- as.data.frame(as.list(values), stringsAsFactors = FALSE)
  identity <- cbind(file = basename(path), identity, stringsAsFactors = FALSE)
  rownames(identity) <- NULL
  identity
}

#' Report what a set of identity keywords says about the design
#'
#' One value across every file tells a scientist about the run. A value per file
#' is an identifier. A keyword with two values and five files splits the folder,
#' and that split is the only candidate for a grouping that ships with the data.
#' The three cases need different reading, so each row carries its role.
#'
#' A timestamp never groups anything, so it is called out rather than counted as
#' a split.
#'
#' @param identity The table from `DescribeFcsIdentity`, one row per file.
#' @return A `data.frame` with one row per keyword that carries a value, its
#'   role, its number of distinct values, and those values.
#' @examples
#' identity <- data.frame(file = c("a.fcs", "b.fcs", "c.fcs"),
#'                        specimen = c("PBMC", "PBMC", "PBMC_001"),
#'                        tube = c("1", "2", "3"))
#' IdentitySplits(identity)
#' @export
IdentitySplits <- function(identity) {
  timing <- c("started", "ended")
  columns <- setdiff(names(identity), "file")
  files <- nrow(identity)
  rows <- lapply(columns, function(column) {
    values <- identity[[column]]
    present <- values[!is.na(values)]
    if (length(present) == 0) {
      return(NULL)
    }
    distinct <- unique(present)
    role <- if (column %in% timing) {
      "timing"
    } else if (length(distinct) == 1) {
      "constant"
    } else if (length(distinct) == files && files > 1) {
      "identifier"
    } else {
      "grouping"
    }
    data.frame(
      keyword = column,
      role = role,
      # A grouping is what a scientist confirms, so every value is listed. An
      # identifier has one value per file, and eight of them are unreadable.
      files_with_it = length(present),
      distinct_values = length(distinct),
      values = paste(utils::head(distinct, if (role == "grouping") 8 else 3),
                     collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame(keyword = character(0), role = character(0),
                      files_with_it = integer(0), distinct_values = integer(0),
                      values = character(0), stringsAsFactors = FALSE))
  }
  result <- do.call(rbind, rows)
  # A grouping is what a scientist has to confirm, so it reads first.
  order_of_role <- match(result$role,
                         c("grouping", "identifier", "constant", "timing"))
  result[order(order_of_role, result$keyword), , drop = FALSE]
}

#' Say where each marker name came from
#'
#' `PanelMarkers` falls back to the detector name when a file leaves `$PnS`
#' empty. The two kinds of name read the same in a table and mean different
#' things, because a detector name is not an antibody. A recipe that writes a
#' marker into a file has to be able to say which is which.
#'
#' @param panel The output of [DescribeFcsPanel()].
#' @return A `data.frame` with `marker`, `channel` and `source`, where `source`
#'   is `"antibody"` when the file named it and `"detector"` when it did not.
#' @examples
#' panel <- data.frame(channel = c("B515-A", "APC-A"), marker = c("CD3", ""),
#'                     kind = c("stain", "unnamed"), stringsAsFactors = FALSE)
#' PanelMarkerSource(panel)
#' @export
PanelMarkerSource <- function(panel) {
  named <- panel[panel$kind == "stain" & nzchar(panel$marker), , drop = FALSE]
  unnamed <- panel[panel$kind == "unnamed", , drop = FALSE]
  result <- rbind(
    data.frame(marker = named$marker, channel = named$channel,
               source = rep("antibody", nrow(named)), stringsAsFactors = FALSE),
    data.frame(marker = unnamed$channel, channel = unnamed$channel,
               source = rep("detector", nrow(unnamed)),
               stringsAsFactors = FALSE)
  )
  result[!duplicated(result$marker), , drop = FALSE]
}

#' Read a `--markers` argument, which either selects or maps
#'
#' A scientist whose file names its markers picks a few of them. A scientist
#' whose file leaves `$PnS` empty has nothing to pick from, and the only useful
#' thing they can give is the mapping from the detector to the antibody. The
#' same argument covers both, because a pair carries an equals sign.
#'
#' @param argument The raw `--markers` string, or `NULL`.
#' @param available The marker names the panel offers, from `PanelMarkers`.
#' @param channels Every detector name in the panel.
#' @return A `data.frame` with `column`, the name to write, and `from`, the
#'   panel name it came from. `NULL` when `argument` is `NULL`.
#' @examples
#' ParseMarkerArgument("APC-A=CD3,FITC-A=CD4", c("APC-A", "FITC-A"),
#'                     c("FSC-A", "APC-A", "FITC-A"))
#' @export
ParseMarkerArgument <- function(argument, available, channels) {
  if (is.null(argument)) {
    return(NULL)
  }
  asked <- trimws(strsplit(argument, ",")[[1]])
  asked <- asked[nzchar(asked)]
  if (length(asked) == 0) {
    stop("--markers is empty. Give a marker name, or a detector=antibody pair.")
  }
  pairs <- grepl("=", asked, fixed = TRUE)

  mapped <- lapply(asked[pairs], function(item) {
    parts <- trimws(strsplit(item, "=", fixed = TRUE)[[1]])
    if (length(parts) != 2 || !all(nzchar(parts))) {
      stop("A mapping reads detector=antibody. This one does not: ", item)
    }
    if (!parts[1] %in% c(channels, available)) {
      stop("This detector is not in the panel: ", parts[1],
           "\nThe panel carries: ", paste(channels, collapse = ", "))
    }
    data.frame(column = parts[2], from = parts[1], stringsAsFactors = FALSE)
  })

  plain <- asked[!pairs]
  unknown <- setdiff(plain, available)
  if (length(unknown) > 0) {
    stop("These markers are not in the panel: ",
         paste(unknown, collapse = ", "),
         "\nThe panel carries: ", paste(available, collapse = ", "),
         "\nWhen the panel names no antibody, give the mapping instead, ",
         "for example --markers \"", channels[1], "=CD3\".")
  }
  selected <- if (length(plain) > 0) {
    data.frame(column = plain, from = plain, stringsAsFactors = FALSE)
  } else {
    NULL
  }
  result <- do.call(rbind, c(mapped, list(selected)))
  if (anyDuplicated(result$column) > 0) {
    stop("A column is named twice: ",
         paste(unique(result$column[duplicated(result$column)]),
               collapse = ", "))
  }
  result
}

#' Print what the reader reported, capped, and write the full list
#'
#' A 27 colour panel produces 27 distinct progress lines from the spillover
#' matcher, and printing all of them hides the result they sit above. The count
#' and the first few lines are what a reader needs, and the rest belong in a
#' file.
#'
#' @param notes A table with `note` and `times`, from [CollectNotes()].
#' @param bundle The bundle folder, or `NULL` to print without writing.
#' @param limit How many distinct lines to print. Defaults to 5.
#' @return The aggregated table, invisibly.
#' @examples
#' ReportNotes(data.frame(note = "uneven tokens", times = 3L), NULL)
#' @export
ReportNotes <- function(notes, bundle, limit = 5) {
  if (is.null(notes) || nrow(notes) == 0) {
    return(invisible(notes))
  }
  notes <- stats::aggregate(times ~ note, data = notes, FUN = sum)
  notes <- notes[order(-notes$times, notes$note), , drop = FALSE]
  cat("\nThe FCS reader raised ", sum(notes$times), " note(s) over ",
      nrow(notes), " distinct line(s):\n", sep = "")
  shown <- utils::head(notes, limit)
  for (index in seq_len(nrow(shown))) {
    cat("  ", shown$times[index], "x  ", shown$note[index], "\n", sep = "")
  }
  if (nrow(notes) > limit) {
    cat("  and ", nrow(notes) - limit, " more, all in reader_notes.csv\n",
        sep = "")
  }
  if (!is.null(bundle)) {
    WriteBundleTable(bundle, notes, "reader_notes.csv")
  }
  invisible(notes)
}

# Every recipe that draws a random subset uses this seed unless the caller
# names another. robustbase::covMcd draws one, and so does
# flowStats::norm2Filter, which spillover_ng calls while it gates each control.
kCytokitSeed <- 42

#' Set the seed a recipe runs under, and report it
#'
#' Without a seed a spillover matrix moves between runs, the compensated values
#' move with it, and a count in a report cannot be produced again. Two unseeded
#' runs of the OMIP-58 analysis differed by six live lymphocytes of 598880 and
#' swung a fitted cut from 60.97 to 79.76 percent of T cells.
#'
#' @param arguments The parsed arguments, which may carry `seed`.
#' @return The seed that was set, invisibly.
#' @examples
#' SetCytokitSeed(list(seed = "7"))
#' @export
SetCytokitSeed <- function(arguments = list()) {
  seed <- if (is.null(arguments$seed)) {
    kCytokitSeed
  } else {
    value <- suppressWarnings(as.integer(arguments$seed))
    if (is.na(value)) {
      stop("--seed must be a whole number, not ", arguments$seed, ".")
    }
    value
  }
  set.seed(seed)
  invisible(seed)
}

#' Shorten a bundle name so that a chained recipe does not stack timestamps
#'
#' A bundle is named `<recipe>_<label>_<timestamp>`. When one recipe reads
#' another's bundle and uses its name as a label, the new bundle carries two
#' recipe names and two timestamps. This strips the recipe name and the
#' timestamp so that only the study is left.
#'
#' @param name A bundle folder name.
#' @return The label inside it.
#' @examples
#' ShortLabel("gate_my_study_2026_08_18_205932")
#' @export
ShortLabel <- function(name) {
  name <- basename(name)
  name <- sub("^(inspect|compensate|gate|cluster|annotate|proportions|compare|claims)_",
              "", name)
  # The timestamp can be the whole of what is left, so the leading underscore
  # is optional.
  name <- sub("(^|_)[0-9]{4}_[0-9]{2}_[0-9]{2}_[0-9]{6}$", "", name)
  if (nzchar(name)) name else "study"
}
