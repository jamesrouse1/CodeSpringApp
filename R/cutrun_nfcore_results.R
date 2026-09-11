
# -------------------------------------------------------------------------
# nf-core/cutandrun Results Explorer adapters
# -------------------------------------------------------------------------
#
# Keep the existing native CUT&RUN Results Explorer and downstream tools.
# These helpers only translate CodeSpring sample names to nf-core output
# names/locations when the selected project uses nf-core/cutandrun.


cutrun_nfcore_results_paths <- function(project) {
  paths <- cutrun_nfcore_submission_paths(project)

  list(
    run_dir = paths$run_dir,
    output_dir = paths$output_dir,
    mapping_path = paths$mapping_path,

    raw_fastqc = file.path(
      paths$output_dir,
      "01_prealign",
      "pretrim_fastqc"
    ),

    trimmed_fastqc = file.path(
      paths$output_dir,
      "01_prealign",
      "trimgalore",
      "fastqc"
    ),

    alignment = file.path(
      paths$output_dir,
      "02_alignment",
      "bowtie2"
    ),

    bedgraph = file.path(
      paths$output_dir,
      "03_peak_calling",
      "01_bam_to_bedgraph"
    ),

    bigwig = file.path(
      paths$output_dir,
      "03_peak_calling",
      "03_bed_to_bigwig"
    ),

    called_peaks = file.path(
      paths$output_dir,
      "03_peak_calling",
      "04_called_peaks"
    ),

    multiqc = file.path(
      paths$output_dir,
      "04_reporting",
      "multiqc"
    )
  )
}



.cutrun_nfcore_mapping_cache <- new.env(
  parent = emptyenv()
)

cutrun_nfcore_result_mapping <- function(
  project,
  cache_seconds = 15
) {
  path <- cutrun_nfcore_results_paths(
    project
  )$mapping_path

  empty <- data.frame(
    stringsAsFactors = FALSE
  )

  if (
    !file.exists(path) ||
    dir.exists(path)
  ) {
    return(empty)
  }

  key <- normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )

  now <- as.numeric(
    Sys.time()
  )

  if (
    exists(
      key,
      envir = .cutrun_nfcore_mapping_cache,
      inherits = FALSE
    )
  ) {
    cached <- get(
      key,
      envir = .cutrun_nfcore_mapping_cache,
      inherits = FALSE
    )

    age <- now - cached$time

    if (
      is.finite(age) &&
      age >= 0 &&
      age < cache_seconds
    ) {
      return(
        cached$data
      )
    }
  }

  x <- tryCatch(
    utils::read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) empty
  )

  required <- c(
    "codespring_sample",
    "nfcore_group",
    "replicate"
  )

  if (
    !NROW(x) ||
    !all(
      required %in% names(x)
    )
  ) {
    x <- empty
  }

  assign(
    key,
    list(
      time = now,
      data = x
    ),
    envir = .cutrun_nfcore_mapping_cache
  )

  x
}

cutrun_nfcore_sample_meta <- function(
  project,
  sample
) {
  mapping <- cutrun_nfcore_result_mapping(
    project
  )

  if (!NROW(mapping)) {
    return(NULL)
  }

  sample <- trimws(
    as.character(
      sample %||% ""
    )
  )

  hit <- mapping[
    trimws(
      as.character(
        mapping$codespring_sample
      )
    ) == sample,
    ,
    drop = FALSE
  ]

  if (NROW(hit) != 1L) {
    return(NULL)
  }

  hit[1, , drop = FALSE]
}


cutrun_nfcore_result_key <- function(x) {
  gsub(
    "[^a-z0-9]+",
    "",
    tolower(
      as.character(
        x %||% ""
      )
    )
  )
}



.cutrun_nfcore_result_file_cache <- new.env(
  parent = emptyenv()
)

cutrun_nfcore_result_files <- function(
  directory,
  pattern = NULL,
  cache_seconds = 15
) {
  if (
    !dir.exists(directory)
  ) {
    return(character(0))
  }

  directory <- normalizePath(
    directory,
    winslash = "/",
    mustWork = FALSE
  )

  pattern_key <- if (
    is.null(pattern)
  ) {
    ""
  } else {
    as.character(pattern)[[1]]
  }

  key <- paste(
    directory,
    pattern_key,
    sep = "\r"
  )

  now <- as.numeric(
    Sys.time()
  )

  if (
    exists(
      key,
      envir = .cutrun_nfcore_result_file_cache,
      inherits = FALSE
    )
  ) {
    cached <- get(
      key,
      envir = .cutrun_nfcore_result_file_cache,
      inherits = FALSE
    )

    age <- now - cached$time

    if (
      is.finite(age) &&
      age >= 0 &&
      age < cache_seconds
    ) {
      return(
        cached$files
      )
    }
  }

  files <- tryCatch(
    list.files(
      directory,
      pattern = pattern,
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    ),
    error = function(e) character(0)
  )

  if (length(files)) {
    exists_file <- file.exists(
      files
    )

    is_directory <- rep(
      FALSE,
      length(files)
    )

    if (any(exists_file)) {
      is_directory[exists_file] <-
        dir.exists(
          files[exists_file]
        )
    }

    files <- files[
      exists_file &
        !is_directory
    ]
  }

  assign(
    key,
    list(
      time = now,
      files = files
    ),
    envir = .cutrun_nfcore_result_file_cache
  )

  files
}

cutrun_nfcore_match_sample_file <- function(
  project,
  sample,
  files,
  read = "",
  prefer = character(0)
) {
  files <- as.character(files %||% character(0))

  files <- files[
    nzchar(files) &
      file.exists(files) &
      !dir.exists(files)
  ]

  if (!length(files)) {
    return("")
  }

  meta <- cutrun_nfcore_sample_meta(
    project,
    sample
  )

  if (!NROW(meta)) {
    return("")
  }

  group <- trimws(
    as.character(
      meta$nfcore_group[[1]]
    )
  )

  replicate <- trimws(
    as.character(
      meta$replicate[[1]]
    )
  )

  if (
    !nzchar(group) ||
    !nzchar(replicate)
  ) {
    return("")
  }

  # nf-core/cutandrun uses:
  #
  #   <group>_R<replicate>
  #
  # as the canonical per-sample ID.
  #
  # Example:
  #   target__IR_Creb_2wks_R1
  #
  nfcore_id <- paste0(
    group,
    "_R",
    replicate
  )

  basenames <- basename(
    files
  )

  lower_names <- tolower(
    basenames
  )

  lower_id <- tolower(
    nfcore_id
  )

  lower_group <- tolower(
    group
  )

  # First choice: literal nf-core sample ID.
  hit <- grepl(
    lower_id,
    lower_names,
    fixed = TRUE
  )

  # Conservative fallback for any punctuation-normalized filenames.
  if (!any(hit)) {
    file_keys <- vapply(
      lower_names,
      cutrun_nfcore_result_key,
      character(1)
    )

    id_key <- cutrun_nfcore_result_key(
      nfcore_id
    )

    hit <- grepl(
      id_key,
      file_keys,
      fixed = TRUE
    )
  }

  # Last fallback: require both biological group and replicate token.
  if (!any(hit)) {
    group_hit <- grepl(
      lower_group,
      lower_names,
      fixed = TRUE
    )

    rep_patterns <- c(
      paste0("_r", replicate),
      paste0("-r", replicate),
      paste0(".r", replicate),
      paste0("_rep", replicate),
      paste0("_replicate", replicate)
    )

    rep_hit <- Reduce(
      `|`,
      lapply(
        rep_patterns,
        function(pattern) {
          grepl(
            pattern,
            lower_names,
            fixed = TRUE
          )
        }
      )
    )

    hit <- group_hit & rep_hit
  }

  candidates <- files[
    hit
  ]

  if (!length(candidates)) {
    return("")
  }

  candidate_names <- tolower(
    basename(candidates)
  )

  score <- rep(
    0,
    length(candidates)
  )

  # Strong preference for filenames beginning with the exact nf-core ID.
  score[
    startsWith(
      candidate_names,
      lower_id
    )
  ] <- score[
    startsWith(
      candidate_names,
      lower_id
    )
  ] + 1000

  score[
    grepl(
      lower_id,
      candidate_names,
      fixed = TRUE
    )
  ] <- score[
    grepl(
      lower_id,
      candidate_names,
      fixed = TRUE
    )
  ] + 500

  # Read-specific scoring is applied only to the filename portion AFTER
  # <group>_R<replicate>, so replicate R1 is never confused with read R1.
  read <- toupper(
    trimws(
      as.character(read %||% "")
    )
  )

  if (
    read %in% c("R1", "R2")
  ) {
    wanted_read <- if (
      identical(read, "R2")
    ) {
      "2"
    } else {
      "1"
    }

    suffix <- vapply(
      candidate_names,
      function(x) {
        pos <- regexpr(
          lower_id,
          x,
          fixed = TRUE
        )

        if (
          pos[[1]] < 0
        ) {
          return(x)
        }

        substring(
          x,
          pos[[1]] +
            attr(pos, "match.length")
        )
      },
      character(1)
    )

    read_patterns <- if (
      identical(wanted_read, "1")
    ) {
      c(
        "_1_fastqc",
        "_1_val_1_fastqc",
        "_1_val_1",
        "_1."
      )
    } else {
      c(
        "_2_fastqc",
        "_2_val_2_fastqc",
        "_2_val_2",
        "_2."
      )
    }

    read_hit <- Reduce(
      `|`,
      lapply(
        read_patterns,
        function(pattern) {
          grepl(
            pattern,
            suffix,
            fixed = TRUE
          )
        }
      )
    )

    score[read_hit] <-
      score[read_hit] + 250
  }

  prefer <- trimws(
    tolower(
      as.character(
        prefer %||% character(0)
      )
    )
  )

  prefer <- prefer[
    nzchar(prefer)
  ]

  if (length(prefer)) {
    for (term in prefer) {
      score[
        grepl(
          term,
          candidate_names,
          fixed = TRUE
        )
      ] <- score[
        grepl(
          term,
          candidate_names,
          fixed = TRUE
        )
      ] + 50
    }
  }

  best <- which(
    score == max(
      score
    )
  )

  # Never silently choose between equally plausible files.
  if (length(best) != 1L) {
    return("")
  }

  candidates[[best]]
}

cutrun_nfcore_fastqc_path <- function(
  project,
  sample,
  read = c(
    "R1",
    "R2"
  ),
  trimmed = TRUE
) {
  read <- match.arg(read)

  paths <- cutrun_nfcore_results_paths(
    project
  )

  directory <- if (
    isTRUE(trimmed)
  ) {
    paths$trimmed_fastqc
  } else {
    paths$raw_fastqc
  }

  files <- cutrun_nfcore_result_files(
    directory,
    "_fastqc\\.html$"
  )

  cutrun_nfcore_match_sample_file(
    project,
    sample,
    files,
    read = read,
    prefer = c(
      "fastqc"
    )
  )
}


cutrun_nfcore_bam_path <- function(
  project,
  sample
) {
  paths <- cutrun_nfcore_results_paths(
    project
  )

  directories <- c(
    file.path(
      paths$alignment,
      "target",
      "la_duplicates"
    ),
    file.path(
      paths$alignment,
      "target",
      "markdup"
    ),
    file.path(
      paths$alignment,
      "target"
    ),
    paths$alignment
  )

  for (directory in directories) {
    files <- cutrun_nfcore_result_files(
      directory,
      "\\.bam$"
    )

    if (!length(files)) {
      next
    }

    path <- cutrun_nfcore_match_sample_file(
      project,
      sample,
      files,
      prefer = c(
        "la_dedup",
        "markdup"
      )
    )

    if (nzchar(path)) {
      return(path)
    }
  }

  ""
}


cutrun_nfcore_bigwig_path <- function(
  project,
  sample
) {
  paths <- cutrun_nfcore_results_paths(
    project
  )

  files <- cutrun_nfcore_result_files(
    paths$bigwig,
    "\\.(bigwig|bw)$"
  )

  cutrun_nfcore_match_sample_file(
    project,
    sample,
    files
  )
}


cutrun_nfcore_bedgraph_path <- function(
  project,
  sample
) {
  paths <- cutrun_nfcore_results_paths(
    project
  )

  files <- cutrun_nfcore_result_files(
    paths$bedgraph,
    "\\.bedgraph$"
  )

  cutrun_nfcore_match_sample_file(
    project,
    sample,
    files
  )
}


cutrun_nfcore_peak_files <- function(
  project,
  tool = c(
    "SEACR",
    "MACS2"
  )
) {
  tool <- match.arg(tool)

  directory <- cutrun_nfcore_results_paths(
    project
  )$called_peaks

  files <- cutrun_nfcore_result_files(
    directory,
    "\\.(bed|narrowpeak|broadpeak)$"
  )

  if (!length(files)) {
    return(character(0))
  }

  paths_lower <- tolower(files)

  # Consensus peaks are group-level outputs and must never be passed
  # to DiffBind as per-sample peak files.
  files <- files[
    !grepl(
      "consensus",
      paths_lower,
      fixed = TRUE
    )
  ]

  paths_lower <- tolower(files)

  if (
    identical(
      tool,
      "MACS2"
    )
  ) {
    return(
      files[
        grepl(
          "macs",
          paths_lower,
          fixed = TRUE
        ) |
          grepl(
            "narrowpeak",
            paths_lower,
            fixed = TRUE
          ) |
          grepl(
            "broadpeak",
            paths_lower,
            fixed = TRUE
          )
      ]
    )
  }

  seacr <- files[
    grepl(
      "seacr",
      paths_lower,
      fixed = TRUE
    ) |
      grepl(
        "stringent",
        paths_lower,
        fixed = TRUE
      ) |
      grepl(
        "relaxed",
        paths_lower,
        fixed = TRUE
      )
  ]

  if (length(seacr)) {
    return(seacr)
  }

  # nf-core/cutandrun uses SEACR as the primary caller by default.
  # If there are ordinary BED files that clearly are not MACS2,
  # retain them as a conservative SEACR fallback.
  files[
    !grepl(
      "macs",
      paths_lower,
      fixed = TRUE
    ) &
      !grepl(
        "narrowpeak",
        paths_lower,
        fixed = TRUE
      ) &
      !grepl(
        "broadpeak",
        paths_lower,
        fixed = TRUE
      )
  ]
}


cutrun_nfcore_peak_path <- function(
  project,
  sample,
  tool = c(
    "SEACR",
    "MACS2"
  )
) {
  tool <- match.arg(tool)

  files <- cutrun_nfcore_peak_files(
    project,
    tool
  )

  prefer <- if (
    identical(
      tool,
      "SEACR"
    )
  ) {
    c(
      "stringent",
      "seacr"
    )
  } else {
    c(
      "macs"
    )
  }

  cutrun_nfcore_match_sample_file(
    project,
    sample,
    files,
    prefer = prefer
  )
}


cutrun_nfcore_multiqc_path <- function(
  project
) {
  directory <- cutrun_nfcore_results_paths(
    project
  )$multiqc

  expected <- file.path(
    directory,
    "multiqc_report.html"
  )

  if (file.exists(expected)) {
    return(
      normalizePath(
        expected,
        winslash = "/",
        mustWork = FALSE
      )
    )
  }

  files <- cutrun_nfcore_result_files(
    directory,
    "multiqc.*\\.html$"
  )

  if (length(files) == 1L) {
    return(
      normalizePath(
        files[[1]],
        winslash = "/",
        mustWork = FALSE
      )
    )
  }

  ""
}



.cutrun_nfcore_inventory_cache <- new.env(
  parent = emptyenv()
)

cutrun_nfcore_output_inventory <- function(
  project,
  cache_seconds = 15
) {
  paths <- cutrun_nfcore_results_paths(
    project
  )

  key <- normalizePath(
    paths$output_dir,
    winslash = "/",
    mustWork = FALSE
  )

  now <- as.numeric(
    Sys.time()
  )

  if (
    exists(
      key,
      envir = .cutrun_nfcore_inventory_cache,
      inherits = FALSE
    )
  ) {
    cached <- get(
      key,
      envir = .cutrun_nfcore_inventory_cache,
      inherits = FALSE
    )

    age <- now - cached$time

    if (
      is.finite(age) &&
      age >= 0 &&
      age < cache_seconds
    ) {
      return(
        cached$data
      )
    }
  }

  mapping <- cutrun_nfcore_result_mapping(
    project
  )

  if (!NROW(mapping)) {
    empty <- data.frame(
      stringsAsFactors = FALSE
    )

    assign(
      key,
      list(
        time = now,
        data = empty
      ),
      envir = .cutrun_nfcore_inventory_cache
    )

    return(empty)
  }

  rows <- lapply(
    seq_len(
      NROW(mapping)
    ),
    function(i) {
      sample <- as.character(
        mapping$codespring_sample[[i]]
      )

      data.frame(
        sample = sample,

        nfcore_group = as.character(
          mapping$nfcore_group[[i]]
        ),

        replicate = as.character(
          mapping$replicate[[i]]
        ),

        bam = cutrun_nfcore_bam_path(
          project,
          sample
        ),

        bigwig = cutrun_nfcore_bigwig_path(
          project,
          sample
        ),

        bedgraph = cutrun_nfcore_bedgraph_path(
          project,
          sample
        ),

        seacr_peaks = cutrun_nfcore_peak_path(
          project,
          sample,
          "SEACR"
        ),

        macs2_peaks = cutrun_nfcore_peak_path(
          project,
          sample,
          "MACS2"
        ),

        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  result <- do.call(
    rbind,
    rows
  )

  assign(
    key,
    list(
      time = now,
      data = result
    ),
    envir = .cutrun_nfcore_inventory_cache
  )

  result
}

