#' @export
#' @title Load Data
#' @description This function will load LC-MS data from a file.
#' @param input A string of the file path or a data frame of the data to be
#' loaded.
#' @param type A string indicating the type of data to be loaded. This can be
#' "known" or "unknown".
#' @param sheet A string or integer indicating the sheet number of the excel
#' file to be loaded.
#' @param id A string indicating the column name of the compound ID.
#' @return A list of the sample metadata, feature metadata, and the data matrix.
load_data <-
  function(input,
           type = "known",
           sheet = 1,
           id = "Metabolite") {
    if (is.character(input)) {
      df <-
        readxl::read_excel(input,
          sheet = sheet,
          col_names = FALSE
        )
    } else {
      df <- input
    }
    df <- as.data.frame(df)
    df[df == "n/a"] <- NA
    df[df == "NA"] <- NA
    df[df == ""] <- NA
    # remove columns with all NAs
    df <- Filter(function(x) {
      !all(is.na(x))
    }, df)

    # find the occurrence of Metabolite word in the file
    ord <- find_indx(df, word = "Metabolite")
    row <- ord[1]
    col <- ord[2]

    compound_ord <- find_indx(df, word = "Compound_ID")
    hmdb_ord <- find_indx(df, word = "HMDB_ID")
    mz_ord <- find_indx(df, word = "MZ")
    rt_ord <- find_indx(df, word = "RT")
    # separate sample metadata
    col_names <- NA
    if (row > 1) {
      sample_metadata <- df[1:(row - 1), col:dim(df)[2]]
      meta_rownames <- sample_metadata[[1]]
      sample_metadata <-
        sample_metadata[1:dim(sample_metadata)[1], 2:dim(sample_metadata)[2]]
      col_names <- sapply(df[row, (col + 1):dim(df)[2]], as.factor)
      colnames(sample_metadata) <- col_names
      rownames(sample_metadata) <- meta_rownames
    } else {
      col_names <- sapply(df[row, (col + 1):dim(df)[2]], as.factor)
      sample_metadata <- NA
    }
    # separate feature metadata
    f_col_names <- NA
    if (col > 1) {
      feature_metadata <- df[row:dim(df)[1], 1:col]
      f_meta_rownames <- feature_metadata[[compound_ord[2]]]
      feature_metadata <-
        feature_metadata[2:dim(feature_metadata)[1], 1:dim(feature_metadata)[2]]
      f_col_names <- sapply(df[row, 1:col], as.factor)
      colnames(feature_metadata) <- f_col_names
      rownames(feature_metadata) <-
        f_meta_rownames[-1] # ro remove "Compound" from the names
    } else {
      feature_metadata <- NA
    }
    df <- as.data.frame(df)


    # separate known data or generate a ID for unknows in the case of using all
    if (type == "known") {
      if (id == "Metabolite") {
        df <- df[!is.na(df[, col]),]
      } else if (id == "HMDB_ID") {
        df <- df[!is.na(df[, hmdb_ord[2]]),]
        df <- df[!grepl("\\*", df[, hmdb_ord[2]]),]
        df <- df[df[, hmdb_ord[2]] != "",]
      }
    } else if (type == "all") {
      if (id == "Compound_ID") {
        name_col <- df[, compound_ord[2]] # df[,col-1]
      } else {
        if (id == "Metabolite") {
          name_col <- col
        } else if (id == "HMDB_ID") {
          name_col <- hmdb_ord[2]
        }
        df[, col] <- ifelse(is.na(df[, name_col]),
          paste(
            df[, compound_ord[2]],
            "MZ",
            round(as.numeric(df[, mz_ord[2]]), digits = 4),
            "RT",
            round(as.numeric(df[, rt_ord[2]]), digits = 2),
            sep = "_"
          ),
          as.character(df[, name_col])
        )
      }
    }

    # remove redundant ion
    df[is.na(df)] <- "NA"
    if (!is.na(hmdb_ord[2])) {
      df <- df[which(df[, hmdb_ord[2]] != "redundant ion"), ]
      df <- df[df[, hmdb_ord[2]] != "Internal Standard", ]
    }
    df <- df[df[, col] != "Metabolite", ]
    # remove feature start with NH4_
    df <- df[!startsWith(sapply(df[, col], as.character), "NH4_"), ]
    df[df == "NA"] <- NA
    df <- df[!is.na(df[, col]), ]

    if (id == "HMDB_ID") {
      row_names <-
        gsub("\\*", "", as.matrix(df[row:dim(df)[1], hmdb_ord[2]]))
    } else if (id == "Metabolite") {
      row_names <- unlist(sapply(df[row:dim(df)[1], col], as.factor))
    } else if (id == "Compound_ID") {
      row_names <-
        unlist(sapply(df[row:dim(df)[1], compound_ord[2]], as.factor))
    }
    data <- df[row:dim(df)[1], col:dim(df)[2]]
    data <- data[, -1]
    # check if row_names has unique values to be used for row names
    dup_names <- row_names[duplicated(row_names)]
    if (length(dup_names) > 0) {
      print(sprintf("There are duplicated values: %s", dup_names))
    }
    rownames(data) <- row_names
    colnames(data) <- col_names
    data <- as.data.frame(data)
    result <- list()
    # put data in the right ordination rows(observation or samples) columns features
    result$sample_metadata <- as.data.frame(t(sample_metadata))
    result$feature_metadata <- feature_metadata
    data <- as.data.frame(t(data))

    # make sure data frame is numeric
    result$data <- numeric_dataframe(data)
    return(result)
  }


#' @export
combine_QI_TF <- function(QI_file, TF_file, output_name) {
  # read the Progenesis file
  qi_data <- read.table(
    QI_file,
    header = F,
    row.names = NULL,
    sep = ",",
    fill = FALSE,
    comment.char = "",
    check.names = FALSE
  )

  # name the columns
  colnames(qi_data) <- lapply(qi_data[3, ], as.character)

  # find the sart of raw abundance data
  data_indx <- find_indx(qi_data, word = "Raw abundance")

  # remove the first two columns
  qi_data <- qi_data[4:nrow(qi_data), ]

  # use raw abundance
  qi_data <- qi_data[, c(1, 3, 5, data_indx[2]:ncol(qi_data))]

  # re-order columns and remove columns with name ""
  indx_last_sample <-
    grep("Accepted Compound ID", colnames(qi_data))
  cols <- colnames(qi_data)
  cols <- cols[1:indx_last_sample]
  samples_cols <-
    cols[!(cols %in% c(
      "Compound",
      "m/z",
      "Retention time (min)",
      "Accepted Compound ID",
      ""
    ))]
  samples_cols <- samples_cols[order(samples_cols)]
  cols_in_order <-
    c(
      c(
        "Compound",
        "m/z",
        "Retention time (min)",
        "Accepted Compound ID"
      ),
      samples_cols
    )
  qi_data <- qi_data[, cols_in_order]


  ##### read first file for TR profiles #########################
  # 1) on the second tab, I calculate the average retention time

  RT_profile_data <-
    readxl::read_excel(TF_file,
      sheet = 2,
      col_names = FALSE
    )
  colnames(RT_profile_data) <-
    sapply(RT_profile_data[1, ], as.factor)
  RT_profile_data_rownames <-
    sapply(RT_profile_data[, 1], as.factor)
  RT_profile_data <- RT_profile_data[-1, -1]
  rownames(RT_profile_data) <- RT_profile_data_rownames[-1]
  RT_profile_data <- as.data.frame(RT_profile_data)

  ##### Calculate average RT
  RT_profile_data <- omicsArt::numeric_dataframe(RT_profile_data)
  RT_profile_data[RT_profile_data <= 0.0] <- NA
  RT <- colMeans(x = RT_profile_data, na.rm = T)
  RT_profile_data <- rbind(RT = RT, RT_profile_data)

  # 2) I add the RT to the first tab, then copy & transpose to a new tab

  ##### read first file for TR profiles #########################
  intensity_profile_data <-
    readxl::read_excel(TF_file,
      sheet = 1,
      col_names = FALSE
    )
  colnames(intensity_profile_data) <-
    sapply(intensity_profile_data[1, ], as.factor)
  intensity_profile_data_rownames <-
    sapply(intensity_profile_data[, 1], as.factor)
  intensity_profile_data <- intensity_profile_data[-1, -1]
  rownames(intensity_profile_data) <-
    intensity_profile_data_rownames[-1]
  intensity_profile_data <- as.data.frame(intensity_profile_data)

  ##### add the average RT
  intensity_profile_data <-
    omicsArt:::numeric_dataframe(intensity_profile_data)
  intensity_profile_data <- rbind(RT = RT, intensity_profile_data)

  # clean data
  ## remove "Standards" mm, Bpp, FFA, miniMM
  clean_rows <- row.names(intensity_profile_data)
  clean_rows <-
    clean_rows[!grepl("*miniMM*|*FFA*|*Bpp*|*mm*|*Standards*", clean_rows)]
  intensity_profile_clean_data <-
    intensity_profile_data[clean_rows, ]
  tf_data <- as.data.frame(t(intensity_profile_clean_data))
  tf_data[
    ,
    c(
      "Compound",
      "m/z",
      "Retention time (min)",
      "Accepted Compound ID"
    )
  ] <- NA
  tf_data[, "Retention time (min)"] <- tf_data[, "RT"]

  # rename TF metabolites if they are in QI data
  row.names(tf_data) <-
    ifelse(
      row.names(tf_data) %in% row.names(qi_data),
      paste(row.names(tf_data), "TF", sep = "_"),
      row.names(tf_data)
    )

  tf_data[, "Software"] <- "TF"
  tf_data$`Accepted Compound ID` <- rownames(tf_data)
  qi_data[, "Software"] <- "QI"
  combined <- rbind(tf_data[, colnames(qi_data)], qi_data)

  # order the final column names
  cols_in_order <-
    c(
      c(
        "Compound",
        "Software",
        "m/z",
        "Retention time (min)",
        "
      Accepted Compound ID"
      ),
      samples_cols
    )
  combined <- combined[, cols_in_order]

  colnames(combined) <- c(
    c("Compound_ID", "Software", "MZ", "RT", "Metabolite"),
    samples_cols
  )
  combined[combined == ""] <- NA
  combined <-
    with(combined, combined[order(Metabolite, na.last = TRUE), ])
  hs <- openxlsx::createStyle(
    textDecoration = "BOLD",
    fontColour = "#FFFFFF",
    fontSize = 12,
    fontName = "Arial Narrow",
    fgFill = "#4F80BD"
  )
  options("openxlsx.borderColour" = "#4F80BD") ## set default border colour
  xlsx::write.xlsx(combined,
    file = output_name,
    colNames = TRUE
  )
  # borders = "rows",
  # headerStyle = hs)
}
