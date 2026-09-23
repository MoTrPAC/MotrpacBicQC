context("Test misc functions")

test_that("remove empty columns works", {
  L3 <- LETTERS[1:3]
  fac <- sample(L3, 10, replace = TRUE)
  df <- data.frame(x = 1, y = 1:10, fac = fac, stringsAsFactors = FALSE)

  # Add empty column
  df$tremove1 <- ""
  df$tremove2 <- NA

  expect_equal(dim(remove_empty_columns(df, verbose = FALSE))[2], 3)
})

test_that("remove empty rows works", {
  L3 <- LETTERS[1:3]
  fac <- sample(L3, 10, replace = TRUE)
  df <- data.frame(x = 1, y = 1:10, fac = fac, stringsAsFactors = FALSE)

  # Add empty rows
  df[nrow(df)+1,] <- NA
  df[nrow(df)+1,] <- ""

  expect_equal(dim(remove_empty_rows(df, verbose = FALSE))[1], 10)
})

test_that("Check all columns in data table bic_animal_tissue_code", {
  expected_colnames <- c("bic_tissue_code", "bic_tissue_name", "motrpac_tissue_code", "tissue_name_release", "abbreviation", "tissue_hex_colour")
  expect_equal( setequal(expected_colnames, colnames(bic_animal_tissue_code)), TRUE)
  expect_gt(dim(bic_animal_tissue_code)[1], 40)
})

test_that("Whitespace is trimmed in character columns", {
  df <- data.frame(
    a = c("  hello", "world  ", "  foo  "),
    b = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  cleaned <- clean_character_columns(df)
  expect_equal(cleaned$a, c("hello", "world", "foo"))
  expect_equal(cleaned$b, c(1, 2, 3))
})

test_that("Function leaves columns without extra whitespace unchanged", {
  df <- data.frame(
    a = c("apple", "banana", "cherry"),
    b = c("dog", "cat", "bird"),
    stringsAsFactors = FALSE
  )
  cleaned <- clean_character_columns(df)
  expect_equal(cleaned$a, c("apple", "banana", "cherry"))
  expect_equal(cleaned$b, c("dog", "cat", "bird"))
})

test_that("Works with empty data frame", {
  df <- data.frame(
    a = character(0),
    b = numeric(0),
    stringsAsFactors = FALSE
  )
  cleaned <- clean_character_columns(df)
  expect_equal(nrow(cleaned), 0)
})

test_that("valid_sample_types is a non-empty character vector with expected types", {
  expect_true(is.character(valid_sample_types))
  expect_true(length(valid_sample_types) >= 13)
  # Core types that must always be present
  expect_true("Sample" %in% valid_sample_types)
  expect_true("QC-Pooled" %in% valid_sample_types)
  expect_true("QC-Reference" %in% valid_sample_types)
  expect_true("QC-Reference-Male" %in% valid_sample_types)
  expect_true("QC-Reference-Female" %in% valid_sample_types)
  expect_true("QC-Blank" %in% valid_sample_types)
  expect_true("QC-PlateControl" %in% valid_sample_types)
  # No duplicates
  expect_equal(length(valid_sample_types), length(unique(valid_sample_types)))
})


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot_na_percentage(): replacement of inspectdf::inspect_na() %>% show_plot()

test_that("plot_na_percentage computes the % of NA per column in column order", {
  df <- data.frame(zeta = c(1, NA, NA, NA),
                   alpha = c(1, 2, 3, 4),
                   mid = c(NA, NA, 1, 2),
                   txt = c("a", NA, "c", "d"),
                   stringsAsFactors = FALSE)
  p <- plot_na_percentage(df)
  expect_s3_class(p, "ggplot")

  # Rows follow the data frame columns, NOT the % of NA (naniar sorts by default)
  expect_equal(as.character(p$data$variable), colnames(df))
  expect_equal(levels(p$data$variable), colnames(df))
  expect_equal(p$data$pct_miss, c(75, 0, 50, 25))
  expect_equal(p$data$n_miss, c(3L, 0L, 2L, 1L))

  # The subtitle reports the number of columns with missing values
  expect_match(p$labels$subtitle, "4 columns, of which 3 have missing values")
})

test_that("plot_na_percentage text labels can be switched off", {
  df <- data.frame(a = c(1, NA), b = c(NA, NA))
  # one bar layer + label layers (inside/above the bars)
  expect_gt(length(plot_na_percentage(df, text_labels = TRUE)$layers), 1)
  expect_length(plot_na_percentage(df, text_labels = FALSE)$layers, 1)
})

test_that("plot_na_percentage handles edge cases", {
  # No missing values at all, and a fully missing column
  expect_equal(plot_na_percentage(data.frame(a = 1:3, b = 4:6))$data$pct_miss, c(0, 0))
  expect_equal(plot_na_percentage(data.frame(a = c(NA, NA)))$data$pct_miss, 100)
  expect_s3_class(plot_na_percentage(data.frame(a = c(1, NA))), "ggplot")

  expect_error(plot_na_percentage(list(a = 1)), "must be a data.frame")
  expect_error(plot_na_percentage(c(1, NA)), "must be a data.frame")
  expect_error(plot_na_percentage(data.frame()), "no columns")
})

test_that("plot_na_percentage composes with the layers added by the QC plots", {
  # Exactly the chains used in metabolomics/olink (text_labels = FALSE) and
  # proteomics (default) after the replacement of inspectdf
  df <- data.frame(a = c(1, NA), b = c(NA, NA))
  p1 <- plot_na_percentage(df, text_labels = FALSE) + ggplot2::ylim(0, 100) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Prevalence of NAs", subtitle = "prefix", y = "% of NAs in each sample")
  expect_s3_class(p1, "ggplot")
  expect_equal(p1$labels$title, "Prevalence of NAs")
  expect_equal(p1$labels$y, "% of NAs in each sample")

  p2 <- plot_na_percentage(df) + ggplot2::ylim(0, 100) + ggplot2::theme_linedraw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
  expect_s3_class(p2, "ggplot")
  # Both must be printable (this is what ends up in the pdf)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(print(p1))
  expect_no_error(print(p2))
})

test_that("plot_na_percentage on a column subset keeps the original column order", {
  # The proteomics QC plots the required columns of the rii/ratio files with
  # `df[intersect(colnames(df), required_columns)]`, which must preserve the
  # order of the file (as `arrange(match(col_name, colnames(df)))` used to do)
  df <- data.frame(protein_id = c("P1", "P2"), s_b = c(NA, 1), s_a = c(1, NA), extra = c(1, 1))
  required_columns <- c("protein_id", "s_a", "s_b")
  p <- plot_na_percentage(df[intersect(colnames(df), required_columns)])
  expect_equal(as.character(p$data$variable), c("protein_id", "s_b", "s_a"))
  expect_equal(p$data$pct_miss, c(0, 50, 50))
})

test_that("plot_na_percentage works on the bundled metabolomics results", {
  p <- plot_na_percentage(results_named)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), ncol(results_named))
  expect_equal(as.character(p$data$variable), colnames(results_named))
  expect_true(all(p$data$pct_miss >= 0 & p$data$pct_miss <= 100))
})
