### load_metaphlan /////////////////////////////////////////////////////////////

test_that("load_metaphlan loads a valid single-sample profile", {
  # Path to test single-sample file (should exist in inst/extdata/)
  file_path <- system.file("extdata", "single_sample_profile.txt", package = "microEDA")

  # Ensure file exists
  expect_true(file.exists(file_path))

  # Load the profile
  result <- load_metaphlan(file_path)

  # Validate object class
  expect_s4_class(result, "metaphlanProfile")

  # Check mpa_version is non-empty (nzchar is a fast way to find out if elements
  # of a character vector are non-empty strings.)
  expect_true(nzchar(result@mpa_version))

  # Check taxonomic profile is a data.frame with expected columns
  expect_s3_class(result@mpa_taxProfile, "data.frame")
  expect_named(result@mpa_taxProfile, c("clade_name", "single_sample_profile"))
})


test_that("load_metaphlan loads a valid merged profile", {
  file_path <- system.file("extdata", "merged_abundance_table.txt", package = "microEDA")

  expect_true(file.exists(file_path))

  result <- load_metaphlan(file_path)

  expect_s4_class(result, "metaphlanProfile")
  expect_true(nzchar(result@mpa_version))
  expect_s3_class(result@mpa_taxProfile, "data.frame")
  expect_true("clade_name" %in% names(result@mpa_taxProfile))
  expect_true(ncol(result@mpa_taxProfile) > 2) # Merged usually has multiple sample columns
})


test_that("load_metaphlan throws error for non-existent file", {
  expect_warning(
    expect_error(
      load_metaphlan("nonexistent.txt"),
      "Please check that the file is a valid MetaPhlAn file"
    ),
    "No such file or directory"
  )
})


test_that("load_metaphlan warns when no database version found", {
  file_path <- system.file("extdata", "no_version_profile.txt", package = "microEDA")

  expect_warning(
    result <- load_metaphlan(file_path),
    "No MetaPhlAn database version found"
  )

  expect_equal(result@mpa_version, character()) # Should be empty
})


test_that("load_metaphlan stops for malformed single-sample file", {
  file_path <- system.file("extdata", "malformed_single.txt", package = "microEDA")

  expect_true(file.exists(file_path))

  expect_error(
    load_metaphlan(file_path),
    "wrong number of columns to be a single MetaPhlAn profile"
  )
})


test_that("load_metaphlan stops for non-Metaphlan file", {
  file_path <- system.file("extdata", "wrong_format.txt", package = "microEDA")

  expect_true(file.exists(file_path))

  expect_warning(
    expect_error(
      load_metaphlan(file_path),
      "Please check that the file is a valid merged MetaPhlAn file"
    ),
    "No MetaPhlAn database version found"
  )
})


test_that("load_metaphlan renames ID to clade_name in metaphlan2 merged profiles", {
  file_path <- system.file("extdata", "merged_v2_format.txt", package = "microEDA")

  expect_true(file.exists(file_path))

  expect_warning(
    result <- load_metaphlan(file_path),
    "No MetaPhlAn database version found"
  )

  expect_s3_class(result@mpa_taxProfile, "data.frame")
  expect_true("clade_name" %in% names(result@mpa_taxProfile))
  expect_false("ID" %in% names(result@mpa_taxProfile))
})

### join_mpa_profiles //////////////////////////////////////////////////////////

# Helper: Create mock metaphlanProfile objects
make_mock_mpa_profile <- function(version, taxa_df) {
  new("metaphlanProfile",
    mpa_version = version,
    mpa_taxProfile = taxa_df
  )
}


test_that("Valid profiles with same DB version are joined correctly", {
  df1 <- data.frame(
    clade_name = c(
      "k__Bacteria|p__Proteobacteria",
      "k__Bacteria|p__Actinobacteria"
    ),
    Sample_1 = c(60, 40), Sample_2 = c(25, 75)
  )
  df2 <- data.frame(
    clade_name = c(
      "k__Bacteria|p__Firmicutes",
      "k__Bacteria|p__Proteobacteria",
      "k__Bacteria|p__Bacteroidota"
    ),
    Sample_3 = c(45, 15, 40)
  )

  prof1 <- make_mock_mpa_profile("#mpa_test_v1.0", df1)
  prof2 <- make_mock_mpa_profile("#mpa_test_v1.0", df2)

  result <- join_mpa_profiles(prof1, prof2)

  expect_s4_class(result, "metaphlanProfile")
  expect_equal(result@mpa_version, "#mpa_test_v1.0")
  expect_equal(colnames(result@mpa_taxProfile), c("clade_name", "Sample_1", "Sample_2", "Sample_3"))
  expect_equal(nrow(result@mpa_taxProfile), 4) # All clades present
  expect_true(all(!is.na(result@mpa_taxProfile))) # No NAs, replaced by 0.0
})


test_that("Profiles from different DB versions throw an error", {
  prof1 <- make_mock_mpa_profile("v3.0", data.frame(clade_name = "A", Sample_1 = 100))
  prof2 <- make_mock_mpa_profile("v4.0", data.frame(clade_name = "B", Sample_2 = 100))

  expect_error(join_mpa_profiles(prof1, prof2), "Profiles from different versions")
})


test_that("Non-metaphlanProfile inputs are dropped with warning", {
  prof1 <- make_mock_mpa_profile("v3.0", data.frame(clade_name = "A", Sample_1 = 1))
  invalid_obj <- list(some_data = 123)

  expect_warning(
    result <- join_mpa_profiles(prof1, invalid_obj),
    'not of class "metaphlanProfile"'
  )
  expect_s4_class(result, "metaphlanProfile")
})


test_that("Empty mpa_version causes an error", {
  prof1 <- make_mock_mpa_profile("v3.0", data.frame(clade_name = "A", Sample_1 = 1))
  prof2 <- make_mock_mpa_profile(character(0), data.frame(clade_name = "B", Sample_2 = 1))

  expect_error(join_mpa_profiles(prof1, prof2), "empty mpa_version")
})


test_that("join_mpa_profiles warns and disambiguates duplicate sample names", {
  df1 <- data.frame(
    clade_name = c(
      "k__Bacteria|p__Proteobacteria",
      "k__Bacteria|p__Actinobacteria"
    ),
    Sample_1 = c(60, 40), Sample_2 = c(25, 75)
  )
  df2 <- data.frame(
    clade_name = c(
      "k__Bacteria|p__Firmicutes",
      "k__Bacteria|p__Proteobacteria",
      "k__Bacteria|p__Bacteroidota"
    ),
    Sample_1 = c(45, 15, 40)
  )

  prof1 <- make_mock_mpa_profile("#mpa_test_v1.0", df1)
  prof2 <- make_mock_mpa_profile("#mpa_test_v1.0", df2)

  expect_warning(
    result <- join_mpa_profiles(prof1, prof2),
    "Duplicate sample names detected: Sample_1"
  )

  expect_s4_class(result, "metaphlanProfile")
  expect_named(result@mpa_taxProfile, c("clade_name", "Sample_1.1", "Sample_2", "Sample_1.2"))
})
