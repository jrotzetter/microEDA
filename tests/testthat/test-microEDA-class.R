test_that("microEDA throws error for invalid input", {
  # Test with NULL input
  expect_error(
    microEDA(NULL),
    "'tax_profile' must be a metaphlanProfile or phyloseq object"
  )

  # Test with wrong class input
  expect_error(
    microEDA(123),
    "'tax_profile' must be a metaphlanProfile or phyloseq object"
  )

  # Test with character input
  expect_error(
    microEDA("not_a_phyloseq_object"),
    "'tax_profile' must be a metaphlanProfile or phyloseq object"
  )
})

### metaphlanProfile ///////////////////////////////////////////////////////////

# Helper: Create mock metaphlanProfile objects
mock_metaphlanProfile <- function(version, taxa_df) {
  new("metaphlanProfile",
    mpa_version = version,
    mpa_taxProfile = taxa_df
  )
}


test_that("microEDA constructor works with metaphlanProfile object", {
  taxa_df <- data.frame(
    clade_name = c(
      "k__Bacteria",
      "k__Bacteria|p__Firmicutes",
      "k__Bacteria|p__Proteobacteria",
      "k__Bacteria|p__Bacteroidota"
    ),
    Sample_1 = c(100, 45.65, 24.35, 30),
    Sample_2 = c(100, 25, 40, 35),
    Sample_3 = c(100, 45, 15, 40)
  )
  version <- "#mpa_test_v1.0"

  mtp <- mock_metaphlanProfile(version, taxa_df)
  me <- microEDA(mtp)

  # Check class inheritance
  expect_s4_class(me, "microEDA")
  expect_s4_class(me, "phyloseq")

  # Check for the presence of 'info' slot
  expect_true("info" %in% slotNames(me))
  expect_type(slot(me, "info"), "list")

  expect_equal(me@info$taxrank, "Phylum")
  expect_equal(me@info$mpa_version, "#mpa_test_v1.0")
  expect_equal(me@info$transforms, "relabund")
})


test_that("throw error with metaphlanProfile object when only one taxonomic rank", {
  taxa_df <- data.frame(
    clade_name = c(
      "k__Bacteria|p__Firmicutes",
      "k__Bacteria|p__Proteobacteria",
      "k__Bacteria|p__Bacteroidota"
    ),
    Sample_1 = c(46, 24, 30),
    Sample_2 = c(25, 40, 35),
    Sample_3 = c(45, 15, 40)
  )
  version <- "#mpa_test_v1.0"

  mtp <- mock_metaphlanProfile(version, taxa_df)

  expect_error(microEDA(mtp), "Only one taxonomic rank found in MetaPhlAn profile")
})

# TODO: Add constructor test for metadata

### phyloseq ///////////////////////////////////////////////////////////////////

# Helper: Create simple mock phyloseq object
mock_phyloseq <- function(taxa_are_rows = TRUE) {
  if (taxa_are_rows) {
    otu <- phyloseq::otu_table(
      matrix(sample(1:100, 20, replace = TRUE),
        nrow = 4,
        dimnames = list(paste0("OTU_", 1:4), paste0("Sample_", 1:5))
      ),
      taxa_are_rows = TRUE
    )
  } else {
    otu <- phyloseq::otu_table(
      t(matrix(sample(1:100, 20, replace = TRUE),
        nrow = 4,
        dimnames = list(paste0("OTU_", 1:4), paste0("Sample_", 1:5))
      )),
      taxa_are_rows = FALSE
    )
  }

  sam <- phyloseq::sample_data(
    data.frame(
      sample_names = paste0("Sample_", 1:5),
      stringsAsFactors = FALSE,
      row.names = paste0("Sample_", 1:5)
    )
  )
  tax <- phyloseq::tax_table(matrix(
    c(
      "Bacteria", "Bacteroidetes", "Bacteroidia",
      "Bacteria", "Firmicutes", "Bacilli",
      "Bacteria", "Proteobacteria", "Gammaproteobacteria",
      "Archaea", "Methanobacteriota", "Methanomicrobia"
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("OTU_", 1:4), c("Kingdom", "Phylum", "Class"))
  ))
  phyloseq::phyloseq(otu, sam, tax)
}


test_that("microEDA can be constructed from phyloseq object", {
  ps <- mock_phyloseq()
  me <- microEDA(ps)

  # Check class inheritance
  expect_s4_class(me, "microEDA")
  expect_s4_class(me, "phyloseq")

  # Check for the presence of 'info' slot
  expect_true("info" %in% slotNames(me))
  expect_type(slot(me, "info"), "list")

  # Check that all phyloseq components are preserved
  expect_s4_class(phyloseq::otu_table(me), "otu_table")
  expect_s4_class(phyloseq::tax_table(me), "taxonomyTable")
  expect_s4_class(phyloseq::sample_data(me), "sample_data")
  # expect_s3_class(phyloseq::phy_tree(me), "phylo")
  expect_error(phyloseq::phy_tree(me), "phy_tree slot is empty.")
  expect_equal(phyloseq::phy_tree(me, errorIfNULL = FALSE), NULL)
})


test_that("microEDA can be constructed from phyloseq object where taxa_are_rows = FALSE", {
  ps <- mock_phyloseq(taxa_are_rows = FALSE)
  expect_equal(phyloseq::taxa_are_rows(ps), FALSE)

  me <- microEDA(ps)
  expect_equal(phyloseq::taxa_are_rows(ps), FALSE)
  expect_equal(phyloseq::taxa_are_rows(me), TRUE)

  # Check class inheritance
  expect_s4_class(me, "microEDA")
  expect_s4_class(me, "phyloseq")

  # Check for the presence of 'info' slot
  expect_true("info" %in% slotNames(me))
  expect_type(slot(me, "info"), "list")

  # Check that all phyloseq components are preserved
  expect_s4_class(phyloseq::otu_table(me), "otu_table")
  expect_s4_class(phyloseq::tax_table(me), "taxonomyTable")
  expect_s4_class(phyloseq::sample_data(me), "sample_data")
  # expect_s3_class(phyloseq::phy_tree(me), "phylo")
  expect_error(phyloseq::phy_tree(me), "phy_tree slot is empty.")
  expect_equal(phyloseq::phy_tree(me, errorIfNULL = FALSE), NULL)

})


test_that("info slot is properly initialized", {
  ps <- mock_phyloseq()
  me <- microEDA(ps)

  # Check that info slot exists
  expect_true("info" %in% slotNames(me))
  expect_type(slot(me, "info"), "list")

  # Check we can add information to the info slot
  me@info$transforms <- c("log", "relabund")

  expect_equal(me@info$taxrank, "Class")
  expect_equal(me@info$transforms, c("log", "relabund"))
})


test_that("phyloseq methods work on microEDA objects", {
  # Load test data
  data("GlobalPatterns", package = "phyloseq")
  me <- microEDA(GlobalPatterns)

  # Test basic phyloseq functions
  expect_equal(phyloseq::nsamples(me), phyloseq::nsamples(GlobalPatterns))
  expect_equal(phyloseq::ntaxa(me), phyloseq::ntaxa(GlobalPatterns))

  # Test that subsetting works
  sample_ids <- rownames(phyloseq::sample_data(me))[1:5] # First 5 samples
  me_subset_samples <- phyloseq::prune_samples(sample_ids, me)

  expect_s4_class(me_subset_samples, "microEDA")
  expect_equal(phyloseq::nsamples(me_subset_samples), 5)

  me_subset_taxa <- phyloseq::prune_taxa(phyloseq::taxa_names(me)[1:5], me) # First 5 taxa

  expect_s4_class(me_subset_taxa, "microEDA")
  expect_equal(phyloseq::ntaxa(me_subset_taxa), 5)
})


test_that("constructed microEDA object preserves data", {
  me <- microEDA(GlobalPatterns)

  # Check that the data is the same as original
  expect_equal(sum(phyloseq::otu_table(me)), sum(phyloseq::otu_table(GlobalPatterns)))

  # Check that sample names are preserved
  expect_equal(phyloseq::sample_names(me), phyloseq::sample_names(GlobalPatterns))

  # Check that taxonomic names are preserved
  expect_equal(phyloseq::taxa_names(me), phyloseq::taxa_names(GlobalPatterns))
})

# TODO: Add constructor test for metadata
