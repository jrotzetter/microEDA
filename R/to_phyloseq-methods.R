#' Convert an object to a phyloseq object
#'
#' S4 generic to convert objects into a `phyloseq` object.
#' Currently, methods for `metaphlanProfile` and `microEDA` objects are implemented.
#'
#' @param x An object to be converted to `phyloseq` class.
#' @param ... Additional arguments passed to methods (e.g., `metadata`, `sample_column`).
#' @return An object of class `phyloseq`.
#' @seealso \code{\link[phyloseq]{phyloseq}}, \code{\link{microEDA}}, \code{\link{metaphlanProfile-class}}
#' @name to_phyloseq
#' @export
#' @examples
#' # Convert a microEDA object to phyloseq
#' data(enterotype, package = "phyloseq")
#' me <- microEDA(enterotype)
#' is(me, "microEDA")
#' ps <- to_phyloseq(me)
#' is(ps, "phyloseq")
#'
#' # Convert a metaphlanProfile object to phyloseq
#' mpa <- microEDA::merged_metaphlan_profiles
#' is(mpa, "metaphlanProfile")
#' ps <- to_phyloseq(mpa)
#' is(ps, "phyloseq")
setGeneric("to_phyloseq", function(x, ...) standardGeneric("to_phyloseq"))


#' Convert microEDA object to phyloseq
#'
#' @param x An object of class `microEDA`.
#' @export
#' @rdname to_phyloseq
setMethod(
  "to_phyloseq", "microEDA",
  function(x) {
    ps <- new("phyloseq",
      otu_table = x@otu_table,
      tax_table = x@tax_table,
      sam_data = x@sam_data,
      phy_tree = x@phy_tree,
      refseq = x@refseq
    )
    validObject(ps)
    return(ps)
  }
)


#' Convert metaphlanProfile object to phyloseq
#' @param metadata Optional `data.frame` containing sample metadata.
#' @param sample_column `Character` string specifying the column in `metadata` with sample IDs.
#' @rdname to_phyloseq
#' @aliases to_phyloseq,metaphlanProfile-method
#' @export
setMethod(
  "to_phyloseq",
  signature(x = "metaphlanProfile"),
  function(x, metadata = NULL, sample_column = NULL) {
    # Use internal constructor to build microEDA object
    me <- .metaphlanConstructor(tax_profile = x, metadata = metadata, sample_column = sample_column)
    if (!validObject(me)) stop("Constructed microEDA object is invalid.")

    # Extract components for phyloseq
    ps <- phyloseq::phyloseq(
      otu_table = me@otu_table,
      tax_table = me@tax_table,
      sam_data = me@sam_data
    )

    # Validate and return
    validObject(ps)
    return(ps)
  }
)
