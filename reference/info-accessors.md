# Get or set the info slot from a microEDA object

Accessor methods for either retrieving or setting the info slot in a
[microEDA](https://jrotzetter.github.io/microEDA/reference/microEDA-class.md)
object. Additionally, each expected element of the info slot can be
retrieved or set individually. Valid keys for the `microEDA` info slot
can be displayed with `infoKeys()`, while filled info fields of a
`microEDA` object can be displayed with `infoFields()`.

## Usage

``` r
info(object)

info(object) <- value

infoKeys()

infoFields(object)

taxrank(object)

taxrank(object) <- value

transforms(object)

transforms(object) <- value

# S4 method for class 'microEDA'
mpa_version(object)

# S4 method for class 'microEDA'
mpa_version(object) <- value

filtered_taxa(object)

filtered_taxa(object) <- value

filter_history(object)

filter_history(object) <- value
```

## Arguments

- object:

  A
  [microEDA](https://jrotzetter.github.io/microEDA/reference/microEDA-class.md)
  object.

- value:

  Only used in replacement methods (setters).

  For `info<-`:

  :   A `list` containing metadata. Valid metadata keys are: 'taxrank',
      'transforms', 'filter_history', 'mpa_version'.

  For `taxrank<-`:

  :   A `character` string specifying the taxonomic rank (e.g.,
      "Species", "Genus").

  For `transforms<-`:

  :   A `character` vector of applied transformation(s) (e.g., c("TSS",
      "log2")).

  For `mpa_version<-`:

  :   A `character` string specifying the MetaPhlAn version (e.g.,
      "#mpa_vJun23_CHOCOPhlAnSGB_202307").

  For `filtered_taxa<-`:

  :   A `list` with otu_table and tax_table (as matrices) of features
      that were filtered out.

  For `filter_history<-`:

  :   A named list specifying used filter parameters. Expected elements
      include:

      - `min_abundance`: numeric threshold for minimum abundance.

      - `min_prevalence`: numeric threshold for minimum prevalence.

      - `group_var`: optional grouping variable for stratified
        filtering.

      - `abundance_criterion`: filtering criterion applied to abundance
        values.

      - `group_requirement`: filter passing requirement for group-wise
        filtering.

      - `keep_filtered`: logical indicating whether filtered-out
        features were preserved in `@info$filtered_taxa`.

      - `n_removed`: integer indicating the number of features removed.

      See
      [`filter_features`](https://jrotzetter.github.io/microEDA/reference/filter_features.md)
      for details on filtering arguments.

## Value

- For `info()`: :

  A `list` containing metadata.

- For `infoKeys()`: :

  Lists valid keys for `microEDA` info slot.

- For `infoFields()`: :

  Lists available field(s) in `object` info slot.

- For `taxrank()`: :

  A `character` string with the lowest taxonomic rank.

- For `transforms()`: :

  A `character` vector of applied transformations.

- For
  [`mpa_version()`](https://jrotzetter.github.io/microEDA/reference/mpa_version.md):
  :

  A `character` string with the used the MetaPhlAn database version.

- For `filtered_taxa()`: :

  A `list` with otu_table and tax_table of features that were filtered
  out.

- For `filter_history()`: :

  A `list` of used filter parameters.

Replacement methods (`<-`) will return a `microEDA` object with the
updated corresponding slots.

## Details

`taxrank`: represents the lowest available taxonomic rank.

`transforms`: represents transformation(s) that were applied to the
abundance data (otu_table).

`mpa_version`: represents the MetaPhlAn database version that was used
to create the profile.

`filtered_taxa`: represents the otu_table and tax_table of features that
were filtered out.

`filter_history`: represents filters that were applied to the abundance
data (otu_table).

## See also

[`mpa_version,metaphlanProfile-method`](https://jrotzetter.github.io/microEDA/reference/mpa_version.md)
