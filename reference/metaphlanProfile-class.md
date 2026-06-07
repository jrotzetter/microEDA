# A class to represent a MetaPhlAn profile

This class stores database version information and taxonomic profile
data from MetaPhlAn output.

## Slots

- `mpa_version`:

  `Character` string indicating the MetaPhlAn database version used.

- `mpa_taxProfile`:

  A `data.frame` containing the taxonomic profile (relative abundances
  of clades).
