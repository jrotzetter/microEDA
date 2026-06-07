# Convert an object to a phyloseq object

S4 generic to convert objects into a `phyloseq` object. Currently,
methods for `metaphlanProfile` and `microEDA` objects are implemented.

## Usage

``` r
to_phyloseq(x, ...)

# S4 method for class 'microEDA'
to_phyloseq(x)

# S4 method for class 'metaphlanProfile'
to_phyloseq(x, metadata = NULL, sample_column = NULL)
```

## Arguments

- x:

  An object of class `microEDA`.

- ...:

  Additional arguments passed to methods (e.g., `metadata`,
  `sample_column`).

- metadata:

  Optional `data.frame` containing sample metadata.

- sample_column:

  `Character` string specifying the column in `metadata` with sample
  IDs.

## Value

An object of class `phyloseq`.

## See also

[`phyloseq`](https://rdrr.io/pkg/phyloseq/man/phyloseq.html),
[`microEDA`](https://jrotzetter.github.io/microEDA/reference/microEDA-class.md),
[`metaphlanProfile-class`](https://jrotzetter.github.io/microEDA/reference/metaphlanProfile-class.md)

## Examples

``` r
# Convert a microEDA object to phyloseq
data(enterotype, package = "phyloseq")
me <- microEDA(enterotype)
is(me, "microEDA")
#> [1] TRUE
ps <- to_phyloseq(me)
is(ps, "phyloseq")
#> [1] TRUE

# Convert a metaphlanProfile object to phyloseq
mpa <- microEDA::merged_metaphlan_profiles
is(mpa, "metaphlanProfile")
#> [1] TRUE
ps <- to_phyloseq(mpa)
is(ps, "phyloseq")
#> [1] TRUE
```
