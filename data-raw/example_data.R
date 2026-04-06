## code to prepare example datasets goes here

path_single_profile <- system.file("extdata",
  "single_sample_profile.txt",
  package = "microEDA",
  mustWork = TRUE
)

single_metaphlan_profile <- load_metaphlan(path_single_profile)

usethis::use_data(single_metaphlan_profile, overwrite = TRUE)


path_merged_profiles <- system.file("extdata",
  "merged_abundance_table.txt",
  package = "microEDA",
  mustWork = TRUE
)

merged_metaphlan_profiles <- load_metaphlan(path_merged_profiles)

usethis::use_data(merged_metaphlan_profiles, overwrite = TRUE)


# Source for merged MetaPhlAn2 profiles data:
# https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Tutorial-(Long-Running)
# path_v2_profiles <- system.file("extdata",
#   "merged_v2_format.txt",
#   package = "microEDA",
#   mustWork = TRUE
# )
#
# merged_v2_profiles <- load_metaphlan(path_v2_profiles)
#
# usethis::use_data(merged_v2_profiles, internal = TRUE)
