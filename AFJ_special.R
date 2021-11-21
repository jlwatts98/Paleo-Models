##### Climate Change Special AFJ #####
source("Header.R")

# load in New Zealand data
nzferns = readr::read_csv("nz_cyatheales/nz_cyatheales.csv")
# look at species
species = unique(nzferns$species)
species
# split dataframe into individual species
indspecies = split(nzferns, nzferns$species)

# remove species below a certain threshold for occurrences
finalspecies = indspecies[sapply(indspecies, function(x) dim(x)[1]) > 200] |>
  purrr::map(~ rename(., lat = decimalLatitude)) |> # rename decimallongitude and decimallatitude
  purrr::map(~ rename(., lon = decimalLongitude))
finalspecies = lapply(finalspecies,as.data.frame)

# clean coordinates with coordinate clean function
# load predictors
predictors = readr::read_rds("objects/predictors.rds")
cleaned = lapply(finalspecies, clean_pres, predictors, thin = TRUE) 

# write rds object
readr::write_rds(cleaned, "objects/cleaned_pres.rds")

