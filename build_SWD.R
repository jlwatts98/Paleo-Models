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

# read rds object
cleaned_pres = readr::read_rds("objects/cleaned_pres.rds")

# get background points
backgr = lapply(cleaned_pres, get_backgr, 200000,
                    5000, predictors)
# write rds object
readr::write_rds(backgr, "objects/backgr.rds")

# read rds object
backgr = readr::read_rds("objects/backgr.rds")
species = names(backgr)

# prepare SWD object for SDMtune

s1 = prepareSWD(species[[1]], p = cleaned_pres[[1]], a = backgr[[1]], 
           env = predictors)
s2 = prepareSWD(species[[2]], p = cleaned_pres[[2]], a = backgr[[2]], 
           env = predictors)
s3 = prepareSWD(species[[3]], p = cleaned_pres[[3]], a = backgr[[3]], 
                env = predictors)
s4 = prepareSWD(species[[4]], p = cleaned_pres[[4]], a = backgr[[4]], 
                env = predictors)
s5 = prepareSWD(species[[5]], p = cleaned_pres[[5]], a = backgr[[5]], 
                env = predictors)
s6 = prepareSWD(species[[6]], p = cleaned_pres[[6]], a = backgr[[6]], 
                env = predictors)
s7 = prepareSWD(species[[7]], p = cleaned_pres[[7]], a = backgr[[7]], 
                env = predictors)
s8 = prepareSWD(species[[8]], p = cleaned_pres[[8]], a = backgr[[8]], 
                env = predictors)
s9 = prepareSWD(species[[9]], p = cleaned_pres[[9]], a = backgr[[9]], 
                env = predictors)
SWDs = list(s1, s2, s3, s4, s5, s6, s7, s8, s9)

# write object
readr::write_rds(SWDs, "objects/SWDs.rds")
