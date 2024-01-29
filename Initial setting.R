require(demography)
require(ftsa)
require(TSA)

####################################
# Function to read demographic data
####################################

read_demographic_data <- function(state, dataType, type){
  valid_States <- c("ACT", "NSW", "QLD", "WA", "NT", "VIC", "TAS", "SA")
  stopifnot(state %in% valid_States)
  
  # valid_dataType <- c("Births", "Deaths_1x1", "Population", "Exposures_1x1", "Mx_1x1",  "fltper_1x1", "mltper_1x1", "bltper_1x1", "E0per_1x1")
  # stopifnot(dataType %in% valid_dataType)
  
  baseUrl <- "https://aushd.org/assets/txtFiles/humanMortality"
  url <- paste(baseUrl, state, dataType, sep = "/")  # Construct the URL
  
  # Attempt to read data
  tryCatch({
    data <- read.table(url(url), skip = 2, , na.strings = ".", header = TRUE)
    return(data)
  }, error = function(e) {
    cat("Error in reading data:", e$message, "\n")
    return(NULL)
  })
}

valid_States <- c("ACT", "NSW", "QLD", "WA", "NT", "VIC", "TAS", "SA")
# valid_dataType <- c("Births", "Deaths_1x1", "Population.txt", "Exposures_1x1", "Mx_1x1")

results <- list()
dataType <- "Population.txt"
type <- "Population"


# Nested loop to read all data combinations
for (state in valid_States){
  # Read data and store in results list
  results[[paste(state, type, sep = "_")]] <- read_demographic_data(state, dataType, type)
}

female_matrices <- list()
male_matrices <- list()

# Number of rows for each matrix (assuming 111 as in your example)
n_rows <- 111

# Loop through each data frame in the list
for (state in names(results)) {
  # Extract the female and male data and convert them into matrices
  female_matrices[[state]] <- matrix(results[[state]][,3], n_rows)
  male_matrices[[state]] <- matrix(results[[state]][,4], n_rows)
}

level =0.95
method <- "classical"
year = 1971:2021
n_year = length(year)
age = 0:100
n_age = length(age)

for (i in 1:8) {
  female_matrices[[i]][101,] <- colSums(female_matrices[[i]][101:111,])
  female_matrices[[i]] <- female_matrices[[i]][-c(102:111),]
  male_matrices[[i]][101,] <- colSums(male_matrices[[i]][101:111,])
  male_matrices[[i]] <- male_matrices[[i]][-c(102:111),]
  rownames(female_matrices[[i]]) = rownames(male_matrices[[i]]) = age
  colnames(female_matrices[[i]]) = colnames(male_matrices[[i]]) = year
}


au_pop_f <- matrix(0,111,51)
au_pop_m <- matrix(0,111,51)

for (i in 1:8) {
  au_pop_f <- au_pop_f + female_matrices[[i]]
  au_pop_m <- au_pop_m + male_matrices[[i]]
}

au_pop_f[101,] <- colSums(au_pop_f[101:111,])
au_pop_f <- au_pop_f[-c(102:111),]
au_pop_m[101,] <- colSums(au_pop_m[101:111,])
au_pop_m <- au_pop_m[-c(102:111),]
rownames(au_pop_f) = rownames(au_pop_m) = age
colnames(au_pop_f) = colnames(au_pop_m) = year


newborn <- read.csv("ABS_BIRTHS_SUMMARY.csv",header = T)

newborn_2022 <- newborn %>% dplyr:::filter(newborn$TIME_PERIOD..Time.Period==2022)
newborn_2022 <- newborn_2022[,c(2,3,5,6)]
colnames(newborn_2022) <- c("measure","region","time_period","value")
region_map <- list(
  "1: New South Wales" = "NSW",
  "2: Victoria" = "VIC",
  "3: Queensland" = "QLD",
  "4: South Australia" = "SA",
  "5: Western Australia" = "WA",
  "6: Tasmania" = "TAS",
  "7: Northern Territory" = "NT",
  "8: Australian Capital Territory" = "ACT",
  "AUS: Australia" = "AUS"
)

# Replace full names with abbreviations
newborn_2022$region <- sapply(newborn_2022$region, function(x) region_map[[x]])


# Filter the dataset for the valid states and required measures
filtered_data <- newborn_2022 %>%
  dplyr:::filter(region %in% valid_states, measure %in% c("1: Births", "4: Male births"))


# Calculate the male ratio for each state
sex_ratio <- filtered_data %>%
  group_by(region) %>%
  summarise(
    Total_Births = sum(value[measure == "1: Births"]),
    Male_Births = sum(value[measure == "4: Male births"]),
    Male_Ratio = Male_Births / Total_Births
  )

sex_ratio$Female_Ratio <- 1-sex_ratio$Male_Ratio
sex_ratio <- as.data.frame(sex_ratio)

# Convert the region column to a factor with levels in the desired order
sex_ratio$region <- factor(sex_ratio$region, levels = valid_states)

# Arrange the dataframe based on the ordered factor
sex_ratio <- dplyr::arrange(sex_ratio, region)


