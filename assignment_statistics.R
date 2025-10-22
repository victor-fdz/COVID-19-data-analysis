## UAB / MSc in BIOINFORMATICS
## CORE BIOINFORMATICS / STATISTICS 
## ASSIGNMENT / COVID-19 DATA ANALYSIS

## --- LIBRARIES --- ##
library(ggplot2)
library(dplyr)

## --- DATA IMPORT AND CLEANING --- ##
file_path <- "covidworldtable.txt"
df_covid <- read.table(file_path, header = TRUE, sep = "\t")

# Convert numeric-like columns' values to actual numeric
cols2convert <- names(df_covid)[!names(df_covid) %in% c("Country_Other")] # All columns unless country info
df_covid[cols2convert] <- lapply(df_covid[cols2convert], function(x) { # lapply applies a function to all the values of the vector (column)
  round(as.numeric(gsub(",", "", x)), 0) # Delete commas, convert it into numeric and round it to 0 decimals
})

# Fix country names with encoding issues, remove 'World' row and add a new CFR (Case Fatality Rate) column
df_covid$Country_Other[df_covid$order==92] <- "Reunion" # Problem: accent mark in the e
df_covid$Country_Other[df_covid$order==155] <- "Curazao" # Problem: z was not properly encoded
df_covid <- df_covid[!df_covid$Country_Other %in% c("World"),]
df_covid$CFR <- (df_covid$Total_Deaths/df_covid$Total_Cases)*100

# Descriptive statistics
summary(df_covid)



## --- FUNCTION TO SAVE PLOTS WITH THE SAME DIMENSIONS --- ##
save_plot <- function(plot, file_name) {
  ggsave(filename = file_name, plot = plot, width = 7, height = 5, dpi = 300)
}

## --- DATA VISUALIZATION: DISTRIBUTION OF VARIABLES --- ##

# Distribution of cases per million
p1 <- ggplot(df_covid, aes(x = Tot_Cases.1M_pop)) +
  geom_histogram(bins = 50, fill = "#74ade9", color = "white", na.rm = TRUE) +
  labs(title = "Distribution of COVID-19 Cases per Million",
       x = "Cases per Million", y = "Number of Countries") +
  theme_light()
save_plot(p1, "distrib_casesM.png")

# Distribution of deaths per million
p2 <- ggplot(df_covid, aes(x = Deaths.1M_pop)) +
  geom_histogram(bins = 50, fill = "#d73027", color = "white", na.rm = TRUE) +
  labs(title = "Distribution of COVID-19 Deaths per Million",
       x = "Deaths per Million", y = "Number of Countries") +
  theme_light()
save_plot(p2, "distrib_deathsM.png")

# Distribution of tests per million
p3 <- ggplot(df_covid, aes(x = Tests.1M_pop)) +
  geom_histogram(bins = 50, fill = "#1a9850", color = "white", na.rm = TRUE) +
  labs(title = "Distribution of COVID-19 Tests per Million",
       x = "Tests per Million", y = "Number of Countries") +
  theme_light()
save_plot(p3, "distrib_testsM.png")

# Distribution of CFR 
p4 <- ggplot(df_covid, aes(x = CFR)) +
  geom_histogram(bins = 50, fill = "#313695", color = "white", na.rm = TRUE) +
  labs(title = "Distribution of COVID-19 Case Fatality Rate",
       x = "CFR (deaths/cases*100)", y = "Number of Countries") +
  theme_light()
save_plot(p4, "distrib_CFR.png")


## --- TOP COUNTRIES VISUALIZATION --- ##

# Function to generalize plotting
make_barplot <- function(df, var, color, title_y) {
  # Select the top 10 countries for the variable
  top_df <- filter(df, !is.na(.data[[var]])) %>% slice_max(.data[[var]], n = 10) # Filter out NA data. %>% = | (pipe)
  # Plot their values
  ggplot(top_df, aes(x = reorder(Country_Other, .data[[var]]), y = .data[[var]])) +
    geom_col(fill = color, na.rm = TRUE) +
    coord_flip() + # Switch coordinates to have the country name at the Y-axis
    labs(x = "Country", y = title_y, title = paste("Top 10 Regions by", title_y)) +
    theme_light()
}

# Apply the function
p5 <- make_barplot(df_covid,"Tot_Cases.1M_pop","#74add1","Cases per Million")
save_plot(p5, "top10_casesM.png")
p6 <- make_barplot(df_covid,"Deaths.1M_pop","#d73027","Deaths per Million")
save_plot(p6, "top10_deathsM.png")
p7 <- make_barplot(df_covid,"Tests.1M_pop","#1a9850","Tests per Million")
save_plot(p7, "top10_testsM.png")
p8 <- make_barplot(df_covid,"CFR","#313695","CFR")
save_plot(p8, "top10_CFR.png")


## --- CORRELATION ANALYSIS --- ##

# Function to generalize plotting
make_corrplot <- function(df, var1, var2, color, title_x, title_y) {
  n <- sum(complete.cases(df[[var1]], df[[var2]])) # To know the final n used in each case
  ggplot(df, aes(x = .data[[var1]], y = .data[[var2]])) +
    geom_point(alpha = 0.6, color = color) +
    geom_smooth(method = "lm", se = FALSE, color = "#d7191c", linetype = "dashed") + # Plot linear model
    scale_x_log10() + scale_y_log10() + # Logarithmic scale to better observe the whole data range
    labs(
      title = paste("Relation between", var1, "and", var2),
      subtitle = paste("n =", n),
      x = title_x, y = title_y) +
    theme_minimal()
}

# 1. Cases VS Deaths
corr_cases_deaths <- cor.test(df_covid$Tot_Cases.1M_pop, # Variable 1
                              df_covid$Deaths.1M_pop, # Variable 2
                              method = "spearman", # Spearman correlation
                              use = "complete.obs") # Use only rows that have values for both columns
cat("Cases VS Deaths — rho:", corr_cases_deaths$estimate,"p-value:", corr_cases_deaths$p.value, "\n") # Cat enables Python-like printing syntax
p9 <- make_corrplot(df_covid, "Tot_Cases.1M_pop", "Deaths.1M_pop", "#2b83ba", "Cases per Million (log10)", "Deaths per Million (log10)")
save_plot(p9, "casesVsDeaths.png")

# 2. Tests VS Cases
corr_tests_cases <- cor.test(df_covid$Tests.1M_pop, df_covid$Tot_Cases.1M_pop, method = "spearman", use = "complete.obs")
cat("Tests VS Cases — rho:", corr_tests_cases$estimate,"p-value:", corr_tests_cases$p.value, "\n")
p10 <- make_corrplot(df_covid, "Tests.1M_pop", "Tot_Cases.1M_pop", "#1a9850", "Tests per Million (log10)", "Cases per Million (log10)")
save_plot(p10, "testsVsCases.png")

# 3. Tests VS Deaths
corr_tests_deaths <- cor.test(df_covid$Tests.1M_pop, df_covid$Deaths.1M_pop, method = "spearman", use = "complete.obs")
cat("Tests VS Deaths — rho:", corr_tests_deaths$estimate,"p-value:", corr_tests_deaths$p.value, "\n")
p11 <- make_corrplot(df_covid, "Tests.1M_pop", "Deaths.1M_pop", "#d7191c", "Tests per Million (log10)", "Deaths per Million (log10)")
save_plot(p11, "testsVsDeaths.png")

# 4. Population VS Cases per Million
corr_pop_case <- cor.test(df_covid$Population, df_covid$Tot_Cases.1M_pop,method = "spearman", use = "complete.obs")
cat("Population VS Cases per Million — rho:", corr_pop_case$estimate,"p-value:", corr_pop_case$p.value, "\n")
p12 <- make_corrplot(df_covid, "Population", "Tot_Cases.1M_pop", "#74ade9", "Population (log10)", "Cases per Million (log10)")
save_plot(p12, "populationVsCases.png")
