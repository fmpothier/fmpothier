# Spatial Autocorrelation Tutorial
### Frankie Pothier
### October 20th, 2024



## Introduction

Spatial autocorrelation is a critical concept in geographic analysis, focusing on the spatial patterns of attributes to determine whether they are clustered, dispersed, or randomly distributed. This concept relates to Tobler's First Law of Geography, which states that "everything is related to everything else, but closer things are more similar than things further away" (Tobler, 1970). When Tobler's First Law of Geography applies, we observe positive spatial autocorrelation, characterized by clusters of similar values in close proximity. Conversely, when the law does not hold, we may find that nearby locations exhibit more dissimilar values, leading to negative spatial autocorrelation, or a dispersed pattern. Thus, spatial autocorrelation occurs when observation and their neighbouring observation are related to one another (F. Dormann et al., 2007; Tobler, 1970). 

To assess spatial autocorrelation, we will utilize the Moran’s I method. This approach extends Pearson's coefficient into the spatial domain by employing a spatial weights matrix, which captures the relationships between observations based on their locations (Westerholt, 2023). Specifically, Moran’s I tests the degree of similarity or difference between a given value and its neighboring values.

The power of spatial autocorrelation analysis becomes particularly evident when examining census data. The census serves as a crucial resource, offering comprehensive insights into the population of Canada. Conducted every five years, the census provides a wealth of information about demographic factors such as age, sex, household structure, languages spoken, and income. Importantly, this information is anonymized, ensuring that it cannot be linked to individual households (Story - Understanding Census Data, n.d.).

Census data is especially well-suited for testing spatial autocorrelation because it has a diverse array of attributes tied to geographic locations. By examining these attributes, we can see if certain demographic or socioeconomic factors show different spatial patterns, which helps us understand how they relate to each other in different areas.

In this tutorial, we will learn how to analyze spatial autocorrelation using R and census data. We will guide you through the process of setting up your R workspace, loading the necessary libraries, cleaning the data, and running and analyzing both Global and Local Moran’s I for the city of Kamloops. Our focus will be on examining the variables of median total income and the percentage of individuals with knowledge of the French language.



### Libraries in R 

In R, libraries are collections of functions designed to perform specific tasks (“Packages and libraries,” 2020). There are two types of libraries: base libraries, which come pre-installed with R, and ‘add-on’ libraries, which can be downloaded separately (“Packages and libraries,” 2020). The base libraries enable R to function properly, while ‘add-on’ libraries extend R’s capabilities by providing additional, pre-defined functions. These allow users to avoid writing every function from scratch, greatly speeding up their workflow. Since R is an open-source language, anyone can create a library, leading to thousands of available libraries. To use one of these ‘add-on’ libraries, it first needs to be downloaded to your computer using a command like ‘install.packages("library_name")’, which installs all the necessary components. After that, you need to load the library into your program to make its functions available using the command library(“library_name”).

For this analysis of spatial autocorrelation using census data, we’ll be using several ‘add-on’ libraries: Knitr, tmap, spdep, e1071, sf, st, and raster. The Knitr library is useful for generating PDFs, graphs, tables, and other output. One common function we’ll use from Knitr is kable, which helps format tables to display our results. The tmap library is ideal for creating maps, with functions like tm_shape, tm_polygons, and tmap_arrange that allow us to plot spatial data, create polygon layers, and organize multiple map layers. The spdep and e1071 libraries are used for statistical analysis. Spdep focuses on spatial statistics, while e1071 covers general statistics and machine learning algorithms. The sf and st libraries are essential for working with shapefiles and other spatial data formats, while the raster library is specialized for handling raster data. Throughout this tutorial, we’ll explore how to use each of these libraries to enhance our analysis. Below shows how to add these libraries.


```{r Libraries, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}

#Install packages if not already installed:
#install.packages("knitr")
#install.packages("tmap")
#install.packages("spdep")
#install.packages("e1071")
#install.packages("sf")
#install.packages("st")
#install.packages("raster")

#Load in libraries:
library("knitr")
library("tmap")
library("spdep")
library("e1071")
library("sf")
library("st")
library("raster")

```

### Working Directory and data files 

Just like making sure the necessary libraries are properly set up, configuring your working directory is essential. The working directory is the folder where all your source code and data are saved (Phillips, 2017), and it’s where all the necessary files for your project should be stored. For the R program to find these files the path to them must be set.

In this tutorial, we’ll be working with two data files from our working directory: a CSV file and a shapefile. The CSV file contains the census data needed for this analysis, while the shapefile holds the spatial data. A CSV (Comma-Separated Values) file stores tabular data where each column is separated by a comma. A shapefile, on the other hand, is a type of file used to store geographic features, such as polygons, lines, or points (ArcMap, n.d.).

To load our data into R, we need to specify the file paths. First, we’ll load the CSV file into a variable called csv using the ‘read.csv’ function from R’s base libraries, which reads the file and stores it as a data frame in our program. After that, we’ll load the shapefile using the ‘st_read’ function from the ‘sf’ library to bring the spatial data into the program.


```{r Read in data, echo=TRUE, eval=TRUE, warning=FALSE}

#From working directory upload the csv data set
csv <- read.csv("~/Downloads/classes/2024/Fall/GEOG 418/assignments/ass3/Assignment3_data/ucgsJQnBVLvP_data.csv") 

#From working directory upload the shp file 
shp <- st_read("~/Downloads/classes/2024/Fall/GEOG 418/assignments/ass3/Assignment3_data/lda_000b16a_e.shp") 

```

### Cleaning the data  

Before moving on to the analysis, it’s important to clean the data to ensure it’s ready for use. Data cleaning is crucial because it helps us avoid any errors or inconsistencies that could affect the results. For example, we want to make sure there are no problematic values that could disrupt our analysis.

We’ll start by renaming the columns in the CSV file so that it’s clear what each column represents. To do this, we’ll create a list (or vector in R) of new column names, making sure they’re in the correct order. Once we have the vector, we can use the 'colnames' function to apply the new names to the data frame. After renaming the columns, we’ll remove any unwanted rows. In this case, we’ll eliminate any rows where the ID has fewer than eight digits.

For this analysis, we’ll focus on a single city—Kamloops. We’ll create a subset of the data that includes only the rows relevant to Kamloops. Once we’ve created this subset, we’ll calculate the percentage of people who speak French in each census tract. This percentage will be added as a new column called PercFrench, which will allow us to map and analyze the percentage of French speakers across different census tracts.


```{r Clean data, echo=TRUE, eval=TRUE, warning=FALSE}

#Create new colemn names for the csv file
cols <- c("GEO UID", "Province code", "Province name", "CD code",
          "CD name", "DA name", "Population", "Land area", 
          "Median total income", "Income Sample Size", "French Knowledge", 
          "Language Sample Size")

#Apply those names to dataframe
colnames(csv) <- cols

#Add column to count number of ID charactors
csv$len <- nchar(csv$`GEO UID`)

#Remove IDs with less than 8 numbers
csv_clean <- subset(csv, csv$len == 8)

#Merge spatial and aspatial data
census_DAs <- merge(shp, csv_clean, 
                    by.x = "DAUID", 
                    by.y = "GEO UID", 
                    all.x = TRUE)

#Subset for Kamloops (or the census needed)
Municp <- subset(census_DAs, census_DAs$CMANAME == "Kamloops")

#Calucate the precent of people that speck french per census track
Municp$PercFrench <- (Municp$`French Knowledge`/Municp$`Language Sample Size`)*100

```

The final step in cleaning the data is to remove any rows that contain missing values. If left in, these missing values can distort the results of the analysis or even prevent it from working properly, so it is crucial to remove them. To do this, we can eliminate any rows with NA or zero values in the columns of interest. Since we are testing for spatial autocorrelation of median total income and knowledge of French in the city of Kamloops, we need to remove any rows with NA or zero values in these specific columns.

The code below demonstrates how to achieve this using the 'which' and 'is.na; functions. The 'is.na' function checks whether a value is NA, and we use the negation operator (!) in front of it because we want to select rows that do not have NA values. The 'which' function then creates a new data frame consisting of only the rows that do not have NA values in the columns of interest.


```{r NA Remove, echo=TRUE, eval=TRUE, warning=FALSE}

#Remove Income NA
Income_noNA <- Municp[which(!is.na(Municp$`Median total income`)),]

#Remove French NA
French_noNA <- Municp[which(!is.na(Municp$`PercFrench`)),]

```

### Descriptive Statistics 

To start our analysis, we will compute some descriptive statistics on the variables of interest: median total income and the percentage of people with French language knowledge. This will help us understand key aspects of the data, such as its distribution.

First, we will compute the mean and standard deviation using built-in R functions. Then, we will calculate the skewness of the data using the ‘skewness’ function from the ‘e1071’ library. Once all the descriptive statistics are calculated, we can add them to a data frame to ensure the results are properly structured. This will allow us to create a table using the ‘kable’ function to view the results (see Table 1).


```{r DescriptiveStats, echo=TRUE, eval=TRUE, warning=FALSE}

#Calculate descriptive stats for Income
meanIncome <- mean(Income_noNA$`Median total income`)
stdevIncome <- sd(Income_noNA$`Median total income`)
skewIncome <- skewness(Income_noNA$`Median total income`)

#Calculate descriptive stats for French
meanFrench <- mean(French_noNA$`PercFrench`)
stdevFrench <- sd(French_noNA$`PercFrench`)
skewFrench <- skewness(French_noNA$`PercFrench`)

#Create dataframe for display in table
data <- data.frame(Variable = c("Income", "French Language"),
                   Mean = c(round(meanIncome,2), round(meanFrench,2)),
                   StandardDeviation = c(round(stdevIncome,2), round(stdevFrench,2)),
                   Skewness = c(round(skewIncome,2), round(skewFrench,2)))

#Produce table
kable(data, caption = paste0("Descriptive statistics for median total income and 
                             the percentage of people with French language 
                             knowledge in the ", 2016, " Census Tract of 
                             Kamloops (Government of Canada, 2021)"))

```

### Creating Census Dissemination Areas Maps 

Next, the census dissemination areas in Kamloops will be mapped to illustrate the median total income and the percentage of people with knowledge of the French language per census tract. This mapping will help us understand how this data are spatially distributed across Kamloops.

To start, a map showing the median total income for each census tract will be created using the ‘tmap’ library. The ‘tm_shape’ function will specify the spatial data being plotted, which is the ‘Income_noNA’ spatial object. Following that, the ‘tm_polygons’ function will define various aspects of the map, including the specific column from the spatial object being mapped, the title, style, color palette, border line weight, and the color for NA values. For this analysis, natural breaks (Jenks) will be used to categorize the data, identifying significant gaps between data values. The color palette will feature shades of red, with light red representing low values and dark red indicating high values.

Once these attributes are set, the ‘tm_layout’ function will be applied to position the legend on the map. A similar process will be followed to create a map showing the percentage of people with knowledge of French.

After both maps are prepared, the ‘tmap_arrange’ function will arrange them side by side for display. This function allows specification of the number of columns and rows for the layout. Since the goal is to have the maps adjacent, the number of columns will be set to 2 and the number of rows to 1, resulting in a single row containing both maps.


```{r StudyArea, echo=TRUE, eval=TRUE, warning=FALSE, fig.cap="The Maps above show the Kamloops census dissemination areas. With the map on the left showing median total income in red, with light red being low values and dark red being high values. The Map on the right is showing the percentage of people with knowledge of French."}

#Choose a colour palette

#uncommened below line to help find colour palettes
#tmaptools::palette_explorer()

#Map median Income
map_Income <- tm_shape(Income_noNA) + 
  tm_polygons(col = "Median total income", 
              title = "Median total income", 
              style = "jenks", 
              palette = "Reds", n = 6,
              border.alpha = 0,
              colorNA = "grey") +
  tm_layout(legend.position = c("RIGHT", "TOP"))

#Map French Knowledge
map_French <- tm_shape(French_noNA) + 
  tm_polygons(col = "PercFrench", 
              title = "Percentage with \n French Knowledge", 
              style = "jenks", 
              palette = "Reds", n = 6,
              border.alpha = 0,
              colorNA = "grey") +
  tm_layout(legend.position = c("RIGHT", "TOP"))

#Print maps side by side
tmap_arrange(map_Income, map_French, ncol = 2, nrow = 1)

```

## Neighbourhood matrix 

The neighborhood matrix is used to describe which polygons share a common border, allowing spatial autocorrelation to be examined in the study area (Moraga, 2023). In the neighborhood matrix, weights are assigned to indicate whether areas are neighbors (Moraga, 2023).

In R, the spdep library can be used to create a neighborhood matrix and perform spatial autocorrelation analysis (Bivand, 2002). The function poly2nb enables users to create a neighborhood matrix and define what constitutes a neighbor. There are two main ways to define a neighbor: the first is the queen type, which states that a neighbor needs only one shared boundary point (Moraga, 2023). The second way is the rook type, which states that a neighbor needs two shared boundary points (Moraga, 2023).

By using the poly2nb function, both kinds of neighborhood matrices can be easily created. If the queen variable is set to TRUE, a neighborhood matrix based on the queen definition of a neighbor will be created; if this variable is not set, it will default to the rook definition of a neighbor (Bivand, 2002). The code below demonstrates how to create rook and queen-defined neighborhood matrices using the median total income and percentage of people with knowledge of the French language per census tract in Kamloops.


```{r Neighbours, echo=TRUE, eval=TRUE, warning=FALSE}

#Income Neighbours - Queens weight
Income.nb <- poly2nb(Income_noNA)
Income.rook <- nb2lines(Income.nb, coords=st_coordinates(st_centroid(Income_noNA)))
crs(Income.rook) <- crs(Income_noNA)

#Income Neighbours - Rooks weight
Income.nb2 <- poly2nb(Income_noNA, queen = FALSE)
Income.queen <- nb2lines(Income.nb2, coords=st_coordinates(st_centroid(Income_noNA)))
crs(Income.queen) <- crs(Income_noNA)

#French Neighbours - Queens weight
French.nb <- poly2nb(French_noNA)
French.queen <- nb2lines(French.nb, coords=st_coordinates(st_centroid(French_noNA)))
crs(French.queen) <- crs(French_noNA)

#French Neighbours - Rooks weight
French.nb2 <- poly2nb(French_noNA, queen = FALSE)
French.rook <- nb2lines(French.nb2, coords=st_coordinates(st_centroid(French_noNA)))
crs(French.rook) <- crs(French_noNA)

```

The maps below illustrate the links between neighbors and how these links differ based on the type of neighborhood matrix used: queen or rook. The map in green displays the queen neighborhood links for the census tracts in Kamloops, while the blue map shows the rook neighborhood links (Figure 2). The map on the far right highlights the differences between the rook and queen neighborhood links (Figure 2). As seen in these maps, the neighborhood links differ only slightly (Figure 2).

The ‘tmap’ library was utilized to create the maps presented in Figure 2. The first step involved mapping the census tracts in Kamloops. The ‘tm_shape’ function was used to specify the spatial object being mapped, followed by the ‘tm_borders’ function to define the color of the census tract borders. Once this was complete, the ‘tm_shape’ function was used again to specify that the rook/queen neighborhood links were being mapped, with the ‘tm_lines’ function defining the color and line weight of the links. The addition signs are used to combine all components of the map into one cohesive map. Once each map is created, the ‘tmap_arrange’ function is used to arrange and print the maps.

```{r Neighboursmap, echo=TRUE, eval=TRUE, warning=FALSE, fig.cap="The maps below show the neighbor links created using a rook weight matrix (middle), a queen weight matrix (left), and both matrices together on the right for the Kamloops census dissemination areas using the median total income."}

#Make queens map
IncomeQueen <- tm_shape(Income_noNA) + tm_borders(col='lightgrey') + 
  tm_shape(Income.queen) + tm_lines(col='green')

#Make rooks map
IncomeRook <- tm_shape(Income_noNA) + tm_borders(col='lightgrey') + 
  tm_shape(Income.rook) + tm_lines(col='blue', lwd = 2)

#Make combined map
IncomeBoth <- tm_shape(Income_noNA) + tm_borders(col='lightgrey') + 
  tm_shape(Income.rook) + tm_lines(col='blue', lwd = 2) +
  tm_shape(Income.queen) + tm_lines(col='green', lwd = 1)

#Print maps in a three pane figure
tmap_arrange(IncomeQueen, IncomeRook, IncomeBoth, ncol = 3, nrow = 1)

```

The code below uses the ‘nb2listw’ function from the ‘spdep’ library to convert a neighborhood object into a spatial weights object (Bivand, 2002). The nb2listw function takes the neighborhood object, the weighting style, and the ‘zero.policy’ variable as inputs (Bivand, n.d.).
The weighting style defines how the weight matrix is created. There are six different styles that can be specified with the ‘nb2listw’ function: “W”, “B”, “C”, “U”, “minmax”, and “S” (Bivand, n.d.). The “W”, “C”, and “B” styles are the most commonly used.

The “W” style uses a row-standardized weighting scheme, meaning that the sum of the weights for all neighbors adds up to 1, and each neighbor receives an equal share. For instance, if there are four neighbors, each would receive a weight of 0.25. Since this is the style used in this tutorial, the style argument is set to “W”.

The “B” style represents binary weighting, where each neighbor gets a weight of 1 and non-neighbors get a weight of 0. The “C” style uses a globally standardized weighting scheme, where all neighbors across the study area receive the same weight.

The ‘zero.policy’ variable determines how the function handles polygons with no neighbors (zero-length weight vectors). If ‘zero.policy’ is set to FALSE, the function will return an error when encountering a polygon with no neighbors. If it is set to TRUE, the algorithm will allow the inclusion of polygons with no neighbors, assigning them a weight of zero in the weight matrix (Bivand, n.d.).

To learn more about how to use the ‘nb2listw’ function, you can find a detailed description [here](https://r-spatial.github.io/spdep/reference/nb2listw.html).

The code below converts the queen neighborhood object for median total income and the percentage of people with knowledge of the French language in Kamloops (per census tract) into spatial weights objects using the “W” style and with ‘zero.policy’ set to TRUE.


```{r Final weights, echo=TRUE, eval=TRUE, warning=FALSE}

#Create Income weights matrix
Income.lw <- nb2listw(Income.nb, zero.policy = TRUE, style = "W")

#Create French weights matrix
French.lw <- nb2listw(French.nb, zero.policy = TRUE, style = "W")

head(Income.lw[["weights"]])[c(1:3)]

```


## Global Moran’s I

Global Moran’s I is a spatial statistical method used to test for spatial autocorrelation in a given study area. It examines how similar or different values are compared to their neighboring values. Global Moran’s I does this by measuring how much a given location's value deviates from the mean and comparing that to how much its neighbors' values deviate from the mean. It extends the Pearson correlation coefficient to the spatial domain by incorporating spatial weights, which describe the relationships between neighbors (Westerholt, 2023). The formula below is used to calculate Global Moran’s I.

$$
I = \frac{\sum_{i=1}^n\sum_{j=1}^nW_{i,j}(x_i - \bar{x})(x_j - \bar{x})}{(\sum_{i=1}^n\sum_{j=1}^nW_{i,j})\sum_{i=1}^n(x_i - \bar{x})^2}
$$

The formula above shows that for each value $x$, $x_i$ represents the value at location $i$, and $x_j$ is the value at a neighboring location. Moran’s I compares each value $x_i$ to the mean of the dataset and does the same for its neighbor $x_j$, capturing how much each value deviates from the overall average. These deviations are multiplied by the spatial weight $W_{i,j}$, which represents the strength of the spatial relationship between locations $i$ and $j$. The spatial weights account for how close or connected the locations are. The denominator standardizes the output, ensuring that values of Moran’s I are consistent. Higher values of Moran's I indicate strong positive spatial autocorrelation (similar values cluster together), while lower values suggest strong negative spatial autocorrelation (dissimilar values cluster together).


```{r Global Morans I, echo=TRUE, eval=TRUE, warning=FALSE}

#Calculate Global Moran's I for Income
miIncome <- moran.test(Income_noNA$`Median total income`, Income.lw, zero.policy = TRUE)

#Extract Global Moran's I results for Income
mIIncome <- miIncome$estimate[[1]]
eIIncome <- miIncome$estimate[[2]]
varIncome <- miIncome$estimate[[3]]

#Calculate Global Moran's I for French
miFrench <- moran.test(French_noNA$PercFrench, French.lw, zero.policy = TRUE)

#Extract Global Moran's I results for French
mIFrench <- miFrench$estimate[[1]]
eIFrench <- miFrench$estimate[[2]]
varFrench <- miFrench$estimate[[3]]

```

The results for the Global Moran’s I of median total income and the percentage of people with knowledge of the French language per census tract in Kamloops are presented in Table 2. The Global Moran’s I for median total income is 0.6127 (see Table 2), indicating a strong positive spatial autocorrelation. This suggests that areas with similar income levels tend to cluster together in Kamloops, meaning high-income areas are often located near other high-income areas, while low-income areas are situated adjacent to other low-income regions.

In contrast, the Global Moran’s I for the percentage of people with knowledge of the French language is 0.3716 (see Table 2). This also reflects a positive spatial autocorrelation, although it is weaker than that observed for median total income.

The expected Moran’s I provides a baseline value that would occur in the absence of spatial autocorrelation, reflecting a random distribution of attributes (see Table 2). Since the Global Moran’s I values for both median total income and the percentage of people with knowledge of the French language exceed the expected Moran’s I, this reinforces the presence of positive spatial autocorrelation for both variables within the census tracts of Kamloops.


```{r Morans I table, echo=FALSE, eval=TRUE, warning=FALSE}

# Create a data frame with the extracted results
moran_results <- data.frame(
  Statistic = c("Global Moran's I", "Expected Moran's I", "Variance"),
  Income = c(mIIncome, eIIncome, varIncome),
  French = c(mIFrench, eIFrench, varFrench)
)

# Display the table
kable(moran_results, caption = "Global Moran's I Results for Median Income and 
      French Knowledge")

```

The code below creates a function called moran.range that finds the maximum and minimum values of Global Moran’s I. This function is then used to determine the minimum and maximum Global Moran’s I for the Median Total Income in the census tracts of Kamloops, resulting in a maximum of 1.0537 and a minimum of -0.7151. Analyzing the range of Moran’s I across a given study area can help provide a deeper understanding of the spatial autocorrelation in the area.

```{r Global Morans Range, echo=TRUE, eval=TRUE, warning=FALSE}

#Function to calculate the range of global Moran's I
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}

#Calculate the range for the Income variable
range <- moran.range(Income.lw)
minRange <- range[1]
maxRange <- range[2]

```

The results from the Global Moran’s I analysis reveal that both median total income and the percentage of people with knowledge of the French language exhibit positive spatial autocorrelation. Notably, the strength of the positive spatial autocorrelation for median total income is greater than that for the percentage of people with knowledge of French.

To assess whether this positive spatial autocorrelation is significantly different from random, we can employ a z-test. The first step in conducting a z-test is to establish the null and alternative hypotheses. In this analysis, the null hypothesis states that both median total income and the percentage of people with knowledge of the French language in Kamloops are spatially random, indicating the absence of spatial autocorrelation. Conversely, the alternative hypothesis asserts that spatial autocorrelation exists for both variables in Kamloops.

For this study, we set the significance level at 95%, which corresponds to an $\alpha$ value of 0.05. Consequently, if the z-score falls below -1.96 or exceeds +1.96, we can reject the null hypothesis and conclude that there is statistically significant spatial autocorrelation present.

The z-score for spatial autocorrelation for both median total income and the percentage of people with knowledge of the French language in Kamloops is calculated in the code below. This code subtracts the expected Moran's I from the global Moran's I and then divides it by the square root of the variance. The results can be found in Table 3.


```{r Global Morans ZScore, echo=TRUE, eval=TRUE, warning=FALSE}

#Calculate z-test for Income
zIncome <- (mIIncome - eIIncome) / (sqrt(varIncome))

#Calculate z-test for French
zFrench <- (mIFrench - eIFrench) / (sqrt(varFrench))

```

```{r Z-score table, echo=FALSE, eval=TRUE, warning=FALSE}

# Create a data frame with the extracted results
Z_results <- data.frame(
  Variables = c("Income", "French Knowledge"),
  Z_scores = c(zIncome, zFrench))

# Display the table
kable(Z_results, caption = "Z-score's for Median Income and 
      French Knowledge")

```

The resulting z-scores for spatial autocorrelation are 12.5769 for median total income and 7.7564 for the percentage of people with knowledge of the French language in Kamloops (see Table 3). These values demonstrate that both median total income and the percentage of people with knowledge of the French language exhibit statistically significant positive spatial autocorrelation.

## Local spatial autocorrelation

Local spatial autocorrelation differs from global spatial autocorrelation as it focuses on each observation and its neighboring values, rather than providing a single value to summarize the spatial autocorrelation across the entire study area (Rey et al., 2023). This allows us to explore the spatial distribution of the attributes of individual observations more deeply (Rey et al., 2023). The local spatial autocorrelation statistic also helps identify areas with strong spatial patterns, such as clustering or outliers that differ from their neighbors (Bivand et al., 2009). Calculating local spatial autocorrelation can help determine where clustering occurs or where values are dispersed.

Local Moran’s I is used to calculate local spatial autocorrelation. Like Global Moran’s I, it measures spatial autocorrelation, but the key difference is that Local Moran’s I provides a specific measure for each individual location, rather than summarizing the entire study area with a single value. The formula below is used to calculate Local Moran’s I for a given area. It does this by finding the difference between an observation $i$ and the mean, and then multiplying it by the sum of the differences between its neighbors and the mean.


$$
I_i = \frac{x_i - \bar{x}}{S_i^2}\sum{_{j=1}^n}W_{i,j}(x_j - \bar{x})\space \space where \space \space S_i^2 = \frac{\sum_{i=1}^n (x_i - \bar{x})^2}{n-1} 
$$

To compute Local Moran’s I in R, the ‘localmoran’ function from the ‘spdep’ library (Bivand, 2002) is used. This function requires two inputs: a numeric vector of observations and a spatial weights object, both of which must have the same length. In this tutorial, the numeric vectors contain values for median total income and the percentage of people with knowledge of the French language per census tract in Kamloops, while the spatial weights object represents the matrix created earlier. After running the ‘localmoran’ function, we can extract the results from its output for further analysis.

```{r Local Morans I, echo=TRUE, eval=TRUE, warning=FALSE}

#Calculate LISA test for Income
lisa.testIncome <- localmoran(Income_noNA$`Median total income`, Income.lw)

#Extract LISA test results for Income
Income_noNA$Ii <- lisa.testIncome[,1]
Income_noNA$E.Ii<- lisa.testIncome[,2]
Income_noNA$Var.Ii<- lisa.testIncome[,3]
Income_noNA$Z.Ii<- lisa.testIncome[,4]
Income_noNA$P<- lisa.testIncome[,5]

#Calculate LISA test for Income
lisa.testFrench <- localmoran(French_noNA$PercFrench, French.lw)

#Extract LISA test results for Income
French_noNA$Ii <- lisa.testFrench [,1]
French_noNA$E.Ii<- lisa.testFrench [,2]
French_noNA$Var.Ii<- lisa.testFrench [,3]
French_noNA$Z.Ii<- lisa.testFrench [,4]
French_noNA$P<- lisa.testFrench [,5]

```

The code below maps the results of the local Moran’s I to visualize the data. We utilize the 'tmap' library, similar to the approach used for the maps above.

```{r MappingLocalMoransI, echo=TRUE, eval=TRUE, warning=FALSE, fig.cap="Local spatial autocorrelation z-scores for Kamloops census dissemination areas: median total income (left) and percentage of people with knowledge of French (right)."}

#Map LISA z-scores for Income
map_LISA_Income <- tm_shape(Income_noNA) +
  tm_polygons(col = "Z.Ii",
              title = "Local Moran's I Z-Scores: median total income",
              style = "fixed",
              border.alpha = 0.1,
              midpoint = NA,
              colorNA = NULL,
              breaks = c(min(Income_noNA$Z.Ii),-1.96,1.96,max(Income_noNA$Z.Ii)),
              palette = "-RdBu", n = 3)+
  tm_compass(position=c("left", "top"))+
  tm_scale_bar(position=c("left", "bottom"))+
  tm_legend(position = c("right", "top"))

#Map LISA z-scores for French
map_LISA_French <- tm_shape(French_noNA) +
  tm_polygons(col = "Z.Ii",
              title = "Local Moran's I Z-Scores: French language",
              style = "fixed",
              border.alpha = 0.1,
              midpoint = NA,
              colorNA = NULL,
              breaks = c(min(French_noNA$Z.Ii),-1.96,1.96,max(French_noNA$Z.Ii)),
              palette = "-RdBu", n = 3)+
  tm_compass(position=c("left", "top"))+
  tm_scale_bar(position=c("left", "bottom"))+
  tm_legend(position = c("right", "top"))

#Plot maps in a 2 pane figure
tmap_arrange(map_LISA_Income, map_LISA_French, ncol = 2, nrow = 1)

```

The map on the left of Figure 3 displays the z-scores from the local Moran’s I test for median total income per census track in Kamloops. This map highlights a few census tracks with statistically significant spatial autocorrelation, marked in red, while most tracks do not exhibit significant autocorrelation, shown in gray. This indicates that these few areas are responsible for the global Moran’s I demonstrating a significant positive spatial autocorrelation.

On the right side of Figure 3, we see the results for the z-scores from the local Moran’s I test for the percentage of people with knowledge of the French language per census track in Kamloops. This map reveals several areas with significant spatial autocorrelation, indicated in red and blue.

While these maps effectively visualize the results of the local Moran’s I, scatterplots of local spatial autocorrelation can provide even deeper insights. The code below is utilized to create scatterplots for local spatial autocorrelation for both median total income and the percentage of people with knowledge of the French language per census track in Kamloops. To accomplish this, we employ the ‘moran.plot’ function from the ‘spdep’ library (Bivand, 2002), which facilitates the creation of scatter plots for spatial data.

```{r MoransIScatter, echo=TRUE, eval=TRUE, warning=FALSE, fig.cap= "Moran's I scatter plot for median total income in Kamloops, 2016."}

#Create Moran's I scatter plot for Income
moran.plot(Income_noNA$`Median total income`, Income.lw, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Median Total Income ($)", 
           ylab="Spatially Lagged Median Total Income ($)", quiet=NULL)

```


```{r MoransIScatter2, echo=TRUE, eval=TRUE, warning=FALSE, fig.cap= "Moran's I scatter plot for percentage of people with knowledge of French in Kamloops, 2016."}

#Create Moran's I scatter plot for French
moran.plot(French_noNA$PercFrench, French.lw, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Respondants with knowledge of French (%)", 
           ylab="Spatially Lagged knowledge of French (%)", quiet=NULL)

```

When examining the Moran’s I scatter plot for median total income in Kamloops in 2016, we observe statistically significant values, represented by diamonds, clustered in Quadrant 1 and Quadrant 3 of the graph (see figure 4). Points in Quadrant 1 indicate areas where both the observed value and its neighbors are high, suggesting positive spatial autocorrelation, meaning that high values tend to cluster together. Furthermore, points in Quadrant 3 represent areas where both the observed value and their neighboring values are low, which also indicates positive spatial autocorrelation. Thus, the analysis of the Moran’s I scatter plot for median total income in Kamloops in 2016 suggests that positive spatial autocorrelation is present.

Similarly, when examining the Moran’s I scatter plot for the percentage of people with knowledge of the French language in Kamloops in 2016, we observe statistically significant values primarily clustering in Quadrant 1 (see figure 5). This clustering indicates that both the observed values and their neighbors possess knowledge of the French language, suggesting positive spatial autocorrelation. Therefore, the analysis of the Moran’s I scatter plot for the percentage of people with knowledge of the French language in Kamloops in 2016 also suggests that positive spatial autocorrelation is present.

## Summary

In this tutorial, we explored how to use R and its libraries to compute spatial autocorrelation for median total income and the percentage of individuals with knowledge of the French language in Kamloops in 2016. Our analysis revealed statistically significant positive spatial autocorrelation for both median total income and the percentage of people fluent in French. The Local Moran's I statistic demonstrated that areas with low median total income are clustered together, as are areas with high total income. Similarly, the Local Moran's I for the percentage of French speakers in each census tract showed significant positive spatial autocorrelation, indicating that individuals proficient in French tend to live near others with the same linguistic skills. This finding suggests the presence of a small French-speaking community forming in Kamloops.

Future research could delve deeper into the underlying factors contributing to these patterns and apply these methods to subsequent census data to examine whether this relationship has changed over time.

## References

Bivand, R. (n.d.). Spatial weights for neighbours lists. Github.Io. Retrieved October 20, 2024, from https://r-spatial.github.io/spdep/reference/nb2listw.html

Bivand, R. (2002). spdep: Spatial Dependence: Weighting Schemes, Statistics [Data set]. In CRAN: Contributed Packages. The R Foundation.

Bivand, R., Müller, W. G., & Reder, M. (2009). Power calculations for global and local Moran’s. Computational Statistics & Data Analysis, 53(8), 2859–2872. https://doi.org/10.1016/j.csda.2008.07.021

F. Dormann, C., M. McPherson, J., B. Araújo, M., Bivand, R., Bolliger, J., Carl, G., G. Davies, R., Hirzel, A., Jetz, W., Daniel Kissling, W., Kühn, I., Ohlemüller, R., R. Peres-Neto, P., Reineking, B., Schröder, B., M. Schurr, F., & Wilson, R. (2007). Methods to account for spatial autocorrelation in the analysis of species distributional data: a review. Ecography, 30(5), 609–628. https://doi.org/10.1111/j.2007.0906-7590.05171.x

Government of Canada, & Canada, S. (2001, January 15). Census of Population. Statcan.Gc.Ca. https://www12.statcan.gc.ca/census-recensement/index-eng.cfm

Moraga, P. (2023). Spatial statistics for data science: Theory and practice with R. Chapman & Hall/CRC.

Packages and libraries. (2020, April 23). Introduction to R. https://hbctraining.github.io/Intro-to-R-flipped/lessons/04_introR_packages.html

Phillips, N. (2017). YaRrr! The Pirate’s Guide to R. APS Observer, 30. https://bookdown.org/ndphillips/YaRrr/the-working-directory.html

Rey, S., Arribas-Bel, D., & Wolf, L. J. (2023). Geographic data science with python. Chapman & Hall/CRC.

Story - understanding Census data. (n.d.). Peelregion.Ca. Retrieved October 19, 2024, from https://data.peelregion.ca/pages/story-understanding-census-data

Tobler, W. R. (1970). A computer movie simulating urban growth in the Detroit region. Economic Geography, 46, 234. https://doi.org/10.2307/143141

Westerholt, R. (2023). A simulation study to explore inference about global Moran’s I with random spatial indexes. Geographical Analysis, 55(4), 621–650. https://doi.org/10.1111/gean.12349

What is a shapefile?—ArcMap. (n.d.). Arcgis.com. Retrieved October 20, 2024, from https://desktop.arcgis.com/en/arcmap/latest/manage-data/shapefiles/what-is-a-shapefile.htm
