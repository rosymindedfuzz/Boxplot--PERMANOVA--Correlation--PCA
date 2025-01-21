#Geochemical Database of City X in Aotearoa New Zealand (by M. Banta and Y. Wang)
#Our dataset consists of elemental concentrations (including trace metals) obtained from Depths A and B in a certain city (City X) in Aotearoa New Zealand. 
#These concentrations correspond to certain land covers in the city.

#Our objectives are:
#1. To create a boxplot of our element(s) of interest.
#2. To know if the land cover has an effect on metal concentrations, by performing PERMANOVA.
#3. To visualize the relationship among our variables for Depth A, by performing Principal Component Analysis.


#Part 1 - Setup

#Clear existing data
rm(list=ls())

#Specify the path name
setwd("C:/Users/rosymindedfuzz")

#Install the necessary packages
install.packages("rio") #This package gives us a nice code that will combine all the sheets in our Excel spreadsheet.
install.packages("dplyr") # This package uses verbs like “select” and “filter” to manipulate our data
install.packages("reshape2") #This package allows us to reshape our data from wide format to long format and vice versa                
install.packages("ggplot2")  #This package lets us visualize our data
install.packages("vegan") #This package allows us to do multivariate analysis such as PERMANOVA.
install.packages("factoextra") #PCA analysis
install.packages("corrplot") #Correlation plot
install.packages("ggfortify")  #PCA plot    

#Load the “rio” library
library(rio)

#Import our spreadsheet
geochem <- import_list("Geochemical Database of City X.xlsx" , rbind=TRUE)

#Inspect data
head(geochem)
View (geochem)

#Part 2 - Data Cleaning

#Convert our dataset into a data frame
geochem.df = data.frame(geochem)

#Inspect data frame
str(geochem.df)

#Our element of interest for this section is Cd (cadmium), because aside from being toxic to human health even in low concentrations, its presence is due to anthropogenic sources.

#Check normality
hist(geochem.df$Cd_ppm)

#Our data is skewed to the right (positively skewed).

#Attempt to transform

#Square root transformation
geochem.t<-sqrt(geochem.df$Cd_ppm)
hist(geochem.t)

#Log transformation
geochem.t<-log(geochem.df$Cd_ppm)
hist(geochem.t)

#Inverse transformation
geochem.t<-1/(geochem.df$Cd_ppm)
hist(geochem.t)

#It seems like none of our transformations can produce plots that show normal distribution.
#We'll just use geochem.df and assume that our data does not produce a normal distribution.

#Part 3 - Summary Statistics

#For our next exercise, we're only interested in metal concentrations that were obtained using ICP-MS

#Load the "dplyr" library
library(dplyr)

#Select range of columns that we want to keep
geochem.df.1 <- geochem.df%>% select(Depth, (Ag_ppm:Zr_ppm))
View(geochem.df.1)

#Produce summary statistics
summary(geochem.df.1)

#Part 4 - Creating a boxplot of elemental concentrations

#Load "ggplot2" library
library(ggplot2)

#Reshape data. Currently, our data is in the wide format. We want to reshape it into a long format.
library(reshape2)      	 
geochem.df.1.long <- melt(geochem.df.1, id = "Depth")                 	 
View(geochem.df.1.long)

#Make separate boxplots for each metal, but make them all appear in one page
p<-ggplot(geochem.df.1.long, aes(x = variable, y = value, color = Depth)) +  
  geom_boxplot()
p + facet_wrap( ~ variable, scales="free")

#Plot element of interest: Cd (cadmium)

#Make a basic boxplot based on ggplot.
boxplot.Cd<-ggplot(data = geochem.df, aes(x=Landcover.Name_2012, y=Cd_ppm, fill=Depth)) +
  geom_boxplot()
boxplot.Cd

#Did the Cd level exceed the govt limit of 5 ppm?
#Now let's polish it a bit

boxplot.Cd.p<- boxplot.Cd +
  labs(title="Distribution of Cd by Land Cover", x="Land cover", y="Concentration (ppm)") + #change the labels for the axes
  theme(axis.text.x=element_text(angle=90, hjust=1))+   #change the angle of the names of the land cover
  coord_cartesian(ylim = c(0,0.35)) +               	#rescale
  geom_hline(yintercept=0.01, linetype="dashed", color = "red")  #add the method detection limit

plot(boxplot.Cd.p)

#Try creating a boxplot for the other elements.

#Part 5 - PERMANOVA

#Load the “vegan” library
library(vegan)

#Isolate the elements that we want to analyze
geochem.metals <- geochem.df%>% select(Ag_ppm:Zr_ppm)
View(geochem.metals)


#Perform PERMANOVA

#PERMANOVA is a non-parametric alternative to MANOVA. 
#It uses a permutation-based technique, where the distribution of the test statistic is assessed by repeatedly relabelling observational units. 
#It doesn’t require normality and homogeneity of variances. This is useful to us since our data is highly skewed (a common occurrence in ecological data), as shown in our histogram.

#For this particular exercise, we are going to use the Bray-Curtis method of calculating the distance matrix. #Other methods also exist such as Euclidean but Bray-Curtis is most robust and is the default in R.

#Sources: https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/
#https://sites.ualberta.ca/~ahamann/teaching/renr690/labs/Lab5.pdf) 

geochem.p <- adonis2(geochem.metals ~ Landcover.Name_2012, data = geochem.df, method="bray", permutations =999) #you may change the number of permutations, but the defaul in R is permutations = 999
geochem.p

#Since p-value < 0.05, it suggests that Landcover.Name_2012 is statistically significant. 
#There is evidence to reject the null hypothesis.
#We will conduct a posthoc test using pairwise comparisons for PERMANOVA because we want to know which groups are different.

#Install "EcolUtils" packages. It's not in CRAN so we need to download it from GitHub.
devtools::install_github("GuillemSalazar/EcolUtils")   
library(EcolUtils)

geochem.df$Landcover.Name_2012=as.factor(geochem.df$Landcover.Name_2012)
geochem.pairwise<-adonis.pair(vegdist(geochem.metals), geochem.df$Landcover.Name_2012)
geochem.pairwise

#Filter values with P.value.corrected < 0.05
geochem.filter <- filter(geochem.pairwise, P.value.corrected < 0.05)
View(geochem.filter)


#Part 6 - Principal Component Analysis (PCA)

#We want to do a PCA of Depth A so we have to remove values in Depth B.
geochem.a <- geochem.df[-c(349:1028), ]
View(geochem.a)

#Select elements of interest. We want to look at metals that are non-essential and toxic even at low levels.
colnames(geochem.a)
geochem.met <- geochem.a[, c("Hg_ppm", "Cd_ppm", "Pb_ppm", "Sb_ppm", "Tl_ppm", "U_ppm", "Landcover.Name_2012")]

#Correlation matrix  
res1 <- cor(geochem.met[1:6]) 
res1

#Prepare a corrplot
library(corrplot)
corrplot(res1) 


#Perform PCA
res.pca <- prcomp(geochem.met[1:6], scale = TRUE) #performing a PCA on the first 6 columns and variables are to be standardized
res.pca
summary(res.pca)


#Plot our principal components
library(ggplot2) 
library(ggfortify)
autoplot(res.pca) #a scatter plot of the PCA results using the default settings.

autoplot(res.pca, data = geochem.met, colour = 'Landcover.Name_2012') #the data points colored according to the Landcover.Name_2012 variable

autoplot(res.pca, data = geochem.met, colour = 'Landcover.Name_2012', loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3) #with the PCA loadings plotted as arrows. The arrows are colored blue, labeled, and with a size of 3 .


autoplot(res.pca, data = geochem.met, colour = 'Landcover.Name_2012',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") # add the x and y axes passing through the origin


# PCA result
# Eigenvalues
library(factoextra)
eig.val <- get_eigenvalue(res.pca) #The eigenvalues represent the amount of variance explained by each PC in the PCA
eig.val

# Results for Variables (Features)
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Results for Individuals (Observations)
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation
