library(readxl)
library(tidyverse)

#transforming data complied by Jennifer and I
classdata <- read_excel("Desktop/2020hsmcResearch/classdata.xlsx")
names(classdata)[5] <- "Rating"
classdata <- classdata %>% filter(Category == 4, Rating != -1, Minute != -1, Second != -1)
classdata <- transform(classdata, Category = as.numeric(classdata$Category))
View(classdata)

#transforming data complied by Zeliha
dataframe <- read_excel("Desktop/2020hsmcResearch/dataframe.xlsx")
dataframeRatings <- read_excel("Desktop/dataframeRatings.xlsx")
names(dataframe)[3] <- "Minute"
names(dataframe)[4] <- "Second"
names(dataframe)[5] <- "Category"
dataframe <- dataframe %>% select(-LessonCode) %>% 
  left_join( dataframeRatings, by = "Lesson") %>% 
  filter(Category == 4, Rating != -1)
View(dataframe)

#creating data frame with the number of teacher moves for each classroom
classdata.full <- full_join(classdata, dataframe) 
X<- data.frame(table(classdata.full$Lesson))
names(X)[1] <- "Lesson"
View(X)

#creating the data frame with the rating of each classroom
classdata.partial <- classdata.full %>% distinct(Lesson, .keep_all = TRUE)
w <- select(classdata.partial, Lesson, Rating)
View(w)

#combining dataframes
master <- inner_join(X, w, by = "Lesson") %>% mutate(Rating.f = factor(Rating))

summary(model1 <- glm(Freq~Rating, family = poisson(), data = master))
summary(model2 <- glm(Freq~Rating.f, family = poisson(), data = master))
betas = coef(model2)
exp(betas[1]+betas[2:4])
round (predict(model2,type = "response"), 2) 
