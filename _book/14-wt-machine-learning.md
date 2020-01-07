# Walkthrough 8: Predicting Students' Final Grades Using Machine Learning Methods {#c14}

## Background

One area of interest is the delivery of online instruction, which is becoming more prevalent: in 2007, over 3.9 million U.S. students were enrolled one or more online courses [@allen2008]. With the dawn of online learning comes an abundance of new educational tools to facilitate that learning. Indeed, online learning interfaces are used to facilitate the submission of assignments and quizzes in courses in which students and instructor meet face-to-face, but these interfaces are also used in fully online courses to deliver all instruction and assessment. 

In a face-to-face classroom, an educator might count on behavioral cues to help them effectively deliver instruction. However, one constraint of online education is that educators do not have access as readily to the behavioral cues that can be essential to effective face-to-face instruction. For example, in a face-to-face classroom, cues such as a student missing class repeatedly or many students seeming distracted during a lecture can trigger a shift in the delivery of instruction. While technology is rapidly developing, many educators find themselves looking for ways to understand and support students online in the same way that face-to-face instructors would. Educational technology affords unique opportunities to support student success online because it provides new methods of collecting and storing data. 

Indeed, online learning management systems often automatically track several types of student interactions with the system and feed that data back to the course instructor. For example, an instructor might be able to quickly see how many students logged into their course on a certain day, or they might see how long students engaged with a posted video before pausing it or logging out. The collection of this data is met with mixed reactions from educators. Some are concerned that data collection in this manner is intrusive, but others see a new opportunity to support students in online contexts in new ways. As long as data are collected and utilized responsibly, data collection can support student success.

One meaningful perspective from which to consider students' engagement with online courses is related to their motivation to achieve. More specifically, it is important to consider how and why students are engaging with the course. Considering the psychological mechanisms behind achievement is valuable because doing so may help to identify meaningful points of intervention for educators and for researchers and administrators in online *and* face-to-face courses interested in the intersection between behavioral trace measures and students' motivational and emotional experiences in such courses.

In this walkthrough, we examine the educational experiences of students in online science courses at a virtual middle school in order to characterize their motivation to achieve and their tangible engagement with the course in terms of behavioral trace measures. To do so, we use a robust data set, which includes self-reported motivation as well as behavioral trace data collected from a learning management system (LMS) to identify predictors of final course grade. Our work examines the idea of educational success in terms of student interactions with an online science course.

We explore the following four questions:

1. Is motivation more predictive of course grades as compared to other online indicators of engagement?
2. Which types of motivation are most predictive of achievement?
3. Which types of trace measures are most predictive of achievement?
4. How does a random forest compare to a simple linear model (regression)?

## Information about the dataset 

This dataset came from 499 students enrolled in online middle school science courses in 2015-2016. The data were originally collected for use as a part of a research study, though the findings have not been published anywhere yet.

The setting of this study was a public provider of individual online courses in a Midwestern state. In particular, the context was two semesters (Fall and Spring) of offerings of five online science courses (Anatomy & Physiology, Forensic Science, Oceanography, Physics, and Biology), with a total of 36 classes. 

Specific information in the dataset included:

- a pre-course survey students completed about their self-reported motivation in science — in particular, their perceived competence, utility value, and interest
- the time students spent on the course (obtained from the learning management system (LMS), Blackboard
- students' final course grades 
- students' involvement in discussion forums

For discussion board responses, we were interested in calculating the number of posts per student and understanding the emotional tone of the discussion board posts. We used the Linguistic Inquiry and Word Count (LIWC; Pennebaker, Boyd, Jordan, & Blackburn, 2015) tool to calculate the number of posts per student and to categorize the emotional tone (positive or negative) and topics of those posts. Those linguistic categorization was conducted after the data was gathered from the discussion posts, but is not replicated here to protect the privacy of the students' posts. Instead, we present the already-categorized discussion board data, in its ready-to-use format. Thus, in the dataset used in this walkthrough, we will see pre-created variables for the mean levels of students' cognitive processing, positive emotions, negative emotions, and social-related discourse.

At the beginning of the semester, students were asked to complete the pre-course survey about their perceived competence, utility value, and interest. At the end of the semester, the time students spent on the course, their final course grades, and the contents of the discussion forums were collected.

In this walkthrough, we used the R package **caret** to carry out the analyses.

## Selecting an analysis

### Defining our research question
When you begin a new project, there are often many approaches to analyzing data and answering questions you might have about it. Some projects have a clearly defined scope and question to answer. This type of project is characterized by 1) a defined number of variables (data inputs) and 2) specific directional hypotheses. For example, if we are studying the effect of drinking coffee after dinner on ability to quickly fall asleep, we might have a very specific directional hypothesis: we expect that drinking coffee after dinner would decrease the ability to fall asleep quickly. In this case, we might collect data by having some people drink coffee and having other people drink nothing or an herbal tea before bed. We could monitor how quickly people from each group fall asleep. Since we collected data from two clearly defined groups, we can then do a statistical analysis that compares the amount of time it takes to fall asleep for each group. One option would be a test called a t-test, which we could use to see if there is a significant difference in the average amount of minutes to fall asleep for the group. This approach works very well in controlled experimental situations, especially when we can change only one thing at a time (in our coffee example, the only thing we changed was the coffee-drinking behavior of our participants - all other life conditions were held equal for both groups). Rarely are educational data projects as clear-cut and simple.

For this walkthrough, we have many sources of data - survey data, learning management system data, discussion forum data, and academic achievement data as measured by final course grades. Luckily, having too much data is what we call a "good problem." In our coffee example above, we had one really specific idea that we wanted to investigate - does coffee affect time taken to fall asleep? In this walkthrough we have many ideas we are curious to explore: the relationships among motivation, engagement in the course (discussion boards, time spent online in the course site), and academic achievement. If we wanted to tackle a simpler problem, we could choose just one of these relationships. For example, we could measure whether students with high motivation earn higher grades than students with low motivation. However, we are being a bit more ambitious than that here - we are interested in understanding the complex relationships among the different types of motivation. Rather than simply exploring whether A affects B, we are interested in the nuances: we suspect that *many* factors affect B, and we would like to see which of those factors has most relative importance. To explore this idea, we will use a machine learning approach.

### Predictive analytics and machine learning
A buzzword in education software spheres these days is "predictive analytics." Administrators and educators alike are interested in applying the methods long utilized by marketers and other business professionals to try to determine what a person will want, need, or do next. "Predictive analytics" is a blanket term that can be used to describe any statistical approach that yields a prediction. We could ask a predictive model: "What is the likelihood that my cat will sit on my keyboard today?" and, given enough past information about your cat's computer-sitting behavior, the model could give you a probability of that computer-sitting happening today. Under the hood, some predictive models are not very complex. If we have an outcome with two possibilities, a logistic regression model could be fit to the data in order to help us answer the cat-keyboard question. In this chapter, we'll compare a machine learning model to another type of regression: multiple regression. We want to make sure to fit the simplest model as possible to our data. After all, the effectiveness in predicting the outcome is really the most important thing: not the fanciness of the model.
    
Data collection is an essential first step in any type of machine learning or predictive analytics. It is important to note here that machine learning only works effectively when (1) a person selects variables to include in the model that are anticipated to be related to the outcome and (2) a person correctly interprets the model's findings. There is an adage that goes, "garbage in, garbage out." This holds true here: if we do not feel confident that the data we collected are accurate, no matter what model we build, we will not be able to be confident in our conclusions. To collect good data, we must first clarify what it is that we want to know (i.e., what question are we really asking?) and what information we would need in order to effectively answer that question. Sometimes, people approach analysis from the opposite direction - they might look at the data they have and ask what questions could be answered based on that data. That approach is okay - as long as you are willing to acknowledge that sometimes the pre-existing dataset may *not* contain all the information you need, and you might need to go out and find additional information to add to your dataset to truly answer your question.
    
When people talk about "machine learning," you might get the image in your head of a desktop computer learning how to spell. You might picture your favorite social media site showing you advertisements that are just a little too accurate. At its core, what machine learning really is is the process of "showing" your statistical model only some of the data at once, and training the model to predict accurately on that training dataset (this is the "learning" part of machine learning). Then, the model as developed on the training data is shown new data - data you had all along, but hid from your computer initially - and you see how well the model that you developed on the training data performs on this new testing data. Eventually, you might use the model on entirely new data.  

## Information on random forests
For our analyses, we used Random Forest modeling [@breiman2001]. Random forest is an extension of decision tree modeling, whereby a collection of decision trees are simultaneously "grown" and are evaluated based on out-of-sample predictive accuracy [@breiman2001].  Random forest is random in two main ways: first, each tree is only allowed to "see" and split on a limited number of predictors instead of all the predictors in the parameter space; second, a random subsample of the data is used to grow each individual tree, such that no individual case is weighted too heavily in the final prediction. 

One thing about random forest that makes it quite different from other types of analysis we might do is that here, we are giving the computer a large amount of information and asking it to find connections that might not be immediately visible to the naked human eye. This is great for a couple of reasons. First, while humans are immensely creative and clever, we are not immune to biases. If we are exploring a dataset, we usually come in with some predetermined notions about what we think is true, and we might (consciously or unconsciously) seek evidence that supports the hypothesis we privately hold. By setting the computer loose on some data, we can learn that there are connections between areas that we did not expect. We must also be ready for our hypotheses to not be supported! Random forest is particularly well-suited to the research questions explored here because we do not have specific directional hypotheses. Machine learning researchers talk about this as "exploring the parameter space" - we want to see what connections exist, and we acknowledge that we might not be able to accurately predict all the possible connections. Indeed, we expect - and hope - that we will find surprising connections. 

Whereas some machine learning approaches (e.g., boosted trees) would utilize an iterative model-building approach, random forest estimates all the decision trees at once. In this way, each tree is independent of every other tree. Thus, the random forest algorithm provides a robust regression approach that is distinct from other modeling approaches. The final random forest model aggregates the findings across all the separate trees in the forest in order to offer a collection of "most important" variables as well as a percent variance explained for the final model.

500 trees were grown as part of our random forest. We partitioned the data before conducting the main analysis so that neither the training nor the testing data set would be disproportionately representative of high-achieving or low-achieving students. The training data set consisted of 80% of the original data (n = 400 cases), whereas the testing data set consisted of 20% of the original data (n = 99 cases). We built our random forest model on the training data set, and then evaluated the model on the testing data set. Three variables were tried at each node.

Note that the random forest algorithm does not accept cases with missing data, and so we deleted cases listwise if data were missing. This decision eliminated 51 cases from our original data set, to bring us to our final sample size of 499 unique students. If you have a very small dataset with a lot of missing data, the random forest approach may not be well suited for your goals – you might consider a linear regression instead. 

A random forest is well suited to the research questions that we had here because it allows for nonlinear modeling. We hypothesized complex relationships between students' motivation, their engagement with the online courses, and their achievement. For this reason, a traditional regressive or structural equation model would have been insufficient to model the parameter space we were interesting in modeling. Our random forest model had one outcome and eleven predictors. 

One term you will hear used in machine learning is "tuning parameter." People often think of tuning parameters as knobs or dials on a radio: they are features of the model that can be adjusted to get the clearest signal. A common tuning parameter for machine learning models is the number of variables considered at each split [@kuhn2008]; we considered three variables at each split for this analysis.  

The outcome was the final course grade that the student earned. The predictor variables included motivation variables (interest value, utility value, and science perceived competence) and trace variables (the amount of time spent in the course, the course name, the number of discussion board posts over the course of the semester, the mean level of cognitive processing evident in discussion board posts, the positive emotions evident in discussion board posts, the negative emotions evident in discussion board posts, and the social-related discourse evident in their discussion board posts). We used this random forest model to address all three of our research questions.

To interpret our findings, we will  three main factors: (1) predictive accuracy of the random forest model, (2) variable importance, and (3) variance explained by the final random forest model. In this walkthrough, we will the R package **caret** to carry out the analyses.

## Analysis



First, we will load the data. Our data is stored in the `dataedu` package that is part of this book. Within that package, the data is stored as an .rda file. 


```r
#loading the data from the .rda file and storing it as an object named 'data'

data <- dataedu::sci_mo_data
```


It's a good practice to take a look at the data and make sure it looks the way you expect it to look. R is pretty smart, but sometimes we run into issues like column headers being read as datapoints. By using the `glimpse` function from the `dplyr` package, we can quickly skim our data and see whether we have all the right variables and datapoints. Remember that the `dplyr` package loads automatically when we load the `tidyverse` library, so there is no need to call the `dplyr` package separately. Now, we'll glimpse the data.


```r
glimpse(sci_mo_data)
```

```
## Observations: 550
## Variables: 20
## $ pre_int           <dbl> 5.0, 4.2, 4.6, 5.0, 4.0, 4.0, 2.4, 4.2, 4.0, 4.2, 4…
## $ pre_uv            <dbl> 5.000000, 5.000000, 4.333333, 4.333333, 3.666667, 4…
## $ pre_percomp       <dbl> 4.0, 4.5, 4.5, 4.0, 3.5, 4.0, NA, 2.5, 4.5, 3.0, 3.…
## $ time_spent        <dbl> 2167.9669, 1550.4668, 850.7329, 3067.4676, 1800.217…
## $ course_ID         <chr> "AnPhA-S116-01", "FrScA-S116-01", "FrScA-S216-01", …
## $ final_grade       <dbl> 93.753600, 95.890411, 88.356164, 91.337600, 88.7202…
## $ subject           <chr> "AnPhA", "FrScA", "FrScA", "AnPhA", "OcnA", "BioA",…
## $ enrollment_reason <chr> "Learning Preference of the Student", "Course Unava…
## $ semester          <chr> "S116", "S116", "S216", "S216", "S116", "S116", "S2…
## $ enrollment_status <chr> "Approved/Enrolled", "Approved/Enrolled", "Approved…
## $ cogproc           <dbl> 10.127500, 14.702162, 15.465128, 11.697391, 13.9536…
## $ social            <dbl> 7.550435, 5.903571, 5.729250, 6.001304, 7.981111, 1…
## $ posemo            <dbl> 3.905882, 3.278378, 3.544242, 1.150000, 5.439643, 7…
## $ negemo            <dbl> 1.5700000, 0.8097500, 1.0472973, 0.7665714, 0.60095…
## $ n                 <dbl> 25, 40, 42, 26, 14, 6, 19, 39, 21, 15, 30, 46, 17, …
## $ section           <chr> "01", "01", "01", "01", "01", "01", "01", "01", "01…
## $ post_int          <dbl> NA, 2.25, NA, NA, 5.00, NA, NA, NA, NA, NA, NA, 5.0…
## $ post_uv           <dbl> NA, 4.333333, NA, NA, 3.666667, NA, NA, NA, NA, NA,…
## $ post_percomp      <dbl> NA, 2.0, NA, NA, 4.5, NA, NA, NA, NA, NA, NA, 4.0, …
## $ WC                <dbl> 74.24000, 64.30769, 96.28571, 87.77551, 67.53846, 6…
```
Scanning the data we glimpsed, we see that we have 662 observations and 111 variables. Many of these variables - everything below "WC" except the variable "n" -  are related to the text content of the discussion board posts. Our analysis here is not focused on the specifics of the discusion board posts, so we will select just a few variables from the LIWC analysis. If you're interested in learning more about analyzing text, the text analysis walkthrough in this volume will be a good place to start. 

As is the case with many datasets you'll work with in education contexts, there is lots of great information in this dataset - but we won't need all of it. Even if your dataset has many variables, for most analyses you will find that you are only interested in some of them. There are statistical reasons not to include twenty or more variables in a data analysis, and the quick explanation of the reason why is that at a certain point, adding more variables will *appear* to make your analysis more accurate, but will in fact obscure the truth from you. It's generally a good practice to select a few variables you are interested in and go from there. As we discussed above, the way to do this is to start with the research questions you are trying to answer. Since we are interested in data from one specific semester, we'll need to narrow down the data to make sure that we only include datapoints relevant to that semester.

Thus, we will *filter* the data to include only the data from one that semester, and then *select* variables of interest. For each step, we will save over the previous version of the "data" object so that our working environment doesn't get cluttered with each new version of the dataset. Keep in mind that the original data will stay intact, and that any changes we make to it within R will not overwrite that original data (unless we tell R to specifically save out a new file with exactly the same name as the original file). Changes we make within our working environment are all totally reversible. 

Below, we will *filter* to remove all the datapoints from the spring 2017 semester (indicated with a value of "S217" for the "semester" variable). We use the "!" to indicate that we want to keep all datapoints EXCEPT the datapoints that have a value of "S217" for the semester variable. Then, we will *select* only the variables we are interested in: motivation, time spent in the course, grade in the course, subject, enrollment information, positive and negative emotions, cognitive processing, and the number of discussion board posts.


```r
#filtering the data to only include 2016 data. 

data <- 
    data %>% 
    filter(semester != "S217")
```

```
## filter: no rows removed
```

```r
#selecting only the variables we are interested in: 
data <- 
    data %>%
    select(
        pre_int,
        pre_uv,
        pre_percomp,
        time_spent,
        final_grade,
        subject,
        enrollment_reason,
        semester,
        enrollment_status,
        cogproc,
        social,
        posemo,
        negemo,
        n
    )
```

```
## select: dropped 6 variables (course_ID, section, post_int, post_uv, post_percomp, …)
```

### Use of caret

Here, we remove observations with missing data (per our note above about random forests requiring complete cases).


```r
#checking how many rows are in our dataset
nrow(data)
```

```
## [1] 550
```

```r
    #we see that we have 550 rows from spring 2017

#calling the na.omit function to eliminate ANY rows that have ANY missing data
data <- na.omit(data)

#checking whether our na.omit call worked as expected
nrow(data)
```

```
## [1] 494
```

```r
    #after running the code above, we see that we now have 499 rows - this is as we expected
```

First, machine learning methods often involve using a large number of variables. Oftentimes, some of these variables will not be suitable to use: they may be highly correlated with other variables, for instance, or may have very little - or no - variability. Indeed, for the data set used in this study, one variable has the same (character string) value for all of the observations. We can detect this variable and any others using the following function:


```r
#we run the nearZeroVar function to determine if there are variables with NO variability
nearZeroVar(data, saveMetrics = TRUE)
```

```
##                   freqRatio percentUnique zeroVar   nzv
## pre_int            1.500000     2.8340081   FALSE FALSE
## pre_uv             1.474359     2.4291498   FALSE FALSE
## pre_percomp        1.225225     1.4170040   FALSE FALSE
## time_spent         1.250000    65.5870445   FALSE FALSE
## final_grade        1.250000    63.1578947   FALSE FALSE
## subject            1.360544     1.0121457   FALSE FALSE
## enrollment_reason  3.229885     1.0121457   FALSE FALSE
## semester           1.344660     0.6072874   FALSE FALSE
## enrollment_status  0.000000     0.2024291    TRUE  TRUE
## cogproc            1.000000    61.1336032   FALSE FALSE
## social             1.000000    62.1457490   FALSE FALSE
## posemo             2.200000    63.7651822   FALSE FALSE
## negemo             4.200000    61.3360324   FALSE FALSE
## n                  1.240000     9.7165992   FALSE FALSE
```
After conducting our zero variance check, we want to scan the "zeroVar" column to see if any of our variables failed this check. If we see any "TRUE" values for "zeroVar," that means we should look more closely at that variable.

In the nearZeroVar function we just ran, we see a result in the ZeroVar column of "TRUE" for the `enrollment_status` variable.  If we look at `enrollment_status`, we will see that it is "Approved/Enrolled" for *all* of the students.  When we use variables with no variability in certian models, it may cause some problems, and so we remove it first.


```r
#Taking the dataset and re-saving it as the same dataset, but without the enrollment status variable
data <- 
    data %>% 
    select(-enrollment_status)
```

```
## select: dropped one variable (enrollment_status)
```

Note that many times you may wish to pre-process the variables, such as by centering or scaling them. Often the data will come to you in a format that is not ready for immediate analysis, as we have discussed elsewhere in the book. For our current dataset, we could work on pre-processing with code like you will see below. We set this next code chunk up to not run here (if you are viewing the book online), as we will do this analysis with the variables' original values.


```r
#example pre-processing step: manipulating the dataset 'data' so that if a variable is numeric, its format will now be scale
data <- 
    data %>% 
    mutate_if(is.numeric, scale)
```

As another pre-processing step, we want to make sure our text data is in a format that we can easily evaluate. To facilitate that, we want to make character string variables into factors.


```r
#converting the text (character) variables in our dataset into factors
data <- 
    data %>% 
    mutate_if(is.character, as.factor)
```

```
## mutate_if: converted 'subject' from character to factor (0 new NA)
```

```
##            converted 'enrollment_reason' from character to factor (0 new NA)
```

```
##            converted 'semester' from character to factor (0 new NA)
```

Now, we will prepare the **train** and **test** datasets, using the caret function for creating data partitions. Here, the **p** argument specifies what proportion of the data we want to be in the **training** partition. Note that this function splits the data based upon the outcome, so that the training and test data sets will both have comparable values for the outcome. This means that since our outcome is final grade, we are making sure that we don't have either a training or testing dataset that has too many good grades - or too many bad grades. Note the `times = 1` argument; this function can be used to create *multiple* train and test sets, something we will describe in more detail later. Before we create our training and testing datasets, we want to initiate a process called "setting the seed." This means that we are ensuring that if we run this same code again, we will get the same results in terms of the data partition. The seed can be any number that you like - some people choose their birthday or another meaningful number. The only constraint is that wehn you open the same code file again to run in the future, you do not change the number you selected for your seed. This ensures your code is reproducible. In fact, it ensures that anyone who runs the same code file on any computer, anywhere, will get the same result. With that background information, try running the code chunk below.


```r
#First, we set a seed to ensure the reproducibility of our data partition.
set.seed(62020)

#we create a new object called trainIndex that will take 80 percent of the data
trainIndex <- createDataPartition(data$final_grade,
                                  p = .8, 
                                  list = FALSE,
                                  times = 1)

#We add a new variable to our dataset, temporarily:
    #this will let us select our rows according to their row number
    #we populate the rows with the numbers 1:499, in order

data <- 
    data %>% 
    mutate(temp_id = 1:494)
```

```
## mutate: new variable 'temp_id' with 494 unique values and 0% NA
```

```r
#we filter our dataset so that we get only the 
    #rows indicated by our "trainIndex" vector
data_train <- 
    data %>% 
    filter(temp_id %in% trainIndex)
```

```
## filter: removed 96 rows (19%), 398 rows remaining
```

```r
#we filter our dataset in a different way so that we get only  the rows NOT in our "trainIndex" vector
    #adding the ! before the temp_id variable achieves the opposite of what we did in the line of code above

data_test <- 
    data %>% 
    filter(!temp_id %in% trainIndex)
```

```
## filter: removed 398 rows (81%), 96 rows remaining
```

```r
#We delete the temp_id variable from (1) the original data, (2) the portion of the original data we marked as training, and (3) the portion of the original data we marked as testing, as we no longer need that variable

data <- 
    data %>% 
    select(-temp_id)
```

```
## select: dropped one variable (temp_id)
```

```r
data_train <- 
    data_train %>% 
    select(-temp_id)
```

```
## select: dropped one variable (temp_id)
```

```r
data_test <- 
    data_test %>% 
    select(-temp_id)
```

```
## select: dropped one variable (temp_id)
```

Finally, we will estimate the models.

Here, we will use the train function, passing *all* of the variables in the data frame (except for the outcome, or dependent variable, `final_grade`) as predictors.

The predictor variables include three indicators of motivation: interest in the course (pre_int), perceived utility value of the course (pre_uv), and perceived competence for the subject matter (pre_percomp). There are a few predictor variables that help differentiate between the different courses in the dataset: subject matter of the course (subject), reason the student enrolled in the course (enrollment_reason), and semester in which the course took place (semester). We have a predictor variable that indicates the amount of time each student spent engaging with the online learning platform of the course (time_spent).

We also have a number of variables associated with the discussion board posts from the course. Specifically, the variables include the average level of cognitive processing in the discussion board posts (cogproc), the average level of social (rather than academic) content in the discussion board posts (social), the positive and negative emotions evident in the discussion board posts (posemo and negemo), and finally, the number of discussion board posts in total (n). We are using all those variables discussed in this paragraph to predict the outcome of the final grade in the course (final_grade).

Note that you can read more about the specific random forest implementation chosen [here](http://topepo.github.io/caret/train-models-by-tag.html#random-forest). To specify that we want to predict the outcome using every variable except the outcome itself, we use the formulation (outcome ~ .,). R interprets this code as: predict the outcome using all the variables except outcome itself. The outcome always comes before the `~`, and the `.` that we see after the `~` means that we want to use all the rest of the variables. An alternative specification of this model would be to write (outcome ~ predictor1, predictor2). Anything that follows the `~` and precedes the comma is treated as predictors of the outcome.

Here, we set the seed again, to ensure that our analysis is reproducible. This step of setting the seed is especially important due to the "random" elements of random forest, because it's likely that the findings would change (just slightly) if the seed were not set. As we get into random forest modeling, you might notice that the code takes a bit longer to run. This is normal - just think of the number of decision trees that are "growing!"


```r
#setting a seed for reproducibility
set.seed(62020)

#we run the model here
rf_fit <- train(final_grade ~ .,
                data = data_train,
                method = "ranger")

#here, we get a summary of the model we just built
rf_fit
```

```
## Random Forest 
## 
## 398 samples
##  12 predictor
## 
## No pre-processing
## Resampling: Bootstrapped (25 reps) 
## Summary of sample sizes: 398, 398, 398, 398, 398, 398, ... 
## Resampling results across tuning parameters:
## 
##   mtry  splitrule   RMSE      Rsquared   MAE      
##    2    variance    15.40824  0.6062405  11.338374
##    2    extratrees  17.51403  0.5443134  12.555711
##   10    variance    13.46816  0.6535352   9.806955
##   10    extratrees  13.58057  0.6767111   9.952872
##   19    variance    13.55542  0.6443721   9.615483
##   19    extratrees  13.07647  0.6814068   9.541441
## 
## Tuning parameter 'min.node.size' was held constant at a value of 5
## RMSE was used to select the optimal model using the smallest value.
## The final values used for the model were mtry = 19, splitrule = extratrees
##  and min.node.size = 5.
```

We have some results! First, we see that we have 400 samples, or 400 observations, the number in the train data set. No pre-processing steps were specified in the model fitting, but note that the output of `preProcess` can be passed to `train()` to center, scale, and transform the data in many other ways. Next, in our example, a resampling technique has been used. This resampling is not for validating the model (per se), but is rather for selecting the tuning parameters - the options that need to be specified as a part of the modeling. These parameters can be manually provided, or can be estimated via strategies such as the bootstrap resample (or *k*-folds cross validation).

As we interpret these findings, we are looking to minimize the error (RMSE) and maximize the variance explained (rsquared).

It appears that the model with the value of the **mtry** tuning parameter equal to 19 seemed to explain the data best, the **splitrule** being "extratrees", and **min.node.size** held constant at a value of 5. We know this model fits best because the RMSE is the lowest of the options (13.13) and the Rsquared is the highest of the options (.582).

The value of resampling here is that it allows for higher accuracy of the model [@james2013]. Without resampling (bootstrapping or cross-validation), the variance would be higher and the predictive accuracy of the model would be lower.

Let's see if we end up with slightly different values if we change the resampling technique to cross-validation, instead of bootstrap resampling. We set a seed again here, for reproducibility.


```r
set.seed(62020)

train_control <-
    trainControl(method = "repeatedcv",
                 number = 10,
                 repeats = 10)

rf_fit1 <-
    train(final_grade ~ .,
          data = data_train,
          method = "ranger",
          trControl = train_control)

rf_fit1
```

```
## Random Forest 
## 
## 398 samples
##  12 predictor
## 
## No pre-processing
## Resampling: Cross-Validated (10 fold, repeated 10 times) 
## Summary of sample sizes: 358, 359, 358, 358, 358, 358, ... 
## Resampling results across tuning parameters:
## 
##   mtry  splitrule   RMSE      Rsquared   MAE      
##    2    variance    14.84149  0.6245157  10.914358
##    2    extratrees  16.87348  0.5851049  12.238703
##   10    variance    12.59610  0.6771022   9.180824
##   10    extratrees  12.84618  0.6956843   9.452986
##   19    variance    12.52421  0.6749867   8.947893
##   19    extratrees  12.41260  0.6945521   9.067094
## 
## Tuning parameter 'min.node.size' was held constant at a value of 5
## RMSE was used to select the optimal model using the smallest value.
## The final values used for the model were mtry = 19, splitrule = extratrees
##  and min.node.size = 5.
```

When we look at this output, we are looking to see which values of the various tuning parameters were selected. We see at the bottom of the output above that the value of **mtry** was 19, the split rule was "extratrees," and the minimum node size is 5. We let this model explore which value of **mtry** was best and to explore whether extra trees or variance was a better split rule, but we forced all iterations of the model to a minimum node size of five (so that minimum node size value in the output shouldn't be a surprise to us). When we look at the bottom row of the output, it shows the final values selected for the model. We see also that this row has the lowest RMSE and highest Rsquared value, which means it has the lowest error and highest predictive power. 

We won't dive into the specifics of the statistics behind these decisions right now, but next we will try adjusting a few different parts of the model to see whether our performance improves. For a detailed statistical explanation of random forest modeling, including more about **mtry** and tuning a model, please see Chapter 8 in the book "An Introduction to Statistical Learning with Applications in R" [@james2013]. 

What would happen if we do not fix **min.node.size** to five? We're going to let **min.node.size** change and let **mtry** change as well.

Let's create our own grid of values to test for **mtry** and **min.node.size**. We'll stick with the default bootstrap resampling method to choose the best model. We will randomly choose some values to use for **mtry**, including the three that were used previously (2, 10, and 19). Let's try 2, 3, 7, 10, and 19.


```r
#create a grid of different values of mtry, different splitrules, and different minimum node sizes to test
tune_grid <-
    expand.grid(
        mtry = c(2, 3, 7, 10, 19),
        splitrule = c("variance", "extratrees"),
        min.node.size = c(1, 5, 10, 15, 20)
    )

#set a seed
set.seed(62020)

#fit a new model, using the tuning grid we created above
rf_fit2 <-
    train(final_grade ~ .,
          data = data_train,
          method = "ranger",
          tuneGrid = tune_grid)

rf_fit2
```

```
## Random Forest 
## 
## 398 samples
##  12 predictor
## 
## No pre-processing
## Resampling: Bootstrapped (25 reps) 
## Summary of sample sizes: 398, 398, 398, 398, 398, 398, ... 
## Resampling results across tuning parameters:
## 
##   mtry  splitrule   min.node.size  RMSE      Rsquared   MAE      
##    2    variance     1             15.32302  0.6110390  11.256089
##    2    variance     5             15.41644  0.6033990  11.358638
##    2    variance    10             15.58342  0.5963452  11.524727
##    2    variance    15             15.76218  0.5868123  11.680710
##    2    variance    20             15.89639  0.5806568  11.794798
##    2    extratrees   1             17.35985  0.5535581  12.406218
##    2    extratrees   5             17.50843  0.5441410  12.558312
##    2    extratrees  10             17.68886  0.5307434  12.730148
##    2    extratrees  15             17.94445  0.5172112  12.937712
##    2    extratrees  20             18.09510  0.5111653  13.051717
##    3    variance     1             14.52728  0.6280811  10.653547
##    3    variance     5             14.57908  0.6260363  10.738165
##    3    variance    10             14.71290  0.6190406  10.886307
##    3    variance    15             14.87838  0.6115016  11.045861
##    3    variance    20             15.10341  0.5993341  11.247855
##    3    extratrees   1             15.93170  0.6051503  11.382563
##    3    extratrees   5             16.11388  0.5974488  11.558824
##    3    extratrees  10             16.40912  0.5831345  11.826247
##    3    extratrees  15             16.68300  0.5716373  12.062076
##    3    extratrees  20             16.85874  0.5650996  12.213533
##    7    variance     1             13.61613  0.6524226   9.985055
##    7    variance     5             13.70025  0.6472675  10.044141
##    7    variance    10             13.78366  0.6435253  10.142727
##    7    variance    15             13.92684  0.6361402  10.290956
##    7    variance    20             14.05235  0.6292943  10.401780
##    7    extratrees   1             13.98467  0.6691778  10.169688
##    7    extratrees   5             14.14466  0.6619404  10.311411
##    7    extratrees  10             14.35877  0.6554805  10.526819
##    7    extratrees  15             14.58194  0.6469771  10.718176
##    7    extratrees  20             14.78158  0.6406403  10.884044
##   10    variance     1             13.44769  0.6554058   9.767580
##   10    variance     5             13.46283  0.6546588   9.811348
##   10    variance    10             13.55409  0.6498397   9.900068
##   10    variance    15             13.69744  0.6428284  10.012504
##   10    variance    20             13.80312  0.6368826  10.132796
##   10    extratrees   1             13.51164  0.6793488   9.891000
##   10    extratrees   5             13.58482  0.6759328   9.948678
##   10    extratrees  10             13.72827  0.6714187  10.097534
##   10    extratrees  15             13.93028  0.6641939  10.283625
##   10    extratrees  20             14.10505  0.6573098  10.438428
##   19    variance     1             13.49065  0.6473468   9.563349
##   19    variance     5             13.54851  0.6444930   9.608090
##   19    variance    10             13.62576  0.6401074   9.679835
##   19    variance    15             13.73846  0.6343041   9.773467
##   19    variance    20             13.81754  0.6297773   9.847623
##   19    extratrees   1             13.02666  0.6843333   9.489768
##   19    extratrees   5             13.07677  0.6819550   9.531943
##   19    extratrees  10             13.20664  0.6756461   9.658847
##   19    extratrees  15             13.28154  0.6734305   9.751242
##   19    extratrees  20             13.37376  0.6700883   9.857556
## 
## RMSE was used to select the optimal model using the smallest value.
## The final values used for the model were mtry = 19, splitrule = extratrees
##  and min.node.size = 1.
```

The model with the same values as identified before for **mtry** (19) and **splitrule** (extratrees), but with **min.node.size** equal to 1 (not 5, as before) seems to fit best. We know this model fits best because the RMSE is lowest (13.08) and the variance explained is highest (.58) for this model, though the improvement seems to be fairly small relative to the difference the other tuning parameters seem to make. 

While the output above gives us a good summary of the model, we might want to look more closely at what we found with our rf_fit2 model. The code below is a way for us to zoom in and look specifically at the *final* random forest model generated by our rf_fit2.

In the code chunk below, you'll notice we are selecting the "finalModel" output using a `$` operator rather than the familiar `select`. We cannot use dplyr and the tidyverse here because of the structure of the rf_fit2 object - we have stored a random forest model as a model, so it's not a normal dataframe. Thus, we extract with a `$`. We want to select only the final model used, and not worry about the prior iterations of the model.


```r
#Here, we select the "finalModel" output from the rf_fit2 model
rf_fit2$finalModel
```

```
## Ranger result
## 
## Call:
##  ranger::ranger(dependent.variable.name = ".outcome", data = x,      mtry = min(param$mtry, ncol(x)), min.node.size = param$min.node.size,      splitrule = as.character(param$splitrule), write.forest = TRUE,      probability = classProbs, ...) 
## 
## Type:                             Regression 
## Number of trees:                  500 
## Sample size:                      398 
## Number of independent variables:  19 
## Mtry:                             19 
## Target node size:                 1 
## Variable importance mode:         none 
## Splitrule:                        extratrees 
## OOB prediction error (MSE):       154.3324 
## R squared (OOB):                  0.6921864
```
In looking at this output, we see the same parameters we noted above: **mtry** is 19, the node size is 1, and the split rule is extra trees. We can also note the *OOB prediction error (MSE)*, of 168.92, and the proportion of the variance explained, or R squared, of 0.59. As before, we want the error to be low and the variance explained to be high.

Now that we understand how to develop a basic machine learning model, and how to use different tuning parameters (such as node size and the splitting rule), we can explore some other related themes. We might wonder about how we could examine the predictive accuracy of the random forest model we just developed.

## Examining predictive accuracy on the test data set

What if we use the test data set - data not used to train the model? Below, we'll create a new object that uses the rf_fit2 model we developed above. We will put our testing data through the model, and assign the predicted values to a row called "pred." At the same, time, we'll make a row called "obs" that includes the real final grades that students earned. Later, we'll compare these predicted and observed values to see how well our model did.


```r
set.seed(62020)

##Create a new object for the testing data including predicted values 
data_test_augmented <-
    data_test %>%
    mutate(pred = predict(rf_fit2, data_test),
           obs = final_grade)
```

```
## mutate: new variable 'pred' with 96 unique values and 0% NA
```

```
##         new variable 'obs' with 86 unique values and 0% NA
```

```r
###Transform this new object into a dataframe
defaultSummary(as.data.frame(data_test_augmented))
```

```
##       RMSE   Rsquared        MAE 
## 11.2145423  0.7433124  8.1097816
```

We can compare this to the values above to see how our model performs when given data that was not used to train the model. Comparing the RMSE values, we see that the RMSE is about the same when we use the model on the test data as it was on the training data. We get a value of 12.36 on the test data here, and it was 13.08 on the training data. The Rsquared value is 0.71 here, as compared to the 0.58 we got when we passed the training data through rf_fit2 earlier. 

While we might have expected that the model performance would be worse for the testing data as compared to the training data, we actually are seeing marginal improvements here: the model does better with the test data than with the training data. These results suggest to us that the model is fairly robust, as we get comparable - in fact, improved - results when running the model on data it has never "seen" before (the testing data). This is good news!

### Variable importance measures

One helpful characteristic of random forest models is that we can learn about which variables contributed most strongly to the predictions in our model, across all the trees in our forest.

We can examine two different variable importance measures using the **ranger** method in **caret**.

Note that importance values are not calculated automatically, but that "impurity" or "permutation" can be passed to the `importance` argument in `train()`. See more [here](https://alexisperrier.com/datascience/2015/08/27/feature-importance-random-forests-gini-accuracy.html).

We'll re-run the rf_fit2 model with the same specifications as before, but this time we will add an argument to call the variable importance metric.


```r
#set a seed
set.seed(62020)

#specify the same model as earlier in the chapter (rf_fit2) with the addition of the variable importance metric
rf_fit2_imp <-
    train(
        final_grade ~ .,
        data = data_train,
        method = "ranger",
        tuneGrid = tune_grid,
        importance = "permutation"
    )

#extract the variable importance from this new model
varImp(rf_fit2_imp)
```

```
## ranger variable importance
## 
##                                                      Overall
## n                                                   100.0000
## subjectFrScA                                         18.5374
## time_spent                                           13.5623
## cogproc                                               5.3760
## social                                                2.4427
## negemo                                                2.3617
## pre_uv                                                2.3189
## semesterS216                                          2.1824
## pre_int                                               2.1678
## posemo                                                2.0695
## subjectPhysA                                          1.6834
## subjectOcnA                                           1.1885
## pre_percomp                                           0.7249
## enrollment_reasonOther                                0.6607
## enrollment_reasonScheduling Conflict                  0.4559
## enrollment_reasonCredit Recovery                      0.4224
## semesterT116                                          0.4199
## subjectBioA                                           0.3773
## enrollment_reasonLearning Preference of the Student   0.0000
```

Our results here give us a ranked order list of the variables in the order of their importance. Variables that appear at the top of the list are more important, and variables that appear at the bottom of the list are less important in the specification of our final random forest model. Remember that we are predicting final grade in the course, so this list will tell us which factors were most important in predicting final grade in online science courses. It can be a bit hard to visually scan a variable importance list, so we might be interested in doing a data visualization.

We can visualize this variable importance list with `ggplot`. 


```r
varImp(rf_fit2_imp) %>%
    pluck(1) %>%
    rownames_to_column("var") %>%
    ggplot(aes(x = reorder(var, Overall), y = Overall)) +
    geom_col() +
    coord_flip() +
    theme_dataedu()
```

<img src="14-wt-machine-learning_files/figure-html/unnamed-chunk-17-1.png" width="672" />

Cool! We can now visualize which variables are most important in predicting final grade. 

The first thing we notice is that the variable "n" is the most important. This variable indicates how much students write in their discussion posts. The second most important variable is the amount of time students spend in their course. The third most important variable is `subjectFrScA.` This is one of the course subjects: forensic science. Being enrolled in the forensic science course has a large impact on final grade. That would indicate to us that the forensic science course - more than the other science subjects in this dataset - is strongly correlated with students' final course grades. We can keep scanning down the list to see the other variables that were indicated as less and less important for the model's predictions. Variable importance can thus help us to better understand the inner workings of a random forest model.

Overall, there are some subject level differences in terms of how predictive subject is. Biology (`subjectBioA`) shows up pretty far down the list, whereas Physiology is in the middle (`subjPhysA`) and forensic science is towards the top (`subjectFrScA`). What this tells us is that the course students are in seems to have a different effect on final grade, depending on the course. Perhaps grades should be normalized within subject: would this still be an important predictor if we did that? We won't dive into that question here, but you can see how the line of research inquiry might progress as you start to explore your data with a machine learning model.

As a quick statistical note: above, we selected our variable importance method to be "permutation" for our demonstrative example. There are other options available in the `caret` package if you would like to explore those in your analyses.

### Comparing a random forest to a regression

You may be curious about comparing the predictive accuracy of the model to a linear model (a regression). Below, we'll specify a linear model and check out how the linear model performs in terms of predicting the real outcomes. We'll compare this with the random forest model's performance (rf_fit2). Note that we are not actually  re-running our random forest model here, but instead we are just making a dataset that includes the values that the rf_fit2 model predicted as well as the actual rf_fit2 values.


```r
#Make sure all variables stored as characters are converted to factors
data_train_lm <- 
    data_train %>% 
    mutate_if(is.character, as.factor) 
```

```
## mutate_if: no changes
```

```r
#Create a linear regression model, using the same formula approach as in the random forest: ~ .
lm_fit <-
    train(final_grade ~ .,
          data = data_train_lm,
          method = "lm")
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient fit
## may be misleading

## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient fit
## may be misleading

## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient fit
## may be misleading

## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient fit
## may be misleading
```

```r
#Append the predicted values to the training dataset for the linear model, so we can see both the predicted and the actual values
data_train_lm <-
    data_train %>%
    mutate(obs = final_grade,
           pred = predict(lm_fit, data_train_lm))
```

```
## mutate: new variable 'obs' with 273 unique values and 0% NA
```

```
##         new variable 'pred' with 398 unique values and 0% NA
```

```r
#Append the predicted values to the training dataset for the random forest
data_train_randomfor <-
    data_train %>%
    mutate(pred = predict(rf_fit2, data_train),
           obs = final_grade)
```

```
## mutate: new variable 'pred' with 398 unique values and 0% NA
```

```
##         new variable 'obs' with 273 unique values and 0% NA
```

```r
#Summarize, as data frames, the training data with the predicted and the actual values
    #for both the linear model
defaultSummary(as.data.frame(data_train_lm))
```

```
##       RMSE   Rsquared        MAE 
## 14.3559662  0.5879137 10.8437098
```

```r
    #and the random forest
defaultSummary(as.data.frame(data_train_randomfor))
```

```
##      RMSE  Rsquared       MAE 
## 4.6050924 0.9699852 3.3104543
```
Our output will come in the order we wrote the code, so the linear model output shows up above the random forest output.

We can see that the random forest technique seems to perform better than regression. Specifically, the RMSE is lower for the random forest (4.80 as compared to 13.50 for the linear model). Second, the variance explained (`Rsquared`) is much higher in the random forest (0.96 as compared to 0.59 for the linear model). It may be interesting to compare the results from the random forest not to a more straightforward model, such as a regression, but to a more sophisticated model, like one for deep learning. As you expand your skills, you might be curious to do something like that!
