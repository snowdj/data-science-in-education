
# Walkthrough 6: Exploring Relationships Using Social Network Analysis With Social Media Data {#c12}

## Vocabulary

- social network analysis
- Application Programming Interface (API)
- edgelist
- edge
- vertex
- sociogram
- selection model
- influence model

## Introduction

This chapter builds upon the previous chapter on [text analysis](#c11) of #tidytuesday data.

### Background

In the past, if a teacher wanted advice about how to plan a unit or to design a lesson, they would likely turn to a trusted peer in their building or district [@spillane2012]. In the present, though, they are as likely to turn to someone in the professional learning network [@trust2016].

There are a few reasons to be interested in social media. For example, if you work in a school district, you may be interested in who is interacting with the content you share. If you are a researcher, you may wish to investigate what teachers, administrators, and others do through state-based hashtags (e.g., @rosenberg2016). Social media-based data can also be interesting because it provides new contexts for learning to take place, such as learning through informal communities.

In this chapter, we focus on a source of data that could be used to understand how one new community functions. That community, #tidytuesday is one sparked by the work of one of the *Data Science in Education Using R* co-authors, Jesse Mostipak, who created the #r4ds community from which #tidytuesday was created. #tidytuesday is a weekly data visualization challenge. A great place to see examples from past #tidytuesday challenges is an interactive Shiny application (https://github.com/nsgrantham/tidytuesdayrocks). 

In this walkthrough, we focus on a) accessing data on #tidytuesday from Twitter and b) trying to understand the nature of the interactions that take place through #tidytuesday. We note that while we focused on #tidytuesday because we think it exemplifies the new kinds of learning that a data science toolkit allows an analyst to try to understand (through new data sources), we also chose this because it is straightforward to access data from Twitter, and we think you may find other opportunities to analyze data from Twitter in other cases.

### Packages, Data Sources and Import, and Methods

In this chapter, we access data using the rtweet package [@kearney2016]. Through rtweet, it is easy to access data from Twitter as long as one has a Twitter account. We will load the tidyverse and rtweet packages to get started. Here is an example of searching the most recent 1,000 tweets which include the hashtag #rstats. When you run this code, you will be prompted to authenticate your access via Twitter. We will also load other packages that we will be using in this analysis.


```r
library(tidyverse)
library(rtweet)
library(dataedu)
library(randomNames)
library(tidygraph)
library(ggraph)
```


```r
rstats_tweets <- 
  search_tweets("#rstats")
```

As described in [the previous chapter](#c11), the search term can be easily changed to other hashtags - or other terms. To search for #tidytuesday tweets, we can simply replace #rstats with #tidytuesday. 


```r
tidytuesday_tweets <- 
  search_tweets("#tidytuesday")
```

## View Data

We can see that there are *many* rows for the data:


```r
tt_tweets
```

```
#> # A tibble: 4,418 x 90
#>    user_id status_id created_at          screen_name text  source
#>    <chr>   <chr>     <dttm>              <chr>       <chr> <chr> 
#>  1 115921… 11631542… 2019-08-18 18:22:42 MKumarYYC   "Fir… Twitt…
#>  2 107332… 11632475… 2019-08-19 00:33:11 cizzart     "El … Twitt…
#>  3 107332… 11450435… 2019-06-29 18:57:17 cizzart     "Pro… Twitt…
#>  4 107332… 11168648… 2019-04-13 00:45:15 cizzart     "#Ar… Twitt…
#>  5 107332… 11228824… 2019-04-29 15:17:02 cizzart     "Pes… Twitt…
#>  6 107332… 11176387… 2019-04-15 04:00:17 cizzart     "Dat… Twitt…
#>  7 107332… 11245531… 2019-05-04 05:55:32 cizzart     "El … Twitt…
#>  8 107332… 11407021… 2019-06-17 19:25:50 cizzart     "#da… Twitt…
#>  9 107332… 11325299… 2019-05-26 06:12:46 cizzart     "El … Twitt…
#> 10 107332… 11233585… 2019-04-30 22:48:43 cizzart     "Vis… Twitt…
#> # … with 4,408 more rows, and 84 more variables: display_text_width <dbl>,
#> #   reply_to_status_id <chr>, reply_to_user_id <chr>,
#> #   reply_to_screen_name <chr>, is_quote <lgl>, is_retweet <lgl>,
#> #   favorite_count <int>, retweet_count <int>, quote_count <int>,
#> #   reply_count <int>, hashtags <list>, symbols <list>, urls_url <list>,
#> #   urls_t.co <list>, urls_expanded_url <list>, media_url <list>,
#> #   media_t.co <list>, media_expanded_url <list>, media_type <list>,
#> #   ext_media_url <list>, ext_media_t.co <list>, ext_media_expanded_url <list>,
#> #   ext_media_type <chr>, mentions_user_id <list>, mentions_screen_name <list>,
#> #   lang <chr>, quoted_status_id <chr>, quoted_text <chr>,
#> #   quoted_created_at <dttm>, quoted_source <chr>, quoted_favorite_count <int>,
#> #   quoted_retweet_count <int>, quoted_user_id <chr>, quoted_screen_name <chr>,
#> #   quoted_name <chr>, quoted_followers_count <int>,
#> #   quoted_friends_count <int>, quoted_statuses_count <int>,
#> #   quoted_location <chr>, quoted_description <chr>, quoted_verified <lgl>,
#> #   retweet_status_id <chr>, retweet_text <chr>, retweet_created_at <dttm>,
#> #   retweet_source <chr>, retweet_favorite_count <int>,
#> #   retweet_retweet_count <int>, retweet_user_id <chr>,
#> #   retweet_screen_name <chr>, retweet_name <chr>,
#> #   retweet_followers_count <int>, retweet_friends_count <int>,
#> #   retweet_statuses_count <int>, retweet_location <chr>,
#> #   retweet_description <chr>, retweet_verified <lgl>, place_url <chr>,
#> #   place_name <chr>, place_full_name <chr>, place_type <chr>, country <chr>,
#> #   country_code <chr>, geo_coords <list>, coords_coords <list>,
#> #   bbox_coords <list>, status_url <chr>, name <chr>, location <chr>,
#> #   description <chr>, url <chr>, protected <lgl>, followers_count <int>,
#> #   friends_count <int>, listed_count <int>, statuses_count <int>,
#> #   favourites_count <int>, account_created_at <dttm>, verified <lgl>,
#> #   profile_url <chr>, profile_expanded_url <chr>, account_lang <lgl>,
#> #   profile_banner_url <chr>, profile_background_url <chr>,
#> #   profile_image_url <chr>
```

## Process Data

Network data, in general, and network data from Twitter, particularly, requires some processing before it can be used in subsequent analyses. In particular, we are going to create an edgelist, a data structure that is especially helpful for understanding the nature of relationships. 

An edgelist looks like the following, where the sender denotes who is initiating the interaction or relationship, and the receiver is who is the recipient of it:




```
#> # A tibble: 12 x 2
#>    sender              receiver          
#>    <chr>               <chr>             
#>  1 Rodriguez, Danielle Rea, Ramona       
#>  2 el-Bacchus, Sameera Ho, Aaron         
#>  3 el-Bacchus, Sameera Kowalik, Jalynn   
#>  4 Victoria, Kendall   Ho, Aaron         
#>  5 Victoria, Kendall   Rea, Ramona       
#>  6 Victoria, Kendall   al-Ghazal, Badraan
#>  7 Netsanet, Monica    Kowalik, Jalynn   
#>  8 Netsanet, Monica    Martinez, Julio   
#>  9 Netsanet, Monica    al-Ghazal, Badraan
#> 10 al-Youssef, Uqbah   Jimenez, Leahana  
#> 11 Hamill, Lauren      Kowalik, Jalynn   
#> 12 Hamill, Lauren      Jimenez, Leahana
```

In this edgelist, the sender could indicate, for example, someone who nominates someone else (the receiver) as someone they go to for help. The sender could also indicate someone who interacted with the receiver, such as by recognizing one of their tweets with a favorite (or a mention). In the following steps, we will work to create an edgelist from the data from #tidytuesday on Twitter.

### Extracting mentions

Let's extract the mentions. There is a lot going on in the code below; let's break it down line-by-line, starting with the `mutate()`:

- `mutate(all_mentions = str_extract_all(text, regex))`: this line uses a regex, or regular expression, to identify all of the usernames in the tweet (*note*: the regex comes from from [this page](https://stackoverflow.com/questions/18164839/get-twitter-username-with-regex-in-r))
- `unnest(all_mentions)` this line uses a tidyr function, `unnest()` to move every mention to its own line, while keeping all of the other information the same (see more about `unnest()` here: https://tidyr.tidyverse.org/reference/unnest.html)
 

```r
regex <- "@([A-Za-z]+[A-Za-z0-9_]+)(?![A-Za-z0-9_]*\\.)"

tt_tweets <-
  tt_tweets %>%
  # Use regular expression to identify all the usernames in a tweet
  mutate(all_mentions = str_extract_all(text, regex)) %>%
  unnest(all_mentions)
```

Let's put these into their own data frame, called `mentions`.


```r
mentions <-
  tt_tweets %>%
  mutate(all_mentions = str_trim(all_mentions)) %>%
  select(sender = screen_name, all_mentions)
```

### Putting the edgelist together

An edgelist is a common social network analysis data structure that has columns for the "sender" and "receiver" of interactions, or relations. For example, someone "sends" the mention to someone who is mentioned, who can be considered to "receive" it. This will require one last processing step. Let's look at our data as it is now.


```r
mentions
```

```
#> # A tibble: 2,447 x 2
#>    sender  all_mentions    
#>    <chr>   <chr>           
#>  1 cizzart @eldestapeweb   
#>  2 cizzart @INDECArgentina 
#>  3 cizzart @ENACOMArgentina
#>  4 cizzart @tribunalelecmns
#>  5 cizzart @CamaraElectoral
#>  6 cizzart @INDECArgentina 
#>  7 cizzart @tribunalelecmns
#>  8 cizzart @CamaraElectoral
#>  9 cizzart @AgroMnes       
#> 10 cizzart @AgroindustriaAR
#> # … with 2,437 more rows
```

What needs to happen to these to make them easier to work with in an edgelist? One step is to remove the "@" symbol from the columns we created and to save the results to a new tibble, `edgelist`.


```r
edgelist <- 
  mentions %>% 
  # remove "@" from all_mentions column
  mutate(all_mentions = str_sub(all_mentions, start = 2)) %>% 
  # rename all_mentions to receiver
  select(sender, receiver = all_mentions)
```

## Analysis and Results

Now that we have our edgelist, it is straightforward to plot the network. We'll use the {tidygraph} and {ggraph} packages to visualize the data.

### Plotting the network

Because large networks (like this one) can present challenges, it is common to filter them to only include some individuals. Let's explore how many interactions each individual in the network sent.


```r
interactions_sent <- edgelist %>% 
  # this counts how many times each sender appears in the data frame, effectively counting how many interactions each individual sent 
  count(sender) %>% 
  # arranges the data frame in descending order of the number of interactions sent
  arrange(desc(n))

interactions_sent
```

```
#> # A tibble: 618 x 2
#>    sender            n
#>    <chr>         <int>
#>  1 thomas_mock     347
#>  2 R4DScommunity    78
#>  3 WireMonkey       52
#>  4 CedScherer       41
#>  5 allison_horst    37
#>  6 mjhendrickson    34
#>  7 kigtembu         27
#>  8 WeAreRLadies     25
#>  9 PBecciu          23
#> 10 sil_aarts        23
#> # … with 608 more rows
```

618 senders of interactions is a lot! What if we focused on only those who sent more than one interaction?


```r
interactions_sent <- interactions_sent %>% 
  filter(n > 1)
```

That leaves us with only 349, perhaps a more reasonable number.

We now need to filter the edgelist to only include these 349 individuals. The following code simply uses the `filter()` function combined with the `%in%` operator to do this:


```r
edgelist <- edgelist %>% 
  # the first of the two lines below filters to include only senders in the interactions_sent data frame
  # the second line does the same, for receivers
  filter(sender %in% interactions_sent$sender,
         receiver %in% interactions_sent$sender)
```

We'll use the `as_tbl_graph()` function, which (by default) identified the first column as the "sender" and the second as the "receiver." Let's look at the object it creates, too.


```r
g <- 
  as_tbl_graph(edgelist)

g
```

```
#> # A tbl_graph: 267 nodes and 975 edges
#> #
#> # A directed multigraph with 7 components
#> #
#> # Node Data: 267 x 1 (active)
#>   name           
#>   <chr>          
#> 1 dgwinfred      
#> 2 datawookie     
#> 3 jvaghela4      
#> 4 FournierJohanie
#> 5 JonTheGeek     
#> 6 jakekaupp      
#> # … with 261 more rows
#> #
#> # Edge Data: 975 x 2
#>    from    to
#>   <int> <int>
#> 1     1    32
#> 2     1    36
#> 3     2   120
#> # … with 972 more rows
```

We can see that the network now consists of 267 individuals - the 267 who sent more than one interaction.

Next, we'll use the `ggraph()` function:


```r
g %>%
  # we chose the kk layout as it created a graph which was easy-to-interpret, but others are available; see ?ggraph
  ggraph(layout = "kk") +
  # this adds the points to the graph
  geom_node_point() +
  # this adds the links, or the edges; alpha = .2 makes it so that the lines are partially transparent
  geom_edge_link(alpha = .2) +
  # this last line of code adds a ggplot2 theme suitable for network graphs
  theme_graph()
```

<img src="12-wt-social-network-analysis_files/figure-html/unnamed-chunk-16-1.png" width="100%" style="display: block; margin: auto;" />

Finally, let's size the points based on a measure of centrality, typically a measure of how (potentially) influence an individual may be, based on the interactions observed.




```r
g %>% 
  # this calculates the centrality of each individual using the built-in centrality_authority() function
  mutate(centrality = centrality_authority()) %>% 
  ggraph(layout = "kk") + 
  geom_node_point(aes(size = centrality, color = centrality)) +
  # this line colors the points based upon their centrality
  scale_color_continuous(guide = 'legend') + 
  geom_edge_link(alpha = .2) +
  theme_graph()
```

<img src="12-wt-social-network-analysis_files/figure-html/unnamed-chunk-18-1.png" width="100%" style="display: block; margin: auto;" />

There is much more you can do with {ggraph} (and {tidygraph}); check out the {ggraph} tutorial here: https://ggraph.data-imaginist.com/

## Conclusion

In this chapter, we used social media data (from the #tidytuesday hashtag) to prepare and visualize social network data. This is a powerful technique; one that can reveal who is interacting with whom, and one that can begin to suggest why.

<!-- Wondering if we need to do more with the visualizations do begin so that this can be more warrented. -->

Behind these visualizations, though, there are also statistical models and methods that can help to understand what is going on with respect to particular relationships in a network in additional ways.

One way to consider these models and methods is in terms of two *processes* at play in our relationships (cite). These two processes are commonly (though not exclusively) the focus of statistical analyses of networks. In addition to not being exclusive, they do not interact independently: they affect each other reciprocally (Xu, Frank, & Penuel, 2018). They are:

- *Selection*: the processes regarding who chooses to have a relationship with whom
- *Infuence*: the processes regarding how who we have relationships with affects our behavior

While these are complex, they can be studied with the type of data collected from asking people about their relationships (and possibly asking them about or studying their behavior--or measuring some outcome). Happily, the use of these methods has expanded along with R: many of the best tools for studying social networks are in the form of long-standing R packages. Additionally, while there are many potential nuances to studying selection and influence, these are models that can fundamentally be carried out with regression, or the linear model (or extensions of it). We describe these in the *Technical Appendix* for this chapter, as they do not use the tidytuesday dataset and are likely to be of interest to readers only after having mastered preparing and visualizing network data.

## Technical Appendix: Influence and selection models

As noted above, there is much more to understanding interactions, and network analysis, beyond creating edgelists and visualizing network data (through the use of an edgelist). Two processes that are particularly important (and able to be studied with network data using R) are for influence and selection. 

### An example of influence

First, let's look at an example of influence. To do so, let's create three different data frames. Here is what they should, at the end of the process, contain:

- A data frame indicating who the *nominator* and *nominee* for the relation (i.e., if Stefanie says that José is her friend, then Stefanie is the nominator and José the nominee) - as well as an optional variable indicating the weight, or strength, of their relation.
- This data frame and its type can be considered the basis for many types of social network analysis and is a common structure for network data: it is an *edgelist*.
- Data frames indicating the values of some behavior - an outcome - at two different time points.

In this example, we create some example data that can be used to explore questions about how influence works.

Let's take a look at the merged data. What this data now contains is the first data frame, `data1`, with each nominees' outcome at time 1 (`yvar1`). Note that we will find each nominators' outcome at time 2 later on.


```r
data1 <-
  data.frame(
    nominator = c(2, 1, 3, 1, 2, 6, 3, 5, 6, 4, 3, 4),
    nominee = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 6),
    relate = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  )

data2 <-
  data.frame(nominee = c(1, 2, 3, 4, 5, 6),
             yvar1 = c(2.4, 2.6, 1.1, -0.5, -3, -1))

data3 <-
  data.frame(nominator = c(1, 2, 3, 4, 5, 6),
             yvar2 = c(2, 2, 1, -0.5, -2, -0.5))
```

### Joining the data

Next, we'll join the data into one data frame. Note that while this is sometimes tedius and time-consuming, especially with large sources of network data, it is a key step for being able to carry out network analysis - often, even for creating visualiations that are informative.


```r
data <-
  left_join(data1, data2, by = "nominee")

data <-
  data %>% 
  # this makes merging later easier
  mutate(nominee = as.character(nominee)) 

# calculate indegree in tempdata and merge with data
tempdata <- data.frame(table(data$nominee))

tempdata <-
  tempdata %>%
  rename(
    # rename the column "Var1" to "nominee" 
    "nominee" = "Var1", 
    # rename the column "Freq" to "indegree"
    "indegree" = "Freq"
    ) %>% 
  # makes nominee a character data type, instead of a factor, which can cause problems
  mutate(nominee = as.character(nominee))

data <- 
  left_join(data, tempdata, by = "nominee")
```

#### Calculating an exposure term

This is the key step that makes this model - a regression, or linear, model - one that is special. It is creating an exposure term. The idea is that the exposure term "captures" how your interactions with someone, over some period of time (between the first and second time points) impact some outcome. This model accounts for an individual's initial report of the outcome, i.e., their time 1 prior value, so it is a model for *change* in some outcome.


```r
# Calculating exposure
data <-
  data %>% 
  mutate(exposure = relate * yvar1)

# Calculating mean exposure
mean_exposure <-
  data %>%
  group_by(nominator) %>%
  summarize(exposure_mean = mean(exposure))
```

What this data frame - `mean_exposure` - contains is the mean of the outcome (in this case, `yvar1`) for all of the individuals the nominator had a relation with.

As we need a final data set with `mean_exposure`,`degree`, `yvar1`, and `yvar2` added, we'll process the data a bit more.


```r
data2 <-
  data2 %>% 
  # rename nominee as nominator to merge these
  rename("nominator" = "nominee") 

final_data <-
  left_join(mean_exposure, data2, by = "nominator")

final_data <- 
  # data3 already has nominator, so no need to change
  left_join(final_data, data3, by = "nominator") 
```

#### Regression (linear model)

Calculating the exposure term is the most distinctive and important step in carrying out influence models. Now, we can simply use a linear model to find out how much relations - as captured by the influence term - affect some outcome.


```r
model1 <-
  lm(yvar2 ~ yvar1 + exposure_mean, data = final_data)

summary(model1)
```

```
#> 
#> Call:
#> lm(formula = yvar2 ~ yvar1 + exposure_mean, data = final_data)
#> 
#> Residuals:
#>       1       2       3       4       5       6 
#>  0.0295 -0.0932  0.0943 -0.0273 -0.0255  0.0222 
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)     0.1161     0.0345    3.37    0.043 *  
#> yvar1           0.6760     0.0241   28.09  9.9e-05 ***
#> exposure_mean   0.1254     0.0361    3.47    0.040 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.0823 on 3 degrees of freedom
#> Multiple R-squared:  0.998,	Adjusted R-squared:  0.997 
#> F-statistic:  945 on 2 and 3 DF,  p-value: 6.31e-05
```

So, the influence model is used to study a key process for social network analysis, but it is one that is useful, because you can quantify, given what you measure and how you measure it, *the network effect*, something that is sometimes not considered, especially in education (Frank, 2009). It's also fundamentally a regression. That's really it, as the majority of the work goes into calculating the exposure term.

### An example of selection

Selection models are also commonly used - and are commonly of interest not only to researchers but also to administrators and teachers (and even to youth and students). 

Here, we briefly describe a few possible approaches for using a selection model.

At its core, the selection model is a regression - albeit, one that is a generalization of one, namely, a logistic regression (sometimes termed a generalized linear model, because it is *basically* a regression but is one with an outcome that consists just of 0's and 1's). Thus, the most straight-away way to use a selection model is to use a logistic regression where all of the relations (note the `relate` variable in `data1` above) are indicated with a 1. But, here is the important and challenging step: all of the *possible relations* (i.e., all of the relations that are possible between all of the individuals in a network) are indicated with a 0 in an edgelist. Note that, again, an edgelist is the preferred data structure for carrying out this analysis. This step involves some data wrangling, especially the idea of widening or lengthening a data frame.

Once all of the relations are indicated with a 1 or a 0, then a simple linear regression can be used. Imagine that we are interested in whether individuals from the *same* group are more or less likely to interact than those from different groups; same could be created in the data frame based upon knowing which group both nominator and nominee are from:


```r
m_selection <- 
  glm(relate ~ 1 + same, data = edgelist1)
```

While this is a straightforward way to carry out a selection model, there are some limitations to it. Namely, it does not account for individuals who send more (or less) nominations overall--and not considering this may mean other effects, like the one associated with being from the *same* group, are not accurate. A few extensions of the linear model - including those that can use data for which relationships are indicated with weights, not just 1's and 0's, have been developed. 

One type of model extends the logistic regression. It can be used for data that is not only 1's and 0's but also data that is normally distributed . It is the amen package available [here](https://cran.r-project.org/web/packages/amen/index.html).

A particularly common one is an Exponential Random Graph Model, or an ERGM. An R package that makes estimating these easy is available [here](https://cran.r-project.org/web/packages/ergm/index.html). That R package, {ergm}, is part of a powerful and often-used collection of packages, including those for working with network data (data that can begin with an edgelist, but may need additional processing that is challenging to do with edgelist data), {statnet}. A link to the statnet packages is [here](https://statnet.org/).

Trust, T., Krutka, D. G., & Carpenter, J. P. (2016). “Together we are better”: Professional learning networks for teachers. Computers & education, 102, 15-34.
