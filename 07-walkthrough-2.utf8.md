---
title: 'Education Dataset Analysis Pipeline: Walkthrough #2'
output: html_document
---
## Introduction
Gradebooks are nearly ubiquitous throughout K-12 classrooms, whether they exist as standalone Excel files, Google Sheets, or in proprietary software. 

This walkthrough goes through a series of analyses using the data science framework (link), using the sample [Assessment Types - Points](http://web.mit.edu/jabbott/www/excelgradetracker.html) Excel gradebook template from MIT. All data in the sample gradebook have been generated, and do not reflect individual student data.

*** 

## Driving Question and Objectives

***

## Data Import
Setting up our environment (note: how deep do we go into working directories?!)


Importing our data (need to sim data for 25 students)
Check text for object naming conventions, discussion of .csv, .xlsx, versatility of import functions within the `tidyverse`
File naming - issues that can arise from spaces

```r
gradebook <- readxl::read_excel(here("/data/gradebooks", "ExcelGradeTrackerAssessmentTypePoints_SIMDATA_01.xlsx"))
```

```
## New names:
## * `` -> ...1
## * `` -> ...3
## * `` -> ...4
## * `` -> ...5
## * `` -> ...6
## * … and 39 more problems
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```r
summary(cars)
```

```
##      speed           dist       
##  Min.   : 4.0   Min.   :  2.00  
##  1st Qu.:12.0   1st Qu.: 26.00  
##  Median :15.0   Median : 36.00  
##  Mean   :15.4   Mean   : 42.98  
##  3rd Qu.:19.0   3rd Qu.: 56.00  
##  Max.   :25.0   Max.   :120.00
```

## Including Plots

You can also embed plots, for example:

<img src="07-walkthrough-2_files/figure-html/pressure-1.png" width="672" />

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
