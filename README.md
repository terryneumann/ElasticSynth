# ElasticSynth
### Synthetic Controls for Comparitive Case Studies Using Elastic Net Regression (glmnet) for Unit Selection

This is an implementation of a more 'flexible' version of the synthetic controls method for case studies and causal impact evaluation.

The traditional synthetic controls implemented by Abadie, Diamond, and Hainmueller (2010) has 5 constraints placed on the output weights:

1. No Intercept 
2. Weights Sum To One
3. Non-Negativity
4. Exact-Balance
5. Constant Weights

However, it is not necessary in all cases to hold all of these constraints.

This method, published by Doudchenko and Imbens in 2016, is an implementation of Synthetic Controls that only holds two of these constraints:

1. No Intercept
2. Constant Weights

by using <b>glmnet</b> as the underlying optimization.

Install with the following command
```r

devtools::install_github('tdn158/ElasticSynth')


```





