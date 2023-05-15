# Second attempt:

## Ideas

poster CROI2020-PVL-poster-200225 rg-1.pptx,
estimating mean population viral load by age and sex, in different communities and by round of collection.
The idea is to extend Imogen's analysis to multiple round, with a particular focus on age-shifts (changes in the drivers-profile of the epidemic)

More aims in the [google doc](https://docs.google.com/document/d/1w9OYMA2JdW_tzH96oaJgVo3BZTnNb3OmPLuRidxxe-A/edit#)

Possible extensions:
- it would be great if we were able to link PVL to incidence.
- 'normalise' PVL by sexual contact intensities.


## Code

Code shared by Oli at:
- `Phyloscanner.R.utilities/misc/phyloscan.viral.load.project.R`

## Lit review:

- https://academic.oup.com/cid/article/66/8/1254/4662847?login=true
- https://www.nejm.org/doi/full/10.1056/nejmoa1702150



## Posterior predictive check:
From AkiVethari's [github](https://avehtari.github.io/modelselection/)

- Simulate new datasets with posterior distribution and then check whether our data is consistent with that.
- (Need a way to visualise then)
- Same but with ancillary statistics. (ancillary meaning independent of model?)
- Can be used to have a posterior predictive p-value. (ppp-value) 
- but usully just better to compare against distribution.

we are using double use of data, so risks of p-value to not be calibrated. For that consider CV.
(is there an easy way to perform CV in STAN?)




