# DunedinPACE-python

Python implementation of [DunedinPACE](https://elifesciences.org/articles/73420) pace of aging biomarker.
[Original R implementation](https://github.com/danbelsky/DunedinPACE).

## Background
DNA methylation beta values are a measure of the methylation level at each CpG site (regions of DNA where a cytosine nucleotide is followed by a guanine nucleotide). The beta value ranges from 0 to 1, representing a percentage. A value of 0 indicates no methylation at all, while a value of 1 indicates complete methylation. These values can be used to provide a quantitative measure of the level of DNA methylation at a specific CpG site, which can provide useful insights into gene expression, cellular differentiation, and disease states.

In the context of aging, the patterns of DNA methylation change over time. Initial epigenetic clocks (eg Horvath's clock) were able to predict chronological age from methylation levels. Subsequent clocks were trained to predict better measurements of aging/disease rather than  simply one's chronological age. However, these were created with samples from many people over different timepoints in their life incurring noise. DunedinPACE is recognized to have the best aging signal because it uses longitudinal data from the same patients across their lifespan. It is currently recognized as the most predictive clock with high correlations to disease and quality of life metrics. It's used as a measure of impact from age rejuvenation protocols such as Bryan Johnson's [blueprint](https://blueprint.bryanjohnson.co/) / [Rejuvenation Olympics](https://rejuvenationolympics.com/about-us/).

## Algorithm Implementation

```
print(1)
```


