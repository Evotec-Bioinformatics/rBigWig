# rBigWig

A R package for working with bigWig files using the [libBigWig C-library](https://github.com/dpryan79/libBigWig)

# Usage

```R
# The package functions always require the path to the file
filename_bigwig <- system.file("data/test.bw", package="rBigWig")

# Fetch the individual scores in the region
rBigWig::fetch_region(filename_bigwig, chromosome="1", start=50, end=150)
#     Position Score
# 1         50     1
# 2         51     1
# 3         52     1
# ...
# 50        99     1
# 51       100     2
# 52       101     2
# ...
# 100      149     2

# Fetch summary statistics 
rBigWig::fetch_region_means(filename_bigwig, chromosome="1", start=50, end=150, bins=10, as_data_frame=TRUE)
#    Start End Score
# 1     50  60     1
# 2     60  70     1
# 3     70  80     1
# 4     80  90     1
# 5     90 100     1
# 6    100 110     2
# 7    110 120     2
# 8    120 130     2
# 9    130 140     2
# 10   140 150     2

``` 


# Installation

Please make sure to have gcc and make installed. Then simply install the package from Github using [devtools](https://github.com/r-lib/devtools):

```R
devtools::install_github("Evotec-Bioinformatics/rBigWig")
```

# Authors

This package was developed by [Evotec SE, Bioinformatics and Biostatistics](https://github.com/Evotec-Bioinformatics)

- [Manuel Landesfeind](https://github.com/landesfeind)

