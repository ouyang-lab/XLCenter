# XLCenter
Precise detection of cross-linking induced termination sites

## Software Pre-requisites
  * Python (>= v 3.13)
  * [numpy] (>= v 1.15.2)
  * [scipy] (>= v 1.0.0)
  * [pandas] (>= v 0.24.2)
  * [statsmodels] (>= v 0.10.1)
  * [scikit-learn] (>= v 0.20.4)
  
  [numpy]: https://pypi.org/project/numpy/
  [scipy]: https://pypi.org/project/scipy/
  [pandas]: https://pypi.org/project/pandas/
  [statsmodels]: https://pypi.org/project/statsmodels/
  [scikit-learn]: https://pypi.org/project/scikit-learn/

## Installation

### Unzip the gzipped tarball package

```
tar zxvf xlcenter_v1.0.tar.gz
```

## Examples
RBFOX2 eCLIP chromosome 1 bam file in the examples folder

## Usage

To call crosslinking-induced sites
```
python xlcenter.py callsites
    -i examples/RBFOX2_eCLIP_rep1_chr1.bam
    -o output/RBFOX2_eCLIP_rep1_chr1
    -f examples/ucsc.hg19.sizes
    -e eCLIP
    -t 0.9
    -d bum
```

To use a more stringent threshold
```
python xlcenter.py threshold
    -p output/RBFOX2_eCLIP_rep1_chr1.txt
    -b output/RBFOX2_eCLIP_rep1_chr1_0.999.bed
    -t 0.999
    -d beta
```
