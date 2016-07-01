# mirna-pathway-finder
Find pathways with [microRNAs](http://www.ncbi.nlm.nih.gov/pubmed/14744438) (miRNAs) at WikiPathways.

miRNAs target genes to regulate gene expression.  As a component for the Genboree workbench, this library:
1. takes as input a user-specified list of miRNAs and
2. produces as output a list of pathways that contain
   A. one or more miRNAs from the user-specified list, and/or 
   B. one or more genes that are targeted by one or more of the miRNAs from the user-specified list

miRNA to gene mapping is performed by pre-processing, as described [here](https://github.com/nrnb/mirna-pathway-finder/blob/master/wp-mir-table-builder/wp-mir-table-hs-readme.txt).

The output table looks like this:

| *Pathway Title* | *Linkout* | *miRNAs* | *miRNA Targets* |
| --------------- | --------- | -------- | --------------- |
| Sample Pathway | [WP4](http://www.wikipathways.org/wpi/WP4) | 2 | 3 |
| Sample Pathway | [WP4](http://www.wikipathways.org/wpi/WP4) | 2 | 3 |
| Sample Pathway | [WP4](http://www.wikipathways.org/wpi/WP4) | 2 | 3 |

## To install

```
git clone git@github.com:ariutta/mirna-pathway-finder.git
cd mirna-pathway-finder
```

Set up your virtualenv:

(how to install [virtualenv on Mac](http://exponential.io/blog/2015/02/10/install-virtualenv-and-virtualenvwrapper-on-mac-os-x/)) 

```
mkvirtualenv mirna-pathway-finder
```

or if you've already created this virtualenv:

```
workon mirna-pathway-finder
```

Install dependencies:

```
pip install -e .
```

## To run

Check command line argument defaults:

```
python mirnapathwayfinder/__init__.py -h
```

if the defaults are OK, you can then find pathways with miRNAs for a specific node_id:

```
python mirnapathwayfinder/__init__.py hsa-miR-370-3p
```

or provide a file path to a list of node_ids:

```
python mirnapathwayfinder/__init__.py 'tests/test1/input/node-list.txt'
```

Or you can override the defaults, e.g., override default output directory path:

```
python mirnapathwayfinder/__init__.py 'tests/test1/input/node-list.txt' -o 'tests/test1/output-actual/'
```

## Todo
* [x] Get command line arguments working
* [ ] Test
* [ ] Publish
