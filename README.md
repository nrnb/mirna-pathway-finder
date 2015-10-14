# mirna-pathway-finder
Find pathways with miRNAs at WikiPathways.

A component for the Genboree workbench.

## To install

```
git clone git@github.com:ariutta/mirna-pathway-finder.git
cd mirna-pathway-finder
```

Set up your virtualenv:

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
