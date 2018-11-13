# Python API client for WikiPathways webservices.

## How to Use

Clone and `cd` into this repo:

```
git clone https://github.com/wikipathways/wikipathways-api-client-py.git
cd wikipathways-api-client-py
```

Install [virtualenvwrapper](https://virtualenvwrapper.readthedocs.org/en/latest/install.html) and create a new virtualenv (if `pip` is available on your system, the virtualenvwrapper steps are optional):

```bash
mkvirtualenv your_preferred_name
```

Install `wikipathways-api-client-py`:

```
pip install -e .
```

Try one of the examples:

```
python ./examples/example_get_pathway_as.py
```

Now that your virtualenv is created, you can leave it at anytime with:

```bash
deactivate
```

And come back to it later with:

```bash
workon your_preferred_name
```

To list all your virtualenvs:

```bash
ls $WORKON_HOME
```

## Troubleshooting

This set of commands was required for one user to install on Linux:

```
sudo apt-get update
sudo apt-get install aptitude
sudo aptitude install libxslt1-dev libxml2-dev
git clone https://github.com/wikipathways/wikipathways-api-client-py.git
cd wikipathways-api-client-py
pip install -e .
python wikipathways-api-client-py/examples/example_get_pathway_as.py
```

## Supported API calls

If it's checked in the following list, the WikiPathways API call is implemented in this client. If it's not checked, we'd be happy to accept your pull request with code implementing it! You can find more details about the calls in the [WikiPathways API documentation](http://wikipathways.org/index.php/Help:WikiPathways_Webservice/API).

- [x] listOrganisms
- [x] listPathways
- [x] getPathwayInfo
- [ ] getPathwayHistory
- [ ] getRecentChanges
- [ ] login
- [x] getPathwayAs
- [ ] updatePathway
- [x] createPathway
- [x] findPathwaysByText
- [x] findPathwaysByXref
- [ ] findInteractions
- [ ] saveCurationTag
- [ ] removeCurationTag
- [ ] getCurationTags
- [ ] getCurationTagsByName
- [x] getColoredPathway
- [ ] getXrefList
- [ ] findPathwaysByLiterature
- [ ] getOntologyTermsByPathway
- [ ] getOntologyTermsByOntology
- [ ] getPathwaysByOntologyTerm
- [ ] getPathwaysByParentOntologyTerm

## Future Enhancements
- [ ] Add tests
- [ ] Add documentation and build docs with PyDoc
- [ ] Publish on PyPI so it can be installed with ```pip```

