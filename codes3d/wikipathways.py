#!/usr/bin/env python

from wikipathways_api_client import WikipathwaysApiClient


def pathways(gene):
    client = WikipathwaysApiClient()
    kwargs = {'query': gene,
              'organism': 'http://identifiers.org/taxonomy/9606'
             }
    return set((pathway['identifier'], pathway['name']) for pathway
                 in client.find_pathways_by_text(**kwargs))

def shared_pathways(genes):
    client = WikipathwaysApiClient()
    kwargs = {'query': ' AND '.join(genes),
              'organism': 'http://identifiers.org/taxonomy/9606'
             }
    return set((pathway['identifier'], pathway['name']) for pathway
                 in client.find_pathways_by_text(**kwargs))
