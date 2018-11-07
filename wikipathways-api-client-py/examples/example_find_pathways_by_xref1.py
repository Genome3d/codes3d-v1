from wikipathways_api_client import WikipathwaysApiClient

def main():
    wikipathways_api_client_instance = WikipathwaysApiClient()

    # findPathwaysByXref
    kwargs = {
        '@id': [
            'http://identifiers.org/ncbigene/7040',
            'http://identifiers.org/affy.probeset/201746_at'
        ]
    }
    pathways_by_xref = wikipathways_api_client_instance.find_pathways_by_xref(**kwargs)
    print pathways_by_xref 

if __name__ == '__main__':
    main()

