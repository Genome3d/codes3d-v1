from wikipathways_api_client import WikipathwaysApiClient

def main():
    wikipathways_api_client_instance = WikipathwaysApiClient()

    # findPathwaysByXref
    kwargs = {
        'system_codes': ['X', 'L'],
        'identifiers': ['201746_at', '7040']
    }
    pathways_by_xref = wikipathways_api_client_instance.find_pathways_by_xref(**kwargs)
    print pathways_by_xref 

if __name__ == '__main__':
    main()
