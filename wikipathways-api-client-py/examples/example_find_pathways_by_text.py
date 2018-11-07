from wikipathways_api_client import WikipathwaysApiClient

def main():
    wikipathways_api_client_instance = WikipathwaysApiClient()

    # findPathwaysByText
    kwargs = {
        'query': 'apoptosis',
        'organism': 'http://identifiers.org/taxonomy/9606'
    }
    pathways_by_text = wikipathways_api_client_instance.find_pathways_by_text(**kwargs)
    print pathways_by_text 

if __name__ == '__main__':
    main()

