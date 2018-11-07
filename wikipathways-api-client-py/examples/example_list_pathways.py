from wikipathways_api_client import WikipathwaysApiClient

def main():
    wikipathways_api_client_instance = WikipathwaysApiClient()

    pathways_all = wikipathways_api_client_instance.list_pathways(organism = 'Homo sapiens')
    print pathways_all

if __name__ == '__main__':
    main()

