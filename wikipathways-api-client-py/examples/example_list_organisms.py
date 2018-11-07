from wikipathways_api_client import WikipathwaysApiClient

def main():
    wikipathways_api_client_instance = WikipathwaysApiClient()

    # Get organisms
    organisms = wikipathways_api_client_instance.list_organisms()
    print organisms

if __name__ == '__main__':
    main()

