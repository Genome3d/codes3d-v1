from wikipathways_api_client import WikipathwaysApiClient

def main():
    wikipathways_api_client_instance = WikipathwaysApiClient()

    # Get pathway info
    info = wikipathways_api_client_instance.get_pathway_info('WP274')
    print info

if __name__ == '__main__':
    main()

