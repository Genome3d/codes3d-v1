from wikipathways_api_client import WikipathwaysApiClient

def main():
    wikipathways_api_client_instance = WikipathwaysApiClient()

    # Get pathway with desired file format
    kwargs = {
        'identifier': 'WP2062',
        'version': 0,
        'file_format': 'application/vnd.gpml+xml'
    }
    file = wikipathways_api_client_instance.get_pathway_as(**kwargs)
    print file

if __name__ == '__main__':
    main()
