from wikipathways_api_client import WikipathwaysApiClient

def main():
    wikipathways_api_client_instance = WikipathwaysApiClient()

    # Get colored pathway
    kwargs = {
        'identifier': 'WP2062',
        'version': 0,
        'element_identifiers': ["ffffff90","ffffffe5"],
        'colors': ["#0000FF","#0000FF"],
        'file_format': 'image/svg+xml'
    }
    file = wikipathways_api_client_instance.get_colored_pathway(**kwargs)
    print file

if __name__ == '__main__':
    main()

