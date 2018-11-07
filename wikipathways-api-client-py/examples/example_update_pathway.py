###
# Test script for new WikiPathways webservice API
# author: msk (mkutmon@gmail.com)
###

import os
import requests
import getpass
from lxml import etree as ET

##################################
# variables

current_directory = os.path.dirname(os.path.realpath(__file__))
gpml_file = current_directory + '/test.gpml'
wp_id = 'WP4'
base_iri = 'http://pvjs.wikipathways.org/wpi/webservicetest/'
description = 'test update webservice function'

# current revision of the pathway! update always to the latest revision on the pathway page!
revision = "78853"

##################################

# define namespaces
namespaces = {'ns1':'http://www.wso2.org/php/xsd','ns2':'http://www.wikipathways.org/webservice'}

# login
username = raw_input('Username: ')
password = getpass.getpass('Password:')
auth = {'name' : username , 'pass' : password}
login_request = requests.get(base_iri + 'login', params=auth)
dom = ET.fromstring(login_request.text)

authentication = ''
for node in dom.findall('ns1:auth', namespaces):
	authentication = node.text

# read gpml file
f = open(gpml_file, 'r')
gpml = f.read()

# update pathway
update_params = {'pwId': wp_id, 'revision' : revision, 'auth' : authentication, 'username': username}
update_data = {'description' : description, 'gpml' : gpml}
update_pathway_request = requests.post(base_iri + 'updatePathway', params=update_params, data=update_data)
print 'Completed with response:'
print update_pathway_request.status_code
print update_pathway_request.text
