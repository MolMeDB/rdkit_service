from rdkit import Chem
import requests
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import base64
import logging
import datetime
import xml.etree.ElementTree as ET
from urllib.parse import urlparse, unquote


class RDKIT(ThreadingHTTPServer):
    # Server initialization 
    def __init__(self, port):
        ThreadingHTTPServer.__init__(self, ("", port), RDKIT_HANDLER)


# Handler for proccessing HTTP requests
class RDKIT_HANDLER(BaseHTTPRequestHandler):

    # ----- Not used ----
    # def auth(self):
    #     # Authorization info
    #     authorization = self.headers.get('Authorization', None)

    #     if authorization:
    #         auth_type, auth_data = authorization.split()

    #         if auth_type == 'Basic':
    #             auth_user, auth_pswd = str(
    #                 base64.b64decode(auth_data),
    #                 'UTF-8'
    #             ).split(':')

    #             return (
    #                 auth_user == self.server.auth_user
    #                 and auth_pswd == self.server.auth_pswd
    #             )

    #         else:
    #             return False

    #     else:
    #         return False

    # Checks uri validity
    def is_valid_path(self, path):
        valid_paths = [
            "smiles/canonize"
        ]

        return path in valid_paths

    # Canonize SMILES
    def canonizeSmiles(self, params = {}):
        if not "smi" in params:
            print("Parameter 'smi' not found in argument list")
            return []
        
        req_smiles = params["smi"]

        canonized_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(req_smiles))

        return [canonized_smiles]

    # Proccess GET request
    def do_GET(self):

        # Each output should be array
        output = [] 

        # Proccess GET params
        query = urlparse(self.path).query

        if query and "=" in query:
            request_params = dict(qc.split("=") for qc in query.split("&"))
        else:
            request_params = {}

        # UNQUOTE params
        for key in request_params:
            request_params[key] = unquote(request_params[key])

        # Default response status code
        status = 200

        # Parsing URL
        path = self.path
        path = path.strip("?") + "?"
        uri = self.path[:path.find("?")].strip("/")

        # Connection test
        if uri == "test":
            self.send_response(status)
            self.send_header('Content-type','application/json')
            self.end_headers()
            return


        # Check input
        if not self.is_valid_path(uri):
            print("Not valid request. Requested path: " + uri)
            self.send_response(400)
            self.send_header('Content-type','application/json')
            self.end_headers()
            return

        print("Initialized... Requested path: " + uri)

        # Edit in future updates
        if uri == "smiles/canonize":
            output = self.canonizeSmiles(request_params)

        output_data = json.dumps(list(output)).encode()

        # Send header info
        self.send_response(status)

        self.send_header("Content-Length", output_data.__len__())
        self.send_header('Content-type','application/json')
        self.end_headers()

        # Send data
        self.wfile.write(output_data)
        





# Where to listen
port = 9696

RDKIT_SERVER = RDKIT(port)

RDKIT_SERVER.serve_forever()


