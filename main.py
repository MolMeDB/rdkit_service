import requests
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import base64
import logging
import datetime
import xml.etree.ElementTree as ET
from urllib.parse import urlparse, unquote
from services import RDKIT as R


class SERVICE(ThreadingHTTPServer):
    # Server initialization 
    def __init__(self, port):
        ThreadingHTTPServer.__init__(self, ("", port), SERVICE_HANDLER)


# Handler for proccessing HTTP requests
class SERVICE_HANDLER(BaseHTTPRequestHandler):

    # Checks uri validity
    def is_valid_path(self, path):
        valid_paths = [
            "smiles/canonize",
            "3dstructure/generate",
            "2dstructure/generate",
            "makeInchi",
            "general"
        ]

        return path in valid_paths

    def init(self):
        self.RDKIT = R.RDKIT()


    # Proccess GET request
    def do_GET(self):
        self.init()

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
            self.answer("", status)
            return


        # Check input
        if not self.is_valid_path(uri):
            print("Invalid request. Requested path: " + uri)
            self.answer("Invalid request. Requested path: " + uri, 400)
            return

        print("Initialized... Requested path: " + uri)

        try:
            if uri == "smiles/canonize":
                output = self.RDKIT.canonizeSmiles(request_params)

            elif uri == "3dstructure/generate":
                output = self.RDKIT.make3Dstructure(request_params)

            elif uri == "2dstructure/generate":
                output = self.RDKIT.make2Dstructure(request_params)

            # Makes inchi from smiles
            elif uri == "makeInchi":
                output = self.RDKIT.makeInchi(request_params)

            # Get general info about molecule
            elif uri == "general":
                output = self.RDKIT.getGeneralInfo(request_params)

        except Exception as e:
            self.answer(e, 404)
            return

        # Send answer
        self.answer(output, status)


    # Server answer
    def answer(self, response, code):
        if isinstance(response, Exception):
            response = response.args[0]

        if not isinstance(response, list) and not isinstance(response, dict):
            response = [response]

        output_data = json.dumps(response).encode()

        self.send_response(code)

        self.send_header("Content-Length", output_data.__len__())
        self.send_header('Content-type','application/json')
        self.end_headers()

        # Send data
        self.wfile.write(output_data)





# Where to listen
port = 9696

RDKIT_SERVER = SERVICE(port)

RDKIT_SERVER.serve_forever()


