from os import access
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json, base64, shutil
import xml.etree.ElementTree as ET
from urllib.parse import urlparse, unquote
from services import RDKIT as R, ssh


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
            "smiles/allCharges",
            "conformers/toCosmo",
            "3dstructure/generate",
            "2dstructure/generate",
            "makeInchi",
            "general",
            "mmpa/fragment",
            "mol/similarity",
            "mol/fingerprint",
            "cosmo/runningJobs",
            "cosmo/getResults",
            "cosmo/clearFolderStructure",
            "cosmo/uploadSDF",
            "cosmo/checkOptimizationStatus",
            "cosmo/optimizeSDF",
            "cosmo/run"
        ]

        return path in valid_paths

    def init(self):
        self.RDKIT = R.RDKIT()


    # Proccess GET request
    def do_GET(self):
        self.init()
        
        auth = self.headers.get('Authorization')
        if auth and str(auth).strip().startswith("Basic"):
            auth = str(auth).replace("Basic", "").strip()
            auth = base64.b64decode(auth).decode('utf-8')
            
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

        ans = True
        asHTML = False
        asFile = False
        
        if uri.startswith("cosmo"):
            if "server" in request_params:
                server = request_params["server"]
            else:
                server = "zuphux.metacentrum.cz"
                
            # Get username and password from HTTP header
            auth = str(auth).split(":")
            if not len(auth) == 2:
                self.answer("Invalid credentials", 401)
                return
            
            try:
                if not ssh.connect(
                    server, 
                    auth[0], 
                    auth[1], 
                    ignoreSFTP=request_params["ignoreSFTP"] if "ignoreSFTP" in request_params else False
                ):
                    self.answer("Invalid credentials", 401)
                    return
            except Exception as e:
                print(e)
                self.answer("Invalid credentials", 401)

        try:
            if uri == "smiles/canonize":
                output = self.RDKIT.canonizeSmiles(request_params)

            elif uri == "smiles/allCharges":
                output = self.RDKIT.getAllChargeSmiles(request_params)

            elif uri == "conformers/toCosmo":
                output = self.RDKIT.COSMO_conformers(request_params)

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

            # MMPA - Get molecule substructures
            elif uri == "mmpa/fragment":
                ans, output, asHTML = self.RDKIT.mmpaFragment(request_params, self)

            # MMPA - Get molecules similarity
            elif uri == "mol/similarity":
                ans, output = self.RDKIT.computeSimilarity(request_params)

            # MMPA - Get molecule fingerprint
            elif uri == "mol/fingerprint":
                output = self.RDKIT.getFingerprint(request_params)

            # COSMO - Download results
            elif uri == "cosmo/runningJobs":
                output = ssh.getRunningJobs(request_params)
                
            # COSMO - Download results
            elif uri == "cosmo/getResults":
                output = ssh.cosmo_download_results(request_params)
                asFile = isinstance(output, str)

            # COSMO - Clear folder structure => for FORCE re-computation
            elif uri == "cosmo/clearFolderStructure":
                output = ssh.cosmo_clear_folder_structure(request_params)
                
            # COSMO - Upload SDF files to remote server
            elif uri == "cosmo/uploadSDF":
                output = ssh.cosmo_upload_sdf(request_params)
                
            # COSMO - Check state of optimization
            elif uri == "cosmo/checkOptimizationStatus":
                output = ssh.cosmo_check_optimization_status(request_params)
                
            # COSMO - Check state of optimization
            elif uri == "cosmo/optimizeSDF":
                output = ssh.cosmo_optimize_sdf(request_params)
    
            # COSMO - Run COSMO
            elif uri == "cosmo/run":
                output = ssh.runCosmo(request_params)
    
        except Exception as e:
            print(e)
            self.answer(e, 404)
            return

        # Send answer
        if ans:
            self.answer(output, status, asHTML, asFile=asFile)


    # Server answer
    def answer(self, response, code, asHTML = False, asFile = False):
        if isinstance(response, Exception):
            response = response.args[0]
            
        if asFile:
            import os
            if not os.path.exists(response):
                self.send_response(404, "File not found")
                return
            self.send_response(code)

            self.send_header("Content-Length", os.path.getsize(response))
            self.send_header('Content-Disposition','attachment; filename="archive.zip')
            self.send_header('Content-type','application/zip')
            output_data = open(response, 'rb').read()
            
            # Finally, remove the output data
            # Remove whole created folder structure
            root = response.split('/')[0]
            if os.path.isdir(root):
                shutil.rmtree(root)

        elif not asHTML:
            if not isinstance(response, list) and not isinstance(response, dict):
                response = [response]

            output_data = json.dumps(response).encode()

            self.send_response(code)

            self.send_header("Content-Length", output_data.__len__())
            self.send_header('Content-type','application/json')
        else: 
            response = str(response).encode()

            output_data = response
            self.send_response(code)

            self.send_header("Content-Length", output_data.__len__())
            self.send_header('Content-type','text/html')

        self.end_headers()

        # Send data
        self.wfile.write(output_data)


# Where to listen
port = 9696

RDKIT_SERVER = SERVICE(port)

RDKIT_SERVER.serve_forever()


