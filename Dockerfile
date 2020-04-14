FROM mcs07/rdkit:latest

WORKDIR /usr/src/app

RUN apt-get update  && apt-get install -yq --no-install-recommends python3-pip

RUN pip3 install --no-cache-dir requests

COPY . .

EXPOSE 9696

CMD [ "python3", "-u", "./main.py" ]
