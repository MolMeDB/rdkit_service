FROM continuumio/miniconda3:23.9.0-0

WORKDIR /usr/src/app

RUN pip3 install --no-cache-dir requests
RUN pip3 install rdkit
RUN pip3 install --no-cache-dir dimorphite-dl
# RUN pip3 install gssapi
# RUN pip3 install awscli --ignore-installed six
RUN pip3 install setuptools paramiko

COPY . .

EXPOSE 9696

CMD [ "python3", "-u", "./main.py" ]
