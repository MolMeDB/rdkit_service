FROM continuumio/miniconda3

WORKDIR /usr/src/app

RUN pip3 install --no-cache-dir requests
RUN pip3 install rdkit-pypi
RUN pip3 install --no-cache-dir dimorphite-dl

COPY . .

EXPOSE 9696

CMD [ "python3", "-u", "./main.py" ]
