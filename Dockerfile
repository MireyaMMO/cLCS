FROM python:3.11.7-slim
LABEL Maintainer="mireya.mmor@gmail.com"

RUN apt update && apt install -y libgeos-dev libproj-dev libnetcdf-dev git 
RUN pip install shapely==1.8.5.post1 --no-binary shapely


RUN echo "-----------------Install OpenDrift -----------------" &&\
    git clone https://github.com/simonweppe/opendrift.git /source/opendrift 

RUN echo "-----------------Install cLCS -----------------" &&\
    git clone https://github.com/MireyaMMO/cLCS.git /source/cLCS &&\
    cd /source/cLCS &&\
    pip install -r requirements.txt &&\
    pip install --no-deps -e . &&\
    cd /source/opendrift &&\
    pip install --no-deps -e .



