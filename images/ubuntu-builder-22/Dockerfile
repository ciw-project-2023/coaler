FROM ubuntu:22.04

RUN apt update
RUN apt upgrade -y
RUN apt install cmake python3 python3-pip -y
RUN pip3 install conan
RUN conan profile detect
RUN conan remote add doc http://server.conan.corealigner.de
RUN conan remote login doc ciw -p VZKhzh2v5nCnijAS3A8R

WORKDIR /code