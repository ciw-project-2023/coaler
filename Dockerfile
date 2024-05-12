FROM ubuntu:22.04

RUN apt update
RUN apt upgrade -y
RUN apt install git build-essential python3 python3.10-venv -y

RUN git clone https://github.com/ciw-project-2023/coaler.git /code

WORKDIR /code
RUN ./configure
RUN make
RUN make install

RUN mkdir /work
WORKDIR /work

ENTRYPOINT ["/usr/local/bin/coaler"]