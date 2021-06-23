# This Dockerfile constructs a docker image that contains an installation
# of the coinsparse library
#
# Example build:
#   docker build --no-cache --tag richford/coinsparse `pwd`
#
#   (but really, use docker-compose up instead).
#
FROM continuumio/miniconda3:4.9.2-alpine

MAINTAINER Adam Richie-Halford <nben@uw.edu>

RUN mkdir /home/coinsparse
COPY ./setup.py ./setup.cfg ./pyproject.toml ./LICENSE ./README.md \
     /home/coinsparse/
COPY ./coinsparse/coinsparse.py /home/coinsparse/coinsparse/coinsparse.py
# RUN echo "version = '0.1.docker'" > /home/coinsparse/coinsparse/_version.py
# RUN echo "version_tuple = (0, 1, 'docker')" >> /home/coinsparse/coinsparse/_version.py
# RUN echo __version = '0.1.docker' > /home/coinsparse/coinsparse/__init__.py
RUN cd /home/coinsparse && python -m pip install -e . 

RUN mkdir /coinsparse_data
WORKDIR /coinsparse_data

ENTRYPOINT ["coinsparse"]