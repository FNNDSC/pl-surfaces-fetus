FROM fnndsc/civet:2.1.1

COPY . /usr/local/src/

RUN apt-get update -qq && apt-get install -qq python3-pip \
    && pip3 install /usr/local/src

CMD ["surfaces_fetus", "--help"]
