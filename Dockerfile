FROM fnndsc/civet:2.1.1

COPY . /usr/local/src/

# CIVET base image is Ubuntu 18.04.4 LTS (old)
# On non-x86_64 architectures, numpy (a dependency of pybicpl)
# must be compiled from a source distribution.
# python-pip3 package installs pypip 9.0.1, which has a bug
# where it fails to resolve Cython as a build dependency of numpy.
# The cleanest solution here would be to bump pip to version 20.

RUN apt-get update -qq && apt-get install -qq python3-pip \
    && pip3 install --upgrade pip \
    && pip3 install /usr/local/src

CMD ["surfaces_fetus", "--help"]

