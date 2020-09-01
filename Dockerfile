FROM debian:stretch

ENV DEBIAN_FRONTEND=noninteractive

COPY /scripts /scripts

RUN /scripts/install-debian-packages

ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV LC_ALL en_US.UTF-8

#-----------------------------------------------------------------------------------------------------
# R
#-----------------------------------------------------------------------------------------------------

ENV R_VERSION=3.6.0
ENV R_ARCHIVE_URL https://cran.r-project.org/src/base/R-3/R-${R_VERSION}.tar.gz

RUN set -ex \
  && curl -sSL ${R_ARCHIVE_URL} | tar -zxf - \
  && cd R-${R_VERSION} && \
  \
  ./configure --prefix=/usr/local \
              --localstatedir=/var \
              --disable-nls \
              --with-readline \
              --without-x \
              --without-recommended-packages \
              --enable-memory-profiling \
              --enable-R-shlib \
              --with-tcltk \
              --with-tcl-config=/usr/lib/tclConfig.sh \
              --with-tk-config=/usr/lib/tkConfig.sh \
  \
  && make -j $(nproc) \
  && make install \
  && cd src/nmath/standalone \
  && make -j $(nproc) \
  && make install \
  && cd / && rm -rf R-${R_VERSION}

RUN set -ex \
  && echo "R_LIBS_SITE=/usr/local/lib/R/site-library" >> /usr/local/lib/R/etc/Renviron \
  && echo "R_LIBS_USER=/usr/local/lib/R/site-library" >> /usr/local/lib/R/etc/Renviron \
  && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' >> /usr/local/lib/R/etc/Rprofile.site

#-----------------------------------------------------------------------------------------------------
# Build R packages
#-----------------------------------------------------------------------------------------------------

ENV R_PACKAGES \
  ggplot2 cowplot svglite gridExtra RColorBrewer reshape argparse

RUN echo 'withCallingHandlers(install.packages(commandArgs(trailingOnly = TRUE)), warning=function(w) stop(w))' \
      | R --slave --no-restore --args ${R_PACKAGES}

#-----------------------------------------------------------------------------------------------------
# Install bamtools (telseq depenency)
#-----------------------------------------------------------------------------------------------------

ENV BAMTOOLS_VERSION 2.5.1
ENV BAMTOOLS_URL https://github.com/pezmaster31/bamtools/archive/v${BAMTOOLS_VERSION}.tar.gz

RUN set -ex \
   && curl -sSL ${BAMTOOLS_URL} | tar zxf - \
   && cd bamtools-${BAMTOOLS_VERSION} \
   && mkdir build \
   && cd build \
   && cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. \
   && make install \
   && cd ../.. \
   && rm -rf bamtools-${BAMTOOLS_VERSION}

#-----------------------------------------------------------------------------------------------------
# Install telseq
#-----------------------------------------------------------------------------------------------------

ENV TELSEQ_VERSION 0.0.2
ENV TELSEQ_URL https://github.com/zd1/telseq/archive/v${TELSEQ_VERSION}.tar.gz

RUN set -ex \
  && curl -sSL ${TELSEQ_URL} | tar zxf - \
  && cd telseq-${TELSEQ_VERSION}/src \
  && ./autogen.sh \
  && ./configure --with-bamtools=/usr/local \
  && make install \
  && cd ../.. \
  && rm -rf telseq-${TELSEQ_VERSION}

#-----------------------------------------------------------------------------------------------------
# Install telomerecat
#-----------------------------------------------------------------------------------------------------

RUN set -ex \
  && apt-get -q update \
  && apt-get -y install python-pip

RUN set -ex \
  && pip2 install telomerecat pysam PyPDF2

ENV SAMTOOLS_VERSION 1.9
ENV HTSLIB_VERSION 1.9
ENV SAMBLASTER_VERSION 0.1.24

ENV HTSLIB_URL https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
ENV SAMBLASTER_URL https://github.com/GregoryFaust/samblaster/releases/download/v.${SAMBLASTER_VERSION}/samblaster-v.${SAMBLASTER_VERSION}.tar.gz
ENV SAMTOOLS_URL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2

#-------------------------------------------------------------------------
# Install samtools
#-------------------------------------------------------------------------

RUN set -ex \
  && curl -sSL ${HTSLIB_URL} | tar -jxf - \
  && curl -sSL ${SAMTOOLS_URL} | tar -jxf - \
  \
  && make -C htslib-${HTSLIB_VERSION} \
  \
  && HTSDIR=$(pwd)/htslib-${HTSLIB_VERSION} make -C samtools-${SAMTOOLS_VERSION} \
  && /bin/cp -p samtools-${SAMTOOLS_VERSION}/samtools /usr/local/bin/samtools \
  \
  && /bin/rm -rf samtools-${SAMTOOLS_VERSION}

ENV PATH="htslib-{HTSLIB_VERSION}:${PATH}"

RUN set -ex && apt-get install -y zip

#-----------------------------------------------------------------------------------------------------
# Install TelomereHunter
#-----------------------------------------------------------------------------------------------------

RUN set -ex \
  && pip2 install telomerehunter

#-----------------------------------------------------------------------------------------------------
# Install Computal
#-----------------------------------------------------------------------------------------------------

RUN set -ex \
  && git clone https://github.com/lilit-nersisyan/computel.git \
  && chmod +x /computel/bin/bowtie2* \
  && chmod +x /computel/computel.sh \
  && chmod +x /computel/setup_test/* \
  && chmod +x /computel/src/scripts/*

#-----------------------------------------------------------------------------------------------------
# Install qMotif
#-----------------------------------------------------------------------------------------------------

RUN set -ex \
  && wget https://sourceforge.net/projects/adamajava/files/qmotif-1.2.tar.bz2 \
  && mv qmotif* /usr/local/bin \
  && cd /usr/local/bin \
  && tar -xjf qmotif* \
  && cd /

#-----------------------------------------------------------------------------------------------------
# Install Motif_counter
#-----------------------------------------------------------------------------------------------------

RUN echo 'install.packages("BiocManager")' | R --slave --no-restore --args
RUN echo 'BiocManager::install(c("motifcounter"), update=TRUE, ask=FALSE)' | R --slave --no-restore --args


ENV R_PACKAGES \
  seqinr psych

RUN echo 'withCallingHandlers(install.packages(commandArgs(trailingOnly = TRUE)), warning=function(w) stop(w))' \
      | R --slave --no-restore --args ${R_PACKAGES}

RUN set -ex \
  && sed -i 's/ans=""/ans="y"/g' computel/computel.sh \
  && sed -i 's/read ans//g' /computel/computel.sh

RUN set -ex \
  && apt-get install -y bc

#-----------------------------------------------------------------------------------------------------
# Install picard
#-----------------------------------------------------------------------------------------------------

RUN set -ex \
  && wget https://github.com/broadinstitute/picard/releases/download/2.21.3/picard.jar \
  && chmod +x picard.jar \
  && mv picard.jar /usr/local/bin

#-----------------------------------------------------------------------------------------------------
# Make mountpoints
#-----------------------------------------------------------------------------------------------------

RUN mkdir -p \
  /hii/work \
  /labdata \
  /resdata \
  /shares

#-----------------------------------------------------------------------------------------------------
# Default command to execute
#-----------------------------------------------------------------------------------------------------

CMD [ "/bin/bash" ]

# vim: ft=Dockerfile
