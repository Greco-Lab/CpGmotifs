
FROM openanalytics/r-base

MAINTAINER Giovanni Scala "grecolab.fi@gmail.com"

# system libraries of general use
RUN apt-get update && apt-get install -y \
     sudo \
     pandoc \
     pandoc-citeproc \
     libexpat1-dev \
     libcairo2-dev \
     libxt-dev \
     libssl-dev \
     libssh2-1-dev \
     libssl1.0.0 \
     libcurl4-openssl-dev \
     libxml2-dev \
     ghostscript


# basic shiny functionality
RUN R -e "install.packages(c('curl','RCurl','shiny', 'shinyjs', 'shinydashboard', 'readr', 'DT', 'tibble', 'gplots',\
                              'dendextend', 'foreach', 'doParallel', 'XML', 'BiocManager'), repos='https://cloud.r-project.org/')"

# install bioconductor packages
# RUN R -e "BiocManager::install(c('minfi','BSgenome.Hsapiens.UCSC.hg19','IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b2.hg19', 'Biostrings'), version = '3.8')"

RUN R -e "BiocManager::install(c('minfi','BSgenome.Hsapiens.UCSC.hg19','IlluminaHumanMethylation450kanno.ilmn12.hg19', 'IlluminaHumanMethylationEPICanno.ilm10b2.hg19', 'Biostrings'))"

# install perl libraries (needed by meme)
RUN cpan File::Which HTML::PullParser HTML::Template HTML::TreeBuilder JSON XML::Simple XML::Parser::Expat

# install meme suite
RUN wget http://meme-suite.org/meme-software/5.0.5/meme-5.0.5.tar.gz && tar -xvf meme-5.0.5.tar.gz && rm meme-5.0.5.tar.gz && cd meme-5.0.5/ && ./configure --prefix=/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt && make && make install && cd .. && rm -r meme-5.0.5

# dowload motifs meme database
RUN wget http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz && tar -xvf motif_databases.12.18.tgz && rm motif_databases.12.18.tgz

# copy the app to the image
RUN mkdir /root/app
COPY app /root/app

COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/app')"]