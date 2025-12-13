# Set base image
FROM --platform=linux/amd64 rocker/shiny

# My authorship
LABEL maintainer="ehill@iolani.org"
LABEL version="1.0.0"
LABEL description="EPI2MEViz - Reboot for Iolani School"

# Convenience packages
RUN apt update && \
    apt install -y curl git wget nano libglpk-dev
# build-essential libcurl4-openssl-dev libssl-dev libxml2-dev libcairo2-dev libxt-dev libblas-dev liblzma-dev libbz2-dev 

# Install R packages from CRAN
RUN R -e "install.packages(c('shiny','ggplot2','dplyr','stringr','plotly','shinydashboard','shinyBS','htmlwidgets','tidyr','shinyjs','vegan','shinyalert','BiocManager'))"

# Install R packages from Bioconductor
RUN R -e "BiocManager::install(c('phyloseq','microbiome'))"

# Copy app to image
RUN rm -r /srv/shiny-server/*
COPY epi2meviz.R /srv/shiny-server/app.R

# Fix permissions
RUN chown -R shiny:shiny /srv/shiny-server