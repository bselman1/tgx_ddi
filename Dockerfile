FROM rocker/rstudio:4.4.1

RUN  apt-get update && \
     apt-get install --yes --no-install-recommends \
        zlib1g-dev \
     && rm -rf /var/lib/apt/lists/*

# Install R dependencies
RUN install2.r devtools BiocManager Cairo heatmap3 dendextend openxlsx BH plogr RSQLite

# Install BioConductor packages
RUN --mount=type=bind,source=bioc_pkgs.R,target=bioc_pkgs.R Rscript --no-save bioc_pkgs.R
