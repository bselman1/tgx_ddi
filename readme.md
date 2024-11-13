## TGx-DDI Biomarker for DNA Damage Classification
The TGx-DDI biomarker was developed as a toxicogenomics signature to identify chemicals that can cause DNA damage in human cells in culture. Development of this biomarker was as a collaborative effort between the Genomics Committee of the Health and Environmental Sciences Institute (HESI), Health Canada, the University of Ottawa, and Georgetown University. The NIEHS Division of Translational Toxicology (DTT) collaborated on the development of the classification tool.

## Overview
This repo is a containerized version of the TGx-DDI Biomarker for DNA Damage Classification tool (https://cebs.niehs.nih.gov/tgxddi/) with some modifications to better allow running multiple chemicals through the classifier.

### Interactive Use
This use case fits best if you want to make changes to the data and or code in real time and have those changes be reflected in the docker container. This setup will mount a working directory on the host computer into the container so any changes on either side (host or in RStudio) will change the file on the host computer.

#### Building the image
```bash
docker build --tag tgx_ddi .
```

#### Running an interactive RStudio container
```bash
docker run -it --rm -e PASSWORD=test123 -p 8788:8787 --mount type=bind,source="$(pwd)/src",destination=/home/rstudio/project tgx_ddi
```