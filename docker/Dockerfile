# syntax=docker/dockerfile:1
# Build context is docker/, not the research data directory.
# CPU environment for Linux x86_64; Compose sets the matching platform.
FROM rocker/r-ver:4.6.1@sha256:4e6e6696ff54a86625c7f80046478e67014a2184e942b85eeff43db7b992d93e

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_JOBS=4
# GROMACS version used for the study, confirmed from the original runs.
ARG GROMACS_VERSION=2022.2
ARG VINA_VERSION=1.2.7

ENV VIRTUAL_ENV=/opt/egfr \
    PATH=/opt/egfr/bin:/opt/gromacs/bin:${PATH} \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg \
    OMP_NUM_THREADS=4

# Python 3.12 comes from the base image's Ubuntu 24.04 repositories.
# System headers support the Python extensions, R graphics, and Bioconductor.
RUN apt-get update && apt-get install -y --no-install-recommends \
      python3.12 python3.12-dev python3.12-venv \
      build-essential gfortran cmake pkg-config git curl wget ca-certificates \
      libfftw3-dev libblas-dev liblapack-dev libgfortran5 libgomp1 \
      libcurl4-openssl-dev libssl-dev libxml2-dev libbz2-dev liblzma-dev zlib1g-dev \
      libmagick++-dev librsvg2-dev libharfbuzz-dev libfribidi-dev \
      libfontconfig1-dev libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev \
      libicu-dev libglpk-dev libgsl-dev libhdf5-dev libnetcdf-dev libzmq3-dev \
      libgit2-dev libxt-dev libx11-dev libcairo2-dev libv8-dev \
      libxrender1 libxext6 libsm6 libgl1 \
      pandoc graphviz ghostscript fonts-dejavu-core tini \
    && rm -rf /var/lib/apt/lists/* \
    && python3.12 -m venv /opt/egfr

COPY requirements.tex /opt/egfr/requirements.tex
RUN python -m pip install pip==26.2.1 setuptools==84.0.0 wheel==0.48.0 \
    && python -m pip install --index-url https://download.pytorch.org/whl/cpu torch==2.14.0 \
    && python -m pip install --only-binary=acpype -r /opt/egfr/requirements.tex \
    && python -m pip check

# Build the manuscript's GROMACS version. SSE2 permits older x86_64 CPUs.
RUN curl --fail --location --retry 3 \
      "https://ftp.gromacs.org/gromacs/gromacs-${GROMACS_VERSION}.tar.gz" \
      -o /tmp/gromacs.tar.gz \
    && mkdir /tmp/gromacs-source \
    && tar -xzf /tmp/gromacs.tar.gz -C /tmp/gromacs-source --strip-components=1 \
    && cmake -S /tmp/gromacs-source -B /tmp/gromacs-build \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/gromacs \
      -DGMX_BUILD_OWN_FFTW=OFF -DGMX_FFT_LIBRARY=fftw3 \
      -DGMX_GPU=OFF -DGMX_MPI=OFF -DGMX_THREAD_MPI=ON -DGMX_SIMD=SSE2 \
      -DREGRESSIONTEST_DOWNLOAD=OFF \
    && cmake --build /tmp/gromacs-build --parallel "${BUILD_JOBS}" \
    && cmake --install /tmp/gromacs-build \
    && rm -rf /tmp/gromacs-source /tmp/gromacs-build /tmp/gromacs.tar.gz

# AutoDock Vina CLI, including the utility for splitting docking poses.
RUN curl --fail --location --retry 3 \
      "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VINA_VERSION}/vina_${VINA_VERSION}_linux_x86_64" \
      -o /usr/local/bin/vina \
    && curl --fail --location --retry 3 \
      "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VINA_VERSION}/vina_split_${VINA_VERSION}_linux_x86_64" \
      -o /usr/local/bin/vina_split \
    && chmod 755 /usr/local/bin/vina /usr/local/bin/vina_split

# Exact direct R package versions; grid ships with the pinned R runtime.
# Package sources: CRAN and Bioconductor 3.23 (the release for R 4.6).
# haploR is archived on CRAN, so use its latest archived release, 4.0.7.
RUN Rscript - <<'RSCRIPT'
options(repos = c(CRAN = "https://cran.r-project.org"), Ncpus = 4)
install.packages(
  "https://github.com/r-lib/remotes/archive/refs/tags/v2.5.0.tar.gz",
  repos = NULL, type = "source"
)
remotes::install_version("BiocManager", version = "1.30.27", upgrade = "never")
options(repos = BiocManager::repositories(version = "3.23"))
pins <- c(
  abind = "1.4-8", cowplot = "1.2.0", data.table = "1.18.6.1",
  DiagrammeR = "1.0.12", DiagrammeRsvg = "0.1", dplyr = "1.2.1",
  ggnewscale = "0.5.2", ggplot2 = "4.0.3", ggpubr = "1.0.0",
  gtable = "0.3.6", magick = "2.9.1", patchwork = "1.3.2",
  pheatmap = "1.0.13", png = "0.1-9", rsvg = "2.7.0",
  stringr = "1.6.0", survival = "3.8-12", tidyr = "1.3.2",
  tidyverse = "2.0.0", UpSetR = "1.4.1", IRkernel = "1.3.2",
  haploR = "4.0.7"
)
for (pkg in names(pins)) {
  remotes::install_version(pkg, version = pins[[pkg]],
                           dependencies = NA, upgrade = "never")
}
bioc_pins <- c(SummarizedExperiment = "1.42.0", TCGAbiolinks = "2.40.0")
for (pkg in names(bioc_pins)) {
  remotes::install_url(
    sprintf("https://bioconductor.org/packages/3.23/bioc/src/contrib/%s_%s.tar.gz",
            pkg, bioc_pins[[pkg]]),
    dependencies = NA, upgrade = "never"
  )
}
IRkernel::installspec(user = FALSE, sys_prefix = TRUE)
RSCRIPT

# Persist the resolved dependency lists inside the image for later reference.
RUN python -m pip freeze > /opt/egfr/python-installed.txt \
    && Rscript -e 'write.csv(installed.packages()[,c("Package","Version")], "/opt/egfr/r-installed.csv", row.names=FALSE)'

RUN if getent passwd 1000 >/dev/null; then \
      userdel "$(getent passwd 1000 | cut -d: -f1)"; \
    fi \
    && if getent group 1000 >/dev/null; then \
      groupmod --new-name egfr "$(getent group 1000 | cut -d: -f1)"; \
    else groupadd --gid 1000 egfr; fi \
    && useradd --uid 1000 --gid egfr --create-home --shell /bin/bash egfr \
    && mkdir -p /project /home/egfr/.cache \
    && chown -R egfr:egfr /project /home/egfr
USER egfr
WORKDIR /project
EXPOSE 8888
ENTRYPOINT ["/usr/bin/tini", "--"]
# Keep Jupyter's generated-token authentication enabled; never bake a token here.
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--ServerApp.root_dir=/project"]
