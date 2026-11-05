# Image for export_nexus_ogt.py: CNVkit (for `cnvkit.py export nexus-ogt`)
# plus numpy / pandas / matplotlib (for plot_ogt_v3.py, which is baked in).
#
# Build from a folder containing this Dockerfile, plot_ogt_v3.py and
# cytoBand.txt.gz (UCSC, same genome build as your samples):
#
#     docker build -t cnv-ogt:0.9.13 .
#
# Check it:
#     docker run --rm cnv-ogt:0.9.13 cnvkit.py version
#     docker run --rm cnv-ogt:0.9.13 python /opt/cnv/plot_ogt_v3.py -h

FROM python:3.12-slim

# CNVkit pulls in numpy, pandas, matplotlib, scipy and pysam (for VCF reading).
# Pin it to the version your sarek run used if that differs.
ARG CNVKIT_VERSION=0.9.13
RUN pip install --no-cache-dir "cnvkit==${CNVKIT_VERSION}"

# matplotlib needs a writable config dir; the container runs as your UID,
# which has no home directory inside the image.
ENV MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    HOME=/tmp \
    PYTHONDONTWRITEBYTECODE=1

COPY plot_ogt_v3.py /opt/cnv/plot_ogt_v3.py
COPY cytoBand.txt.gz /opt/cnv/cytoBand.txt.gz

WORKDIR /tmp
CMD ["cnvkit.py", "version"]
