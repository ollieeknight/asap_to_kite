FROM python:3.11-slim

LABEL description="ASAP-seq ADT tools: asap_to_kite, kite featuremap"

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
        git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install asap_to_kite from GitHub
RUN pip install --no-cache-dir \
        "git+https://github.com/ollieeknight/asap_to_kite.git"

# Install kite featuremap
RUN git clone --depth 1 https://github.com/pachterlab/kite /opt/kite \
    && sed -i '1s|.*|#!/usr/bin/env python3|' /opt/kite/featuremap/featuremap.py \
    && chmod +x /opt/kite/featuremap/featuremap.py \
    && ln -s /opt/kite/featuremap/featuremap.py /usr/local/bin/featuremap
