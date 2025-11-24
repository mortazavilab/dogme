# Start from a base image with Conda pre-installed
FROM mambaorg/micromamba:1.5.6

# Switch to root user to install tools
USER root
RUN apt-get update && apt-get install -y wget procps && rm -rf /var/lib/apt/lists/*

# Copy our environment file
COPY environment.yml /tmp/environment.yml

# Install all the Conda tools
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# --- Install modkit binary ---
# We install v0.5.1. Note: This is the standalone version, so we don't need 
# external LibTorch environments like the manual installation might.
RUN wget https://github.com/nanoporetech/modkit/releases/download/v0.5.1-rc1/modkit_v0.5.1rc1_u16_x86_64.tar.gz && \
    tar -xzf modkit_v0.5.1rc1_u16_x86_64.tar.gz && \
    # 1. Create a new "home" for modkit
    mkdir -p /usr/local/lib/modkit && \
    # 2. Move ALL the extracted files
    # Note: Check the folder name inside the tarball if version changes
    mv dist_modkit_v0.5.1_8fa79e3/modkit /usr/local/lib/modkit/ && \
    mv dist_modkit_v0.5.1_8fa79e3/models /usr/local/lib/modkit/ && \
    mv dist_modkit_v0.5.1_8fa79e3/docs /usr/local/lib/modkit/ && \
    mv dist_modkit_v0.5.1_8fa79e3/README.md /usr/local/lib/modkit/ && \
    mv dist_modkit_v0.5.1_8fa79e3/LICENCE.txt /usr/local/lib/modkit/ && \
    # 3. Create a "shortcut" (symlink)
    ln -s /usr/local/lib/modkit/modkit /usr/local/bin/modkit && \
    # 4. Clean up
    rm modkit_v0.5.1rc1_u16_x86_64.tar.gz && \
    rm -r dist_modkit_v0.5.1_8fa79e3

# --- Install Dorado and its libraries ---
RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.2.0-linux-x64.tar.gz && \
    tar -xzf dorado-1.2.0-linux-x64.tar.gz && \
    mv dorado-1.2.0-linux-x64/bin/dorado /usr/local/bin/dorado && \
    mkdir -p /usr/local/dorado-lib/ && \
    mv dorado-1.2.0-linux-x64/lib/* /usr/local/lib/ && \
    rm -r dorado-1.2.0-linux-x64*

# --- ENVIRONMENT CONFIGURATION ---

# 1. Set MODKIT Variables (Adapted for Docker Paths)
# We point these to where we installed files in the container (/usr/local/lib/modkit)
ENV MODKITBASE="/usr/local/lib/modkit"
# Note: This points to the default model location inside the container. 
# If you need a specific custom model not included in the release, you must COPY it in.
ENV MODKITMODEL="${MODKITBASE}/models/r1041_e82_400bps_hac_v5.2.0@v0.1.0"

# 2. System Paths & Libraries
# LD_LIBRARY_PATH: Includes Dorado libs (/usr/local/lib) and Dorado extras
# We do NOT need DYLD (Mac only) or LIBTORCH (Standalone binary doesn't need it)
ENV LD_LIBRARY_PATH="/usr/local/lib:/usr/local/dorado-lib:${LD_LIBRARY_PATH}"

# 3. Add Conda to PATH
ENV PATH="/opt/conda/bin:${PATH}"

# (Optional) Ensure python alias exists
RUN ln -sf /opt/conda/bin/python3 /opt/conda/bin/python

# Set the entrypoint
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
