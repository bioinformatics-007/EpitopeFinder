FROM python:3.9-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive

# --------------------------------------------------
# Working directory
# --------------------------------------------------
WORKDIR /Vaxelan_2_0

# --------------------------------------------------
# System dependencies
# --------------------------------------------------
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------
# Python dependencies
# --------------------------------------------------
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# --------------------------------------------------
# Install NCBI BLAST+ (stable HTTPS version)
# --------------------------------------------------
RUN mkdir -p /Vaxelan_2_0/tools/ncbi-blast \
    && wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-x64-linux.tar.gz \
       -O /tmp/blast.tar.gz \
    && tar -xzf /tmp/blast.tar.gz \
       -C /Vaxelan_2_0/tools/ncbi-blast --strip-components=1 \
    && rm /tmp/blast.tar.gz

# Add BLAST to PATH
ENV PATH="/Vaxelan_2_0/tools/ncbi-blast/bin:$PATH"

# --------------------------------------------------
# Copy full project (includes tools/clbtope)
# --------------------------------------------------
COPY . .

# --------------------------------------------------
# Tell VaxElan where ClbTope DB lives
# --------------------------------------------------
ENV CLBTOPE_DB=/Vaxelan_2_0/tools/clbtope/clbtope/Database

# --------------------------------------------------
# Run FastAPI Backend
# --------------------------------------------------
EXPOSE 8000
CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0", "--port", "8000"]

