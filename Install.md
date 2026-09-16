# Install Directions

## Prod ccg-rotation-suite Install Directions
1. Download Miniforge3:
```bash
wget https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-x86_64.sh
sudo bash Miniforge3-24.11.3-2-Linux-x86_64.sh
# Please install to /usr/local/miniforge3 directory
```
3. Install Python requirements/dependencies:
```bash
sudo -H /usr/local/miniconda3/bin/pip install -r requirements.txt
```
4. Install GRCh37 for SigProfilerMatrixGenerator:
```python
sudo -H /usr/local/miniforge3/bin/python << EOF
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh37', rsync=True, bash=True)
EOF
```
5. Install poppler. This is required by pdf2image in requirements.txt for python.
```bash
sudo apt install poppler-utils
```
6. Set the TimeZone
```bash
sudo timedatectl set-timezone America/New_York
```

## Dev ccg-rotation-suite Install Directions
3. Install Python requirements/dependencies:
```bash
sudo -H /usr/local/miniconda3/bin/pip install -r requirements-dev.txt
```

## Sarek Install Directions for Debian 13 on Nextflow4
In addition to installing ccg-rotation-suite on the nextflow4 head node, I also installed:
1. Install [docker](https://docs.docker.com/engine/install/debian/):
    * Install docker apt repo:
    ```bash
    # Add Docker's official GPG key:
    sudo apt update
    sudo apt install ca-certificates curl
    sudo install -m 0755 -d /etc/apt/keyrings
    sudo curl -fsSL https://download.docker.com/linux/debian/gpg -o /etc/apt/keyrings/docker.asc
    sudo chmod a+r /etc/apt/keyrings/docker.asc

    # Add the repository to Apt sources:
    sudo tee /etc/apt/sources.list.d/docker.sources <<EOF
    Types: deb
    URIs: https://download.docker.com/linux/debian
    Suites: $(. /etc/os-release && echo "$VERSION_CODENAME")
    Components: stable
    Signed-By: /etc/apt/keyrings/docker.asc
    EOF

    sudo apt update
    ```
    * Install docker
    ```bash
    sudo apt install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
    ```
    * Add yourself to the docker group
    ```bash/opt/gatk/4.6.1.0/gatk-package-4.6.1.0-local.jar
    sudo usermod -a -G docker [username]
    ```
    * Log out on all terminals and log back in.

2. Install java version 17 or later.
3. Install GATK to: `/opt/gatk/4.6.1.0/gatk-package-4.6.1.0-local.jar`
4. Install nextflow.
5. Install Sarek.