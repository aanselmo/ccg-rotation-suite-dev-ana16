# CCG Rotation Suite
These are various reporting scripts that you have to run when you're running production projects.

## Directions for Sarek and it's scripts
Please see the top level [README for Sarek](sarek_summary_reports/docs/README.md)


## Install Directions for Debian 13 on Nextflow4
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
sudo -H /usr/local/miniforge3/bin/python
Python 3.12.12 | packaged by conda-forge | (main, Oct 22 2025, 23:25:55) [GCC 14.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from SigProfilerMatrixGenerator import install as genInstall
>>> genInstall.install('GRCh37', rsync=False, bash=True)
```
4. Install [docker](https://docs.docker.com/engine/install/debian/):
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
    ```bash
    sudo usermod -a -G docker [username]
    ```
    * Log out on all terminals and log back in.