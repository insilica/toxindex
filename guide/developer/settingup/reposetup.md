---
title: Setting up Toxindex Repo in VM
sidebar_position: 3
---

# Setting up Toxindex Repo in VM

This guide walks you through setting up the ToxIndex repository in a virtual machine environment.

## 1. Install Git

```bash
sudo apt update && sudo apt install -y git
```

## 2. Add SSH Key to GitHub

### Generate SSH Key
```bash
ssh-keygen -t ed25519 -C "kyu@toxindex-prod" -f ~/.ssh/id_ed25519 -N ""
```

### Display Public Key
```bash
cat ~/.ssh/id_ed25519.pub
```

### Add Key to GitHub
1. Go to [github.com](https://github.com/) > Settings > SSH and GPG keys
2. Click "New SSH key"
3. Paste the public key content
4. Authenticate via GitHub mobile app if prompted

### Clone Repository
```bash
git clone git@github.com:insilica/toxindex.git
```

## 3. Install Git LFS

```bash
# Add Git LFS repository
curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash

# Install Git LFS
sudo apt-get install git-lfs

# Initialize Git LFS
git lfs install

# Track LFS files
git lfs track

# Check LFS status
git lfs ls-files
git lfs status

# Pull LFS files
git lfs pull
```

## 4. Install Nix

```bash
# Install Nix
sh <(curl -L https://nixos.org/nix/install) --no-daemon

# Create Nix config directory
mkdir -p ~/.config/nix

# Enable experimental features
echo "experimental-features = nix-command flakes" >> ~/.config/nix/nix.conf

# Source profile
source ~/.profile

# Verify installation
nix --version
```

## 5. Install Docker

```bash
# Update package index
sudo apt update

# Install prerequisites
sudo apt install -y ca-certificates curl gnupg lsb-release

# Add Docker's official GPG key
sudo mkdir -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg

# Add Docker repository
echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# Update package index
sudo apt update

# Install Docker
sudo apt install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# Enable and start Docker
sudo systemctl enable --now docker

# Add user to docker group
sudo usermod -aG docker $USER
```

## Next Steps

After completing these steps:

1. **Log out and log back in** to apply the Docker group changes
2. **Verify Docker installation** with `docker --version`
3. **Test Docker** with `docker run hello-world`
4. **Navigate to the repository** and follow the project-specific setup instructions

## Troubleshooting

### Docker Permission Issues
If you get permission denied errors with Docker:
```bash
# Log out and log back in, or run:
newgrp docker
```

### Git LFS Issues
If LFS files aren't downloading properly:
```bash
git lfs pull --all
```

### Nix Issues
If Nix commands aren't found:
```bash
source ~/.profile
# or restart your terminal
```