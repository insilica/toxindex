---
title: creating a GCP VM
sidebar_position: 2
---

Create a GCP VM in toxindex project

## 

1. Ask kyu for compute admin role in GCP toxindex project
2. Create a VM on console: https://console.cloud.google.com/compute/instancesAdd?authuser=1&invt=Ab6AjA&project=toxindex

## Settings for the VM

1. Set Debian OS 
2. Set machine type: e2-standard-4 to begin with 60GB disk space
3. Set Access scopes to “Allow full access to all Cloud APIs”
4. Set Service account for the VM as “gcs-access” to gain all permissions for future GCP API usages

## Connecting to VM

1. In your local terminal, install gcloud cli in order to connect to GCP VM
2. Sign into GCP using gcloud auth
3. Go to GCP console > VM instances > look for external IP address of the VM you created
4. Add the IP to .ssh/config
5. In your local terminal, run: 

```bash
ssh [vm name]
```

1. In cursor or VScode, install remote-ssh extension
2. connect to VM in your IDE

### Set static IP to avoid change in ip from restart (optional)

The IP of VM changes every restart. You may lose connection to GKE in that case until you add your IP again to GKE.

gcloud compute addresses create toxindex-dev-static-ip --region=us-east4

gcloud compute instances delete-access-config toxindex-dev --zone=us-east4-a --access-config-name="External NAT”