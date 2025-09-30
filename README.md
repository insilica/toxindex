# Toxindex Portal

A collection of toxicology workflows for research and analysis. This portal provides tools and workflows for toxicology research, including hazard prediction, data management, and more.

---

## Table of Contents

- [Project Overview](#project-overview)
- [Services](#services)
- [Development Setup](#development-setup)
- [Production Setup](#production-setup)
- [Troubleshooting](#troubleshooting)
- [Roadmap](#roadmap)
- [Tech Stack](#tech-stack)
- [Contributing](#contributing)
- [License](#license)

---

## Project Overview

Toxindex Portal is designed to streamline toxicology research by providing a unified interface for running workflows, managing data, and visualizing results. It integrates modern web technologies with robust backend services for a seamless research experience.

---

## Services

| Service     | Description                                                      |
|-------------|------------------------------------------------------------------|
| **Postgres**| Database service.                                                |
| **Flyway**  | Runs migrations on the Postgres service at startup.              |
| **Webserver**| Flask web UI for users, projects, and workflow management.      |

---

## Development Setup

add your user to the docker group
sudo usermod -aG docker $USER
and reboot to take it effect.
type "groups" in terminal and see if "docker" is in the list

```sh
# Clean up previous environments (if needed)
rm -rf .venv ~/.cache/uv
unset LD_LIBRARY_PATH
```

```sh
# Enter development environment
nix develop
```

```sh
# Install frontend dependencies
cd frontend && npm install
```

**Development Commands:**

- Start Flask app:
  ```sh
  python -m webserver.app
  ```
- Start redis listener:
  ```sh
  python redis_listener_standalone.py
  ```
- Start Celery worker:
  ```sh
  celery -A workflows.celery_worker worker --loglevel=info -Q celery,probra,openai,raptool,pathway

  celery -A workflows.celery_worker_probra worker --loglevel=info -Q probra
  celery -A workflows.celery_worker_raptool worker --loglevel=info -Q raptool
  celery -A workflows.celery_worker_ranking worker --loglevel=info -Q ranking
  celery -A workflows.celery_worker_plain worker --loglevel=info -Q openai
  ```
- Start frontend:
  ```sh
  cd frontend && npm run dev
  ```

- Run Redis Listener
``` 
python redis_listener_standalone.py
```
=======


## Production Setup (not for local developlment)

### 1. Connect to Server
```sh
ssh kyu-ubuntu30
```

### 2. Expand RAM with Swap (if needed)
```sh
sudo dd if=/dev/zero of=/swapfile bs=1M count=1024
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
echo '/swapfile swap swap defaults 0 0' | sudo tee -a /etc/fstab
```

### 3. RAM Management
```sh
# Clear RAM cache
sudo sync; sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'

# Check RAM usage
free -h
```

### 4. Enter Nix Environment
```sh
cd toxindex && nix develop
```

### 5. Install and Build Frontend
```sh
cd frontend && npm install
npm run build
```

### 6. Install and Configure Nginx
```sh
sudo apt update
sudo apt install nginx
```

- Create Nginx config file:
  ```sh
  sudo nano /etc/nginx/conf.d/toxindex.conf
  ```
  
  Example config:
  ```nginx
  server {
      listen 80;
      server_name 18.118.10.140; # Use 'localhost' for dev

      root /home/ubuntu/toxindex/frontend/dist;
      index index.html;

      location / {
          try_files $uri /index.html;
      }

      location /api/ {
          proxy_pass http://127.0.0.1:6513;
          proxy_set_header Host $host;
          proxy_set_header X-Real-IP $remote_addr;
          proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
          proxy_set_header X-Forwarded-Proto $scheme;
      }

      location /socket.io/ {
          proxy_pass http://127.0.0.1:6513;
          proxy_http_version 1.1;
          proxy_set_header Upgrade $http_upgrade;
          proxy_set_header Connection "upgrade";
          proxy_set_header Host $host;
          proxy_set_header X-Real-IP $remote_addr;
          proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
          proxy_set_header X-Forwarded-Proto $scheme;
      }
  }
  ```

- Test config:
  ```sh
  sudo nginx -t
  ```

- Set permissions to expose frontend static files to Nginx:
  ```sh
  sudo chown -R ubuntu:www-data /home/kyu/toxindex/frontend/dist
  sudo chmod -R 755 /home/kyu/toxindex/frontend/dist
  sudo chmod -R o+rx /home/kyu/toxindex/frontend/dist
  sudo chmod o+x /home/kyu
  ```

- Reload Nginx:
  ```sh
  sudo systemctl reload nginx
  ```

 - Create a systemd Service File for backend
sudo nano /etc/systemd/system/toxindex-backend.service

[Unit]
Description=Toxindex Backend Service
After=network.target

[Service]
# Set your user and working directory
User=kyu
WorkingDirectory=/home/kyu/toxindex

# Activate nix develop environment and run your backend
ExecStart=/usr/bin/bash -c 'source /home/kyu/toxindex/.venv/bin/activate && python -m webserver.app'

# Restart on failure
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target




- Create a systemd Service File for celery worker
sudo nano /etc/systemd/system/toxindex-celery.service

[Unit]
Description=Toxindex Celery Worker
After=network.target

[Service]
User=kyu
WorkingDirectory=/home/kyu/toxindex
Environment="PATH=/home/kyu/.nix-profile/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ExecStart=/usr/bin/bash -c 'cd /home/kyu/toxindex && nix develop --command celery -A workflows.celery_worker worker --loglevel=info'
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target

Reload systemd
sudo systemctl daemon-reload

Start the Service
sudo systemctl start toxindex-backend
sudo systemctl start toxindex-celery

Enable the Service at Boot
sudo systemctl enable toxindex-backend
sudo systemctl enable toxindex-celery

Check Status and Logs

Check status:
sudo systemctl status toxindex-backend
sudo systemctl status toxindex-celery

View logs:
journalctl -u toxindex-backend -f
journalctl -u toxindex-celery -f
---

## Troubleshooting

- If a port is in use:
  ```sh
  lsof -i :8000
  kill <PID>
  ```

---

## Roadmap

The basic workflow:

1. Users create environments.
2. Users execute workflows in their environments.

### ToxIndex RAP Workflow

- Users can predict hazard with toxtransformer.
- Example: User enters prompt: _"Is PFOA hepatotoxic?"_
- System creates a task, adds it to a table, sends it to the workflow, and creates a webhook for activities/results.

---

## Tech Stack

- **Frontend:** React, Tailwind CSS, Vite, TypeScript
- **Backend:** Flask, Eventlet, Celery, Postgres, Flyway
- **Data:** Langchain, Agno, Scikit-learn

---

## Contributing

Contributions are welcome! Please open issues or submit pull requests for improvements.

---

## License

Distributed under the MIT License. See `LICENSE` for more information.


# CI/CD Test - Mon Jul 28 21:57:31 UTC 2025
