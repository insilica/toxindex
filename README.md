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
- Start Celery worker:
  ```sh
  celery -A workflows.celery_worker worker --loglevel=info
  ```
- Start frontend:
  ```sh
  cd frontend && npm run dev
  ```

---

## Production Setup

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
  sudo chown -R ubuntu:www-data /home/ubuntu/toxindex/frontend/dist
  sudo chmod -R 755 /home/ubuntu/toxindex/frontend/dist
  sudo chmod -R o+rx /home/ubuntu/toxindex/frontend/dist
  sudo chmod o+x /home/ubuntu
  ```

- Reload Nginx:
  ```sh
  sudo systemctl reload nginx
  ```

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


