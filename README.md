# Toxindex Portal

A collection of toxicology workflows for research and analysis.

---

## Table of Contents

- [Services](#services)
- [Development Setup](#development-setup)
- [Usage](#usage)
- [Troubleshooting](#troubleshooting)
- [Roadmap](#roadmap)
- [Tech Stack](#tech-stack)
- [Contributing](#contributing)
- [License](#license)

---

## Services

1. **Postgres** – Database service.
2. **Flyway** – Runs at the start of `docker compose up` and performs all necessary migrations on the Postgres service.
3. **Webserver** – Flask web user interface. Features users, projects, and project-specific service views.

---

## Development Setup

```sh
# Clean up previous environments (if needed)
exit
rm -rf .venv
rm -rf ~/.cache/uv
unset LD_LIBRARY_PATH
nix develop
# Update dependencies
nix flake update

# Enter development environment
nix develop

# If c++ binaries not found
unset LD_LIBRARY_PATH && nix develop

# Install frontend dependencies
cd frontend && npm install
```

---

## Usage

### Start Services

- **Production:**

  - connect to ssh
    `ssh kyu-ubuntu30`



  - make swap
sudo dd if=/dev/zero of=/swapfile bs=1M count=1024
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
echo '/swapfile swap swap defaults 0 0' | sudo tee -a /etc/fstab
  - clear RAM
    `sudo sync; sudo echo 3 | sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'`
  - check RAM
    `free -h`
  
  - start React
    `cd toxindex && nix develop && cd frontend && npm install &&  && npm run build`

sudo apt update
sudo apt install nginx

sudo nano /etc/nginx/conf.d/toxindex.conf

server {
    listen 80;
    server_name 18.118.10.140;  # or your EC2 public IP

    root /home/ubuntu/toxindex/frontend/dist;
    index index.html;

    location / {
        try_files $uri /index.html;
    }

    location /api/ {
        proxy_pass http://127.0.0.1:6513/;
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

sudo nginx -t
sudo systemctl reload nginx

sudo chown -R ubuntu:www-data /home/ubuntu/toxindex/frontend/dist
sudo chmod -R 755 /home/ubuntu/toxindex/frontend/dist
sudo chmod -R o+rx /home/ubuntu/toxindex/frontend/dist
sudo chmod o+x /home/ubuntu
sudo systemctl reload nginx

sudo tail -n 50 /var/log/nginx/error.log

  - Start webserver:  
    `gunicorn webserver.app:app --bind 0.0.0.0:8000 --worker-class eventlet`
  - Start Redis listener:  
    `cd toxindex && nix develop && python redis_listener_service.py`
  - Start Celery worker:  
    `cd toxindex && nix develop && celery -A workflows.celery_worker worker --loglevel=info`

terminal #1
gunicorn webserver.app:app --bind 0.0.0.0:8000 --worker-class eventlet

terminal #2
python redis_listener_service.py

terminal #3
celery -A workflows.celery_worker worker --loglevel=info

- **Development:**
  - Start Flask app:  
    `python -m webserver.app`
  - Start Celery worker:  
    `celery -A workflows.celery_worker worker --loglevel=info`
  - Start frontend:  
    `cd frontend && npm run dev`

---

## Troubleshooting

If Gunicorn cannot start (e.g., port 8000 is in use):

```sh
lsof -i :8000
kill <PID>
```

---

## Committing Changes

After making changes:

```sh
git add .
git commit -m "Describe your change"
git push
```

---

## Roadmap

The basic workflow:

1. Users create environments.
2. Users execute workflows in their environments.

### ToxIndex RAP Workflow

- Users can predict hazard with toxtransformer.
- Example: User enters prompt "Is PFOA hepatotoxic?"
- System creates a task, adds it to a table, sends it to the workflow, and creates a webhook for activities/results.

See tasks in `pfoa_tasks.md` for more details.

---

## Tech Stack

- **Frontend:** React, Tailwind CSS, Vite, TypeScript
- **Backend:** Flask, Eventlet, Celery
- **Data:** Langchain, Agno, Scikit-learn, Postgres, Flyway
- **Deployment:** Gunicorn on EC2

---

## Contributing

Please open issues or submit pull requests.

---

## License

TBD



