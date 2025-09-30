---
title: Inspecting Cloud SQL DB Values
sidebar_position: 1
---
**Inspecting DB**

curl -o cloud-sql-proxy https://storage.googleapis.com/cloud-sql-connectors/cloud-sql-proxy/v2.18.0/cloud-sql-proxy.linux.amd64

chmod +x cloud-sql-proxy

./cloud-sql-proxy toxindex:us-east4:toxindex

**Open a new terminal**

sudo apt-get update

sudo apt-get install postgresql-client

psql -h 127.0.0.1 -U postgres -d toxindex

type DB password - ask Kyu

**Now you are in DB** 

- Example) show all messages for a specific task

SELECT * FROM messages WHERE task_id = 'ddd08d18-082a-4404-a03d-6e40532300e3' ORDER BY created_at ASC;