# Toxindex Portal
A collection of toxicology relevent web applications.

## Services
1. postgres
2. flyway - runs at the start of docker compose up and performs all necessary migrations on postgres service.
3. webserver - a flask web user interface, it has users, projects, and projects have views for each services
4. report - generates pdf reports and returns a report view
5. search - a simple search app that lets users look up information on chemicals

## Development
update your /etc/hosts file to include {subdomain}.localhost for each service:
```
127.0.0.1       report.localhost
127.0.0.1       search.localhost
```
1. make sure you have `docker --version` 24.0.0 or higher
2. navigate to the services directory `cd services`
3. get `test.env` from a team member and add it to ./services/test.env
3. start the project `docker compose --env-file test.env up`
4. go to localhost:6513

## How to add a service
1. create a new directory in services
2. add a dockerfile
3. add the service in services/docker-compose.yml
4. add the service to the ./services/webserver/webserver/templates/layout.html in the id="services-{{ project.project_id }}" object. 

## TODO
1. [x] get reports downloading
2. [x] deploy to toxindex.com --- feb 15
3. [x] update the report -------- feb 16
5. [x] talk to hartung ---------- feb 19
4. [ ] write phase 2 report ----- feb 20-29



extra
-- get https working
-- microservices should upload their own (json, pdf, and html) reports
-- store .env files in s3
-- store raw data in s3
