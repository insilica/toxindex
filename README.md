# Toxindex Portal

To run this project:

1. make sure you have `docker --version` 24.0.0 or higher
2. navigate to the services directory `cd services`
3. get `test.env` from a team member and add it to ./services/test.env
3. start the project `docker compose --env-file test.env up`
4. go to localhost:6513

## Services
1. postgres
2. flyway - runs at the start of docker compose up and performs all necessary migrations on postgres service.
3. report - generates pdf reports and returns a report view
4. webserver - a flask web user interface, it has users, projects, and projects have views for each services

## Todo
0. get project creation working with postgres and the webserver views
1. correctly list the projects a user has access to in the left sidebar
2. when a user has no projects, populate the content area with a message to click '+ new project'
3. in services/report when a user creates a report, update the relevant postgres tables accordingly