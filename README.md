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
1. Iterate on services/report by:
    1. docker compose up
    2. create an account
    3. create a project
    4. click on the project and then the report link
    5. now you can edit /services/report and see the changes in real time
    6. work on basic report crud, make the table prettier, report names, delete, etc.