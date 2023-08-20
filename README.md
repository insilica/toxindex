# Toxindex Portal

To run this project:

1. make sure you have `docker --version` 24.0.0 or higher
2. navigate to the services directory `cd services`
3. start the project `docker compose up`
4. go to localhost:6513

## Services
1. postgres
2. flyway - runs at the start of docker compose up and performs all necessary migrations on postgres service.
3. report - generates pdf reports and returns a report view
4. webserver - a flask web user interface, it has users, projects, and projects have views for each services

## Todo
0. better portal for login/registration

1. update flyway with the full db schema
    1. need tables for users, projects, reports

2. create account page

3. create project settings page 
    1. invite other users to project

4. create actual reports
