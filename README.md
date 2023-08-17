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
1. update flyway with the full db schema
    1. need tables for users, projects, chemicals, sources, source_chemicals, project_sources, reports, and project_reports.

1. create the upload service
    1. upload should take a list of inchi identifiers in a file
    2. it should store the file to s3 and reference it in the sources table
    3. it should parse the inchi identifiers, guarantee they are valid, and store them in the chemicals table and source_chemicals table
    4. it should associate the source with the project in the project_source table
    5. it should have a webservice view for uploading the inchi lines file

2. create the overview service
    1. should have a webservice view that just lists the number of uploaded compounds
    2. should have other stuff, get creative.
