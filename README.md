# Toxindex Portal
A collection of toxicology workflows

## Services
1. postgres
2. flyway - runs at the start of docker compose up and performs all necessary migrations on postgres service.
3. webserver - a flask web user interface, it has users, projects, and projects have views for each services

## Development
1. nix develop
2. flask run --app webserver.app --host=0.0.0.0 --port=6513

# TODO
1. build actual chat page
2. fix sidebar styling
3. add streamlit application
4. 