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
1. [x] TL set up cvae microservice ------------- Jan 23
2. [ ] TL update microservice to return property **titles**, categories, reasons_for_categorization, known_values, predicted_values
3. [ ] TL update microservice to allow requesting properties for a specific category

4. JW get the services iframe working with ports other than 80
5. JW get the user registration working with ports other than 80.
6. JW get the new_report form submission to open the report's `/` route rather than opening a new window

7. JW get report generation working with pending reports
   1. generate a dummy microservice response, this can be anything for now, but work in some random delay
   2. handle long delays (1-2 minutes) for report generation, make the UI not suck for that.
   3. return json version of the report as well

8. get chemical search working ---------- Jan 23
9. deploy to toxindex ------------------- Jan 24
10. test with partners ------------------- Jan 25
