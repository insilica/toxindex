# Use the official Flyway image as a base
FROM flyway/flyway:latest

# Copy the migration scripts to the container
COPY ./sql /flyway/sql

# Default command to run Flyway migrate
CMD ["-url=${FLYWAY_URL}", "-user=${FLYWAY_USER}", "-password=${FLYWAY_PASSWORD}", "-locations=filesystem:/flyway/sql", "migrate"]
