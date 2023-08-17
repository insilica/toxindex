#!/bin/sh

flyway -url=$FLYWAY_URL -user=$FLYWAY_USER -password=$FLYWAY_PASSWORD -locations=filesystem:/flyway/sql migrate
