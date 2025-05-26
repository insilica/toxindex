{
  description = "Dev shell with Flyway, Python, PostgreSQL (TCP/IP focus)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        pgSettings = {
          user = "postgres";
          password = "devpassword";
          port = "5433";
          dbName = "toxindex";
          host = "localhost";
          # Define a local socket directory
          socketDirSubPath = "socket"; # Relative to $PWD/.postgres
        };

        python-with-deps = pkgs.python3.withPackages (ps: with ps; [
          flask
          python-dotenv
          stripe
          flask-login
          flask-wtf
          sendgrid
          gunicorn
          psycopg2-binary
          email-validator
          boto3
          openai
          celery
          redis
          flask-socketio
          eventlet
          pandas
        ]);
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            pkgs.flyway
            python-with-deps
            pkgs.postgresql_15
            pkgs.postgresql_jdbc
            pkgs.redis
            pkgs.awscli2
          ];

          shellHook = ''
            echo "Setting up PostgreSQL development environment..."

            export PGDATA="$PWD/.postgres/data"
            # Define and create the explicit socket directory for the server
            export PG_SOCKET_DIR="$PWD/.postgres/${pgSettings.socketDirSubPath}"

            export PGUSER="${pgSettings.user}"
            export PGPASSWORD="${pgSettings.password}"
            export PGHOST="${pgSettings.host}" # For clients to default to TCP/IP
            export PGPORT="${pgSettings.port}"
            export PGDATABASE="${pgSettings.dbName}"
            export PATH="${pkgs.postgresql_15}/bin:$PATH"
            
            PGLOGFILE="$PWD/.postgres/logfile"

            mkdir -p "$PWD/.postgres" # Ensure .postgres parent directory exists
            mkdir -p "$PG_SOCKET_DIR"  # Ensure the socket directory exists before initdb/start

            if [ ! -f "$PGDATA/postgresql.conf" ]; then
              echo "Initializing new PostgreSQL instance in $PGDATA..."
              mkdir -p "$PGDATA" # Ensure PGDATA directory itself exists
              initdb --username="$PGUSER" --pwfile=<(echo "$PGPASSWORD") -D "$PGDATA"

              echo "Configuring PostgreSQL..."
              echo "listen_addresses = '${pgSettings.host}'" >> "$PGDATA/postgresql.conf"
              echo "port = ${pgSettings.port}" >> "$PGDATA/postgresql.conf"
              # Crucially, set a writable unix_socket_directories
              echo "unix_socket_directories = '$PG_SOCKET_DIR'" >> "$PGDATA/postgresql.conf"
            fi

            if ! pg_ctl -D "$PGDATA" status > /dev/null 2>&1; then
              echo "Starting PostgreSQL server (Log: $PGLOGFILE)..."
              if ! pg_ctl -D "$PGDATA" -l "$PGLOGFILE" -w -t 30 start; then
                echo "--------------------------------------------------------------------"
                echo "PostgreSQL FAILED to start. Contents of $PGLOGFILE:"
                cat "$PGLOGFILE"
                echo "--------------------------------------------------------------------"
                echo "Exiting due to PostgreSQL startup failure."
                exit 1
              fi
              echo "PostgreSQL server started successfully."
            else
              echo "PostgreSQL already running."
            fi

            # Create database if it doesn't exist
            # psql will use PGHOST=localhost and PGPORT for TCP/IP connection
            if ! psql -lqt | cut -d \| -f 1 | grep -qw "${pgSettings.dbName}"; then
              echo "Creating database: ${pgSettings.dbName}..."
              createdb "${pgSettings.dbName}"
            else
              echo "Database ${pgSettings.dbName} already exists."
            fi

            # ... (AWS SSO Configuration - unchanged) ...
             if [ -f .aws-profile ]; then
              export AWS_PROFILE=$(cat .aws-profile)
              echo "Using saved AWS_PROFILE=$AWS_PROFILE"
            else
              echo "No saved AWS_PROFILE found."
              read -p "Enter your AWS SSO profile name: " entered_profile
              export AWS_PROFILE=$entered_profile
              echo "$AWS_PROFILE" > .aws-profile
              echo "Saved AWS_PROFILE to .aws-profile"
            fi

            if ! grep -q "\[profile $AWS_PROFILE\]" ~/.aws/config 2>/dev/null; then
              echo "AWS profile '$AWS_PROFILE' not found in ~/.aws/config. Launching 'aws configure sso'..."
              aws configure sso
            fi

            echo "Triggering AWS SSO login for profile '$AWS_PROFILE'..."
            if ! aws sts get-caller-identity --profile "$AWS_PROFILE" >/dev/null 2>&1; then
              aws sso login --profile "$AWS_PROFILE"
            fi

            echo "Running Flyway migrations..."
            flyway \
              -url="jdbc:postgresql://${pgSettings.host}:${pgSettings.port}/${pgSettings.dbName}" \
              -user="${pgSettings.user}" \
              -password="${pgSettings.password}" \
              -locations="filesystem:./flyway/sql" \
              migrate

            echo "Loading application environment script..."
            bash scripts/load_environment.sh "$AWS_PROFILE" "insilica/toxindex+development"

            echo "Starting Redis server..."
            if ! pgrep redis-server > /dev/null; then
              echo "Starting Redis server..."
              redis-server --daemonize yes
            else
              echo "Redis is already running."
            fi


            echo "Development environment ready."
            echo "PostgreSQL is accessible on: ${pgSettings.host}:${pgSettings.port} (TCP/IP)"
            echo "PostgreSQL server socket directory: $PG_SOCKET_DIR"
            echo "Database: ${pgSettings.dbName}, User: ${pgSettings.user}"

            export FLASK_APP=webserver.app
            export FLASK_ENV=development
            export DEBUG=1 
            export PREFERRED_URL_SCHEME=http
            export SERVER_NAME="${pgSettings.host}:6513" 
            source .env
            export PYTHONPATH="$PWD"
          '';
        };

        apps.flask = flake-utils.lib.mkApp {
          drv = pkgs.writeShellScriptBin "flask-dev" ''
            #!${pkgs.bash}/bin/bash
            set -euo pipefail 

            export FLASK_APP=webserver.app
            export FLASK_ENV=development
            export DEBUG=1 
            export PREFERRED_URL_SCHEME=http
            export SERVER_NAME="${pgSettings.host}:6513" 

            export DB_HOST="${pgSettings.host}" # Connect via TCP
            export DB_PORT="${pgSettings.port}"
            export DB_NAME="${pgSettings.dbName}"
            export DB_USER="${pgSettings.user}"
            export DB_PASSWORD="${pgSettings.password}"

            echo "Starting Flask development server on http://${pgSettings.host}:6513..."
            exec flask run --host=0.0.0.0 --port=6513
          '';
        };
      });
}