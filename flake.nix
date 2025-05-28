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

          shellHook = ''source scripts/flake/shellhook.sh'';
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
