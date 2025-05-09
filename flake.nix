{
  description = "Dev shell with Flyway, Python, PostgreSQL";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
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
        ]);
      in {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            pkgs.flyway
            python-with-deps
            pkgs.postgresql_15
          ];
            shellHook = ''
            export PGDATA=$PWD/pgdata
            export PATH=${pkgs.postgresql_15}/bin:$PATH

            if [ ! -d "$PGDATA" ]; then
                echo "Initializing PostgreSQL data directory..."
                pg_ctl init -D "$PGDATA"
            fi

            echo "Starting PostgreSQL..."
            pg_ctl -D "$PGDATA" -o "-p 5432" -w start

            echo "Ensuring devdb exists..."
            createdb devdb || true

            echo "Running Flyway migrations..."
            flyway \
                -url=jdbc:postgresql://localhost:5432/devdb \
                -user=postgres \
                -locations=filesystem:./flyway/sql \
                migrate
            '';
        };
      });
}
