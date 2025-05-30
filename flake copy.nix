{
  description = "Dev shell with Flyway, Python, PostgreSQL (TCP/IP focus)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.11";
    flake-utils.url = "github:numtide/flake-utils";
    raptool = {
      url = "git+ssh://git@github.com/toxindex/RAPtool.git";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    rap = {
      url = "git+ssh://git@github.com/toxindex/RAP.git
      inputs.nixpkgs.follows = "nixpkgs";
    }
  };

  outputs = { self, nixpkgs, flake-utils, raptool }:
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

        # Main Python environment with all required packages
        python-with-deps = pkgs.python3.withPackages (ps: with ps; [
          # Web server packages
          flask
          python-dotenv
          stripe
          requests
          flask-login
          flask-wtf
          sendgrid
          gunicorn
          psycopg2
          email-validator
          boto3
          openai
          celery
          redis
          flask-socketio
          gevent
          gevent-websocket
          pandas
          annotated-types

          # RAPtool packages
          rdkit
          rapidfuzz
          aiohttp
          tenacity
          pyarrow
          fastparquet
          seaborn
          matplotlib
          scikit-learn
          cython
          blosc2

          # Additional packages
          (ps.buildPythonPackage {
            pname = "pubchempy";
            version = "1.0.4";
            src = pkgs.fetchurl {
              url = "https://files.pythonhosted.org/packages/aa/fb/8de3aa9804b614dbc8dc5c16ed061d819cc360e0ddecda3dcd01c1552339/PubChemPy-1.0.4.tar.gz";
              sha256 = "sha256-JOncL8kKsVOydkv4BeUQsUEHAIhPrwUQqefPDWHY7Q4=";
            };
            propagatedBuildInputs = [ ps.requests ];
            doCheck = false;  # Skip tests
          })

          # Add RAPtool as a Python package
          (ps.buildPythonPackage {
            pname = "raptool";
            version = "0.1.0";
            src = raptool.outPath;  # Use the actual source path
            propagatedBuildInputs = [
              ps.rdkit
              ps.rapidfuzz
              ps.aiohttp
              ps.tenacity
              ps.pyarrow
              ps.fastparquet
              ps.seaborn
              ps.matplotlib
              ps.scikit-learn
              ps.cython
              ps.blosc2
            ];
            doCheck = false;  # Skip tests
            format = "setuptools";
            # Install in development mode to ensure all submodules are included
            installPhase = ''
              mkdir -p $out/lib/python3.11/site-packages
              cp -r $src/raptool $out/lib/python3.11/site-packages/
            '';
          })
        ]);

        # Create a development shell with both Python and RAPtool
        devShell = pkgs.mkShell {
          buildInputs = [
            pkgs.flyway
            python-with-deps
            pkgs.postgresql_15
            pkgs.postgresql_jdbc
            pkgs.redis
            pkgs.awscli2
            pkgs.python310Packages.pip
            pkgs.gcc
            pkgs.git
            inputs.rap.packages.${system}.default
          ];

          shellHook = ''
            source scripts/flake/shellhook.sh
            export PYTHONPATH="$PWD:/nix/store/kp1122p01gg4nrxpkizqd8f7ffhd779h-python3-3.11.8-env/lib/python3.11/site-packages:$PYTHONPATH"
          '';
        };
      in
      {
        devShells.default = devShell;

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

        # Expose the Python environment as a package
        packages.default = python-with-deps;
      });
}
