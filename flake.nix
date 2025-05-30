{
  description = "toxindex fullstack devShell with venv and GitHub SSH Python deps";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            pkgs.python3
            pkgs.git
            pkgs.python3.pkgs.venvShellHook
            pkgs.flyway
            pkgs.postgresql_15
            pkgs.postgresql_jdbc
            pkgs.redis
            pkgs.awscli2
            pkgs.gcc
          ];

          shellHook = ''
            # Source your fullstack environment setup (Postgres, Redis, AWS, etc.)
            source scripts/flake/shellhook.sh

            export PYTHONPATH=$PWD:$PYTHONPATH

            if [ -f .env ]; then
              echo "Found .env file"
              set -a
              source .env
              set +a
            else
              echo "Warning: .env file not found"
            fi

            if [ ! -d .venv ]; then
              echo "Creating local venv and installing pip packages..."
              python -m venv .venv
              . .venv/bin/activate
              unset PYTHONPATH
              export LD_LIBRARY_PATH="${pkgs.stdenv.cc.cc.lib}/lib:$LD_LIBRARY_PATH"
              pip install --upgrade pip
              GIT_CLONE_PROTECTION_ACTIVE=false pip install -r requirements.txt
            else
              echo "Using existing .venv"
              . .venv/bin/activate
              unset PYTHONPATH
              export LD_LIBRARY_PATH="${pkgs.stdenv.cc.cc.lib}/lib:$LD_LIBRARY_PATH"
            fi

            echo "Python packages in venv:"
            pip list
          '';
        };
      }
    );
}
