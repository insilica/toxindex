{
  description = "toxindex fullstack devShell with venv and GitHub SSH Python deps";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";
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
            pkgs.python310
            pkgs.uv
            pkgs.git
            pkgs.flyway
            pkgs.postgresql_15
            pkgs.postgresql_jdbc
            pkgs.redis
            pkgs.awscli2
            pkgs.gcc
            pkgs.zlib 
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
              uv venv
              source .venv/bin/activate
              GIT_CLONE_PROTECTION_ACTIVE=false uv pip install -r requirements.txt
            else
              source .venv/bin/activate
            fi
            unset PYTHONPATH
            
            export LD_LIBRARY_PATH="${pkgs.zlib.out}/lib:${pkgs.stdenv.cc.cc.lib}/lib:$LD_LIBRARY_PATH"
            # force eventlet to use the system resolver instead of DNS monkey patching
            export EVENTLET_NO_GREENDNS=yes 
            echo "Python packages in uv venv:"
            uv pip list
          '';
        };
      }
    );
}
