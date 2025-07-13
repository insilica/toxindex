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
            pkgs.python312
            pkgs.uv
            pkgs.git
            pkgs.flyway
            pkgs.postgresql_15
            pkgs.postgresql_jdbc
            pkgs.redis
            pkgs.awscli2
            pkgs.gcc
            pkgs.zlib 
            pkgs.nodejs
          ];

          shellHook = ''
            # Unset NVM to avoid conflicts with Nix-provided node
            unset NVM_DIR
            export PATH=$(echo $PATH | tr ':' '\n' | grep -v '\.nvm' | paste -sd: -)
            export LANG=C.UTF-8
            export LC_ALL=C.UTF-8
            export LD_LIBRARY_PATH="${pkgs.zlib.out}/lib:${pkgs.stdenv.cc.cc.lib}/lib:$LD_LIBRARY_PATH"
            hash -r

            export PYTHONPATH=$PWD:$PYTHONPATH

            if [ -f .env ]; then
              echo "Found .env file"
              set -a
              source .env
              set +a
            else
              echo "Warning: .env file not found"
            fi

            # Set up Python environment
            [ ! -d .venv ] && uv venv --python ${pkgs.python312}/bin/python3.12
            source .venv/bin/activate
            
            # Sync dependencies (idempotent)
            GIT_CLONE_PROTECTION_ACTIVE=false uv sync
            
            # force eventlet to use the system resolver instead of DNS monkey patching
            export EVENTLET_NO_GREENDNS=yes 
            echo "Python packages in uv venv:"
            uv pip list

            # Source your fullstack environment setup (Postgres, Redis, AWS, etc.)
            # This runs after the virtual environment is set up so Python scripts can access packages
            source scripts/flake/shellhook.sh

            # Unset PYTHONPATH after shellhook runs to avoid conflicts
            unset PYTHONPATH

            # Pathway tool below
            # Run or reuse the Blazegraph container
            if ! docker ps -a --format '{{.Names}}' | grep -q '^blazegraph$'; then
              docker run --name blazegraph \
                -d -p 8889:8080 \
                lyrasis/blazegraph:2.1.5
            else
              docker start blazegraph 2>/dev/null || true
            fi

            # Wait for Blazegraph to start listening
            until curl -s http://localhost:8889/bigdata/namespace/kb/sparql -o /dev/null; do
              printf '.'; sleep 1
            done
            echo "Blazegraph up—loading TTL files"

            # Add data from WikiPathways/AOP-Wiki
            # curl -X POST \
            #  -H "Content-Type: text/turtle" \
            #  --data-binary @pathway_data/wikipathways.ttl \
            #  http://localhost:8889/bigdata/namespace/kb/sparql
            # curl -X POST \
            #  -H "Content-Type: text/turtle" \
            #  --data-binary @pathway_data/AOPWikiRDF.ttl \
            #  http://localhost:8889/bigdata/namespace/kb/sparql
          '';
        };
      }
    );
}
