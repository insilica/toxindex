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
            for data_source in wikipathways AOPWikiRDF; do
              curl -X POST \
                -H "Content-Type: text/turtle" \
                --data-binary @pathway_data/${data_source}.ttl \
                http://localhost:8889/bigdata/namespace/kb/sparql
            done
          '';
        };
      }
    );
}
