{
  description = "toxindex fullstack devShell with venv and GitHub SSH Python deps";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.05";
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
            pkgs.google-cloud-sdk
            pkgs.gcc
            pkgs.zlib 
            pkgs.nodejs
            pkgs.gh
            pkgs.libjpeg

            # --- Sphinx docs ---
            pkgs.python312Packages.sphinx
            pkgs.python312Packages.myst-parser
            pkgs.python312Packages.furo
            pkgs.python312Packages.sphinx-autobuild
            pkgs.python312Packages.sphinx-rtd-theme
            pkgs.graphviz

            # --- for rd.kit.Chem.Draw (used in SyGMa) ---
            pkgs.xorg.libX11
            pkgs.xorg.libXext
            pkgs.xorg.libSM
            pkgs.xorg.libXrender
            pkgs.expat
          ];


          shellHook = ''
            # Remove NVM from PATH to avoid conflicts
            export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v '\.nvm' | paste -sd:)
            unset NVM_DIR

            # Set locale and library paths
            export LANG=C.UTF-8
            export LC_ALL=C.UTF-8
            # export LD_LIBRARY_PATH="${pkgs.zlib.out}/lib:${pkgs.stdenv.cc.cc.lib}/lib:$LD_LIBRARY_PATH"
            export LD_LIBRARY_PATH="${pkgs.zlib.out}/lib:${pkgs.stdenv.cc.cc.lib}/lib:${pkgs.xorg.libXrender}/lib:${pkgs.xorg.libX11}/lib:${pkgs.xorg.libXext}/lib:${pkgs.xorg.libSM}/lib:${pkgs.expat}/lib:$LD_LIBRARY_PATH"


            # Ensure GitHub CLI is authenticated
            if ! gh auth status >/dev/null 2>&1; then
              echo "Please run: gh auth login"
              exit 1
            fi

            # Ensure we're using the Nix-provided Python
            export PYTHONPATH="$PWD:$PYTHONPATH"
            hash -r

            # Check Python version before proceeding
            echo "Using Python: $(which python3)"
            echo "Python version: $(python3 --version)"

            # Remove existing .venv if it exists
            if [ -d ".venv" ]; then
              echo "Removing existing .venv..."
              rm -rf .venv
            fi

            # Create new venv with Nix Python explicitly
            echo "Creating new virtual environment with Python 3.12..."
            python3 -m venv .venv
            source .venv/bin/activate

            # Verify we're using the right Python
            echo "Virtual environment Python: $(which python)"
            echo "Virtual environment Python version: $(python --version)"
            
            # Sync dependencies with explicit Python version
            echo "Syncing dependencies..."
            GIT_CLONE_PROTECTION_ACTIVE=false uv sync --python python3
            if [ $? -ne 0 ]; then
              echo "Error: 'uv sync' failed. Please check your dependencies."
              exit 1
            fi

            # Install rdkit first (must be before sygma)
            pip install rdkit
            if [ $? -ne 0 ]; then
              echo "Error: Failed to install rdkit."
              exit 1
            fi

            # Install SyGMa (original tool)
            pip install --no-build-isolation sygma
            if [ $? -ne 0 ]; then
              echo "Error: Failed to install sygma."
              exit 1
            fi

            # Install local RAPtool package in editable mode
            echo "Installing RAPtool in editable mode..."
            pip install -e RAPtool/
            if [ $? -ne 0 ]; then
              echo "Error: Failed to install RAPtool. Please check the RAPtool directory."
              exit 1
            fi

            # Install local pathway_analysis_tool package in editable mode
            echo "Installing pathway_analysis_tool in editable mode..."
            pip install -e pathway_analysis_tool/
            if [ $? -ne 0 ]; then
              echo "Error: Failed to install pathway_analysis_tool. Please check the pathway_analysis_tool directory."
              exit 1
            fi

            # Install local metabolite-sygma package in editable mode
            echo "Installing metabolite-sygma in editable mode..."
            pip install -e sygma/metabolite-sygma/
            if [ $? -ne 0 ]; then
              echo "Error: Failed to install metabolite-sygma. Please check the sygma/metabolite-sygma directory."
              exit 1
            fi

            # Source your fullstack environment setup (Postgres, Redis, AWS, etc.)
            # This runs after the virtual environment is set up so Python scripts can access packages
            # source scripts/flake/shellhook.sh

            # Load environment variables from .env file
            if [ -f ".env" ]; then
              echo "Loading environment variables from .env file..."
              export $(grep -v '^#' .env | xargs)
              echo "Environment variables loaded: REDIS_HOST=$REDIS_HOST, REDIS_PORT=$REDIS_PORT"
            else
              echo "Warning: .env file not found. Redis connection may fail."
            fi

            # Unset PYTHONPATH after shellhook runs to avoid conflicts
            unset PYTHONPATH

            # Pathway tool below
            # Run or reuse the Blazegraph container
            #if ! docker ps -a --format '{{.Names}}' | grep -q '^blazegraph$'; then
            #  docker run --name blazegraph \
            #    -d -p 8889:8080 \
            #    lyrasis/blazegraph:2.1.5
            #else
            #  docker start blazegraph 2>/dev/null || true
            #fi

            # Wait for Blazegraph to start listening
            #until curl -s http://localhost:8889/bigdata/namespace/kb/sparql -o /dev/null; do
            #  printf '.'; sleep 1
            #done
            #echo "Blazegraph up—loading TTL files"

            # install npm dependencies in frontend
            # (cd frontend && npm install)

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
