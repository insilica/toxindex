# Toxindex Portal

## Setup

- Install `nix` via package manager or [nix-installer](https://github.com/DeterminateSystems/nix-installer)
- Ensure nix flake support is enabled

  ```sh
  # ~/.config/nix/nix.conf
  experimental-features = nix-command flakes
  ```

- Install `direnv` ([instructions](https://direnv.net/docs/installation.html))
- Enable `.envrc` file for project directory

  ```console
  $ direnv allow
  ## Your shell will now load the project environment automatically
  ```

- Get a copy of the `.env` file from a team member and add it to project root
- Run the project with `devenv up` (configuration in `flake.nix`)
- Go to http://localhost:6513

## Services

1. postgres
2. flyway - runs on startup and performs all necessary migrations on postgres db
3. report - generates pdf reports and returns a report view
4. webserver - a flask web user interface, it has users, projects, and projects have views for each services

## Todo

1. add a postgrest or postgraphile service as a datastore
2. migrate existing functionality to use the datastore
