{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.05";
    nixpkgs-unstable.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    systems.url = "github:nix-systems/default";
    # devenv.url = "github:cachix/devenv";
    devenv.url = "github:cachix/devenv/v0.6.3";
  };

  nixConfig = {
    extra-trusted-public-keys =
      "devenv.cachix.org-1:w1cLUi8dv3hnoSPGAuibQv+f9TZLr6cv/Hm9XgU50cw=";
    extra-substituters = "https://devenv.cachix.org";
  };

  outputs = { self, nixpkgs, nixpkgs-unstable, devenv, systems, ... }@inputs:
    let forEachSystem = nixpkgs.lib.genAttrs (import systems);
    in {
      devShells = forEachSystem (system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
          jdk = pkgs.openjdk17_headless;
          env = {
            DB_HOST = "localhost";
            DB_PORT = 5415;
            DB_NAME = "toxindex";
            DB_USER = "postgres";
            WEBSERVER_PORT = 6513;
            REPORTS_PORT = 6515;
            __STDCXX_PATH = "${pkgs.stdenv.cc.cc.lib}/lib";
            WEBSERVER_SERVER_NAME = "toxindex.com";
            WEBSERVER_PREFERRED_URL_SCHEME = "https";
          };
        in with env; {
          default = devenv.lib.mkShell {
            inherit inputs pkgs;

            modules = [
              # https://devenv.sh/reference/options/
              {
                inherit env;

                packages = with pkgs; [
                  ## dev
                  bun
                  nodePackages.npm
                  nodejs
                  poetry
                  ## db
                  (flyway.override { jre_headless = jdk; })
                  ## reports
                  wkhtmltopdf-bin
                  ## tools
                  gnugrep
                  gnused
                  fd
                  git
                  nixfmt
                  python311Packages.grip # Markdown previewer
                ];

                languages = {
                  python.enable = true;
                  python.package = pkgs.python311;
                  nix.enable = true;
                  java.enable = true;
                  java.jdk.package = jdk;
                };

                services = {
                  postgres = {
                    enable = true;
                    package = pkgs.postgresql_15;
                    initialDatabases = [{ name = DB_NAME; }];
                    initialScript = "CREATE USER ${DB_USER} SUPERUSER;";
                    port = DB_PORT;
                    listen_addresses = DB_HOST;
                  };
                };

                process.implementation = "process-compose";
                process.process-compose = { port = 9999; };

                scripts = { };

                processes = let
                  env_ready = {
                    check-secrets = {
                      condition = "process_completed_successfully";
                    };
                  };
                  db_ready = {
                    postgres = { condition = "process_started"; };
                    flyway = { condition = "process_completed_successfully"; };
                  };
                in {
                  check-secrets = { exec = "check-secrets"; };

                  flyway = {
                    exec = "flyway migrate";
                    process-compose = {
                      depends_on = {
                        postgres = { condition = "process_started"; };
                      };
                    };
                  };

                  webserver = {
                    exec = "./services/webserver/run.sh";
                    process-compose = { depends_on = db_ready // env_ready; };
                  };

                  report = {
                    exec = "./services/report/run.sh";
                    process-compose = { depends_on = db_ready // env_ready; };
                  };

                };

                pre-commit.hooks = {
                  editorconfig-checker.enable = true;
                  shellcheck.enable = true;
                  shfmt.enable = true;
                  # markdownlint.enable = true;
                  # mdsh.enable = true;
                  nixfmt.enable = true;
                  nil.enable = true;
                  black.enable = true;
                  yamllint.enable = true;
                };
                pre-commit.settings = { nixfmt.width = 80; };

                enterShell = "";
              }
            ];
          };
        });
    };
}
