{
  description = "mmfqcount - FASTQ read counter";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/release-25.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
  }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs {inherit system;};
    in {
      # 📦 Build dein Rust binary
      packages.default = pkgs.rustPlatform.buildRustPackage {
        pname = "mmfqcount";
        version = "0.1.0";

        src = self;

        cargoLock = {
          lockFile = ./Cargo.lock;
        };

        nativeBuildInputs = with pkgs; [
          pkg-config
        ];

        buildInputs = with pkgs; [
          # falls du später z.B. openssl brauchst:
          # openssl
        ];

        meta = with pkgs.lib; {
          description = "FASTQ read counter tool";
          license = licenses.mit;
          platforms = platforms.linux;
          mainProgram = "mmfqcount";
        };
      };

      # ▶ CLI run ohne install
      apps.default = {
        type = "app";
        program = "${self.packages.${system}.default}/bin/mmfqcount";
      };

      # 🧪 Dev environment (wie venv + pytest + IDE)
      devShells.default = pkgs.mkShell {
        packages = with pkgs; [
          cargo
          rustc
          rustfmt
          clippy
          rust-analyzer
          pkg-config
        ];

        shellHook = ''
          echo "mmfqcount dev shell ready"
        '';
      };
    });
}
