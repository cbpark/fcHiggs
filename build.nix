{ pkgs ? (import <nixpkgs> {}), source ? ./., version ? "dev" }:

with pkgs;

stdenv.mkDerivation {
  name = "fcHiggs-${version}";
  src = lib.cleanSource source;
  enableParallelBuilding = true;
  nativeBuildInputs = [ lhapdf ];
  buildInputs = [ lhapdf.pdf_sets.NNPDF23_lo_as_0130_qed ];

  installPhase = ''
    mkdir -p $out/{bin,lib}
    for destdir in bin lib; do
      cp "$destdir"/* $out/$destdir
    done
  '';
}
