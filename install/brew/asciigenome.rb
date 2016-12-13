class Asciigenome < Formula
  desc "Text Only Genome Viewer!"
  homepage "https://github.com/dariober/ASCIIGenome"
  url "https://github.com/dariober/ASCIIGenome/releases/download/v0.1.0/ASCIIGenome-0.6.4.zip"
  sha256 "687ea8b8428ceef49b5a25f0516f9146f08614191f8eb041c577b142e1c91f24"

  def install
    jar = "ASCIIGenome.jar"
    java = share/"java"
    java.install jar
    bin.write_jar_script java/jar, "ASCIIGenome"

    inreplace("ASCIIGenome", /^prefix.*$/, "prefix=#{java}")
    bin.install %w[ASCIIGenome]
  end

  test do
    system "ASCIIGenome", "--version"
  end
end
