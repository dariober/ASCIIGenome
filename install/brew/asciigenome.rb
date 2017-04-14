class Asciigenome < Formula
  desc "Text Only Genome Viewer!"
  homepage "https://github.com/dariober/ASCIIGenome"
  url "https://github.com/dariober/ASCIIGenome/releases/download/v1.5.0/ASCIIGenome-1.5.0.zip"
  sha256 "dd15381292ece25d5774429e7316bb004a35f4e8bfb48debff2977257eae7cb5"

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
