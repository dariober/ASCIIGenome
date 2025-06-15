class Asciigenome < Formula
  desc "Text Only Genome Viewer!"
  homepage "https://github.com/dariober/ASCIIGenome"
  url "https://github.com/dariober/ASCIIGenome/releases/download/v1.19.0/ASCIIGenome-1.19.0.zip"
  sha256 "9f520e4cc47a6fcd9ee48c98488dbd912e5dbdcdfe93d998333ec4bc2c4ffb03"

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
