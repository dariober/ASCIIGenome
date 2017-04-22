class Asciigenome < Formula
  desc "Text Only Genome Viewer!"
  homepage "https://github.com/dariober/ASCIIGenome"
  url "https://github.com/dariober/ASCIIGenome/releases/download/v1.6.0/ASCIIGenome-1.6.0.zip"
  sha256 "d0d989590f3dc9dcb3a0f17f1f99eadb2d507a4f5f79190b9b0ef20ab8cf4ec9"

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
