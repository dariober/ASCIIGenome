class Asciigenome < Formula
  desc "Text Only Genome Viewer!"
  homepage "https://github.com/dariober/ASCIIGenome"
  url "https://github.com/dariober/ASCIIGenome/releases/download/v1.8.0/ASCIIGenome-1.8.0.zip"
  sha256 "bb6ae616604898045ed7dca6d267c506cff9e78968d149344f79d1444d80bd0f"

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
