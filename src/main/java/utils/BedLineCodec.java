/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2015 creation

*/
package utils;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import java.util.regex.Pattern;

public class BedLineCodec extends AsciiFeatureCodec<BedLine> {
  private Pattern tab = Pattern.compile("[\t]");

  public BedLineCodec() {
    super(BedLine.class);
  }

  @Override
  public BedLine decode(String line) {

    if (line.trim().isEmpty()) {
      return null;
    }
    if (BedLine.isBedHeader(line)) return null;

    String[] tokens = tab.split(line);
    if (tokens.length < 2) return null;

    return new BedLine(tokens);
  }

  @Override
  public Object readActualHeader(LineIterator reader) {
    return null;
  }

  @Override
  public boolean canDecode(final String path) {
    return path.toLowerCase().endsWith(".bed");
  }
}
