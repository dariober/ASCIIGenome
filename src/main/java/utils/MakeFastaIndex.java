package utils;

import exceptions.InvalidCommandLineException;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

public class MakeFastaIndex {

  public MakeFastaIndex(String fastaFile) {}

  public List<fastaIndexRecord> makeIndex()
      throws InvalidCommandLineException, FileNotFoundException {

    List<fastaIndexRecord> faidx = new ArrayList<fastaIndexRecord>();

    return faidx;
  }

  private class fastaIndexRecord {}
}
