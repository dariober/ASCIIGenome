package utils;


import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

import exceptions.InvalidCommandLineException;

public class MakeFastaIndex {

    public MakeFastaIndex(String fastaFile) {
    }
    
    public List<fastaIndexRecord> makeIndex() throws InvalidCommandLineException, FileNotFoundException {
        
        List<fastaIndexRecord> faidx= new ArrayList<fastaIndexRecord>();
        
        return faidx;
    }

    private class fastaIndexRecord {
    }
    
}
