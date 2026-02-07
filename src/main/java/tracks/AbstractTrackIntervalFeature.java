package tracks;

import colouring.Xterm256;
import com.google.common.collect.Lists;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.lang3.StringUtils;
import org.broad.igv.bbfile.BBFileReader;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.*;
import java.util.regex.Pattern;

public abstract class AbstractTrackIntervalFeature<T extends IntervalFeature>
        extends AbstractTrack {
    /** For GTF/GFF data: Use this attribute to get the feature names */
    protected TabixReader tabixReader; // Leave *protected* for TrackBookmark to work
    protected BBFileReader bigBedReader;
    protected List<T> intervalFeatureList = new ArrayList<>();

    protected abstract T createFeature(String line)
            throws InvalidGenomicCoordsException;


    @Override
    /**
     * Collect features mapping to the current genomic coordinates and update the list of interval
     * features for the track. Also update the mapping of features to the terminal screen. This method
     * should be called only from within methods that change the features being displayed. E.g.
     * setGc(), which changes the coordinates, or setHideRegex() & setShowRegex() which change the
     * visible features. update() should not change anything other than the list of features and the
     * mapping.
     */
    public void update()
            throws IOException,
            InvalidGenomicCoordsException,
            ClassNotFoundException,
            InvalidRecordException,
            SQLException {

        this.intervalFeatureList =
                this.getFeaturesInInterval(
                        this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());
        for (T ift : this.intervalFeatureList) {
            ift.mapToScreen(this.getGc().getMapping());
        }
    }

    protected List<T> getFeaturesInInterval(String chrom, int from, int to)
            throws IOException, InvalidGenomicCoordsException {
        if (from < 1) {
            System.err.println("from < 1: " + from + "; resetting to 1.");
            from = 1;
        }
        if (from > to) {
            System.err.println(
                    "Invalid coordinates: from: "
                            + from
                            + "; to: "
                            + to
                            + "; Resetting to initial 1-"
                            + Integer.MAX_VALUE);
            throw new InvalidGenomicCoordsException();
        }
        List<T> xFeatures = new ArrayList<>();
        this.removeInvisibleFeatures(xFeatures);
        return xFeatures;
    }

    protected void removeInvisibleFeatures(List<T> iftList)
            throws InvalidGenomicCoordsException, IOException {

        for (int i = 0; i < iftList.size(); i++) {

            T tr = iftList.get(i);
            String x = tr.getRaw();

            boolean showIt = true;
            if (this.getShowRegex() != null
                    && !this.getShowRegex().pattern().equals(Filter.DEFAULT_SHOW_REGEX.getValue())) {
                showIt = this.getShowRegex().matcher(x).find();
            }

            boolean hideIt = false;
            if (this.getHideRegex() != null
                    && !this.getHideRegex().pattern().isEmpty()
                    && !this.getHideRegex().pattern().equals(Filter.DEFAULT_HIDE_REGEX.getValue())) {
                hideIt = this.getHideRegex().matcher(x).find();
            }
            if (!showIt || hideIt) {
                iftList.set(i, null);
            }
        }
        iftList.removeIf(Objects::isNull);

        if (!this.getAwk().isEmpty()) {
            this.getAwk();
            String[] rawLines = new String[iftList.size()];
            for (int i = 0; i < iftList.size(); i++) {
                rawLines[i] = iftList.get(i).getRaw();
            }
            boolean[] passAwk;

            // Awk
            try {
                passAwk = Utils.passAwkFilter(rawLines, this.getAwk());
            } catch (Exception e) {
                System.err.print(
                        Utils.padEndMultiLine(
                                "Error processing awk script.", this.getGc().getUserWindowSize()));
                try {
                    this.setAwk("");
                } catch (ClassNotFoundException | InvalidRecordException | SQLException e1) {
                    e1.printStackTrace();
                }
                throw new IOException();
            }
            for (int i = 0; i < passAwk.length; i++) {
                if (!passAwk[i]) {
                    iftList.set(i, null);
                }
            }
            iftList.removeIf(Objects::isNull);
        }
    }

    /**
     * Return the coordinates of the next feature so that the start coincide with the start of the
     * feature and the end is the start + windowSize.
     */
    public GenomicCoords coordsOfNextFeature(GenomicCoords currentGc, boolean getPrevious)
            throws InvalidGenomicCoordsException, IOException {

        T nextFeature;
        if (getPrevious) {
            nextFeature = this.getPreviousFeature(currentGc.getChrom(), currentGc.getFrom());
        } else {
            nextFeature = this.getNextFeature(currentGc.getChrom(), currentGc.getTo());
        }
        if (nextFeature == null) {
            return currentGc;
        }
        GenomicCoords nextGc =
                new GenomicCoords(
                        Utils.coordinatesToString(
                                nextFeature.getChrom(),
                                nextFeature.getFrom(),
                                nextFeature.getFrom() + currentGc.getGenomicWindowSize() - 1),
                        currentGc.getUserWindowSize(),
                        currentGc.getSamSeqDict(),
                        currentGc.getFastaFile());
        return nextGc;
    }

    protected T getNextFeature(String startChrom, int from)
            throws IOException, InvalidGenomicCoordsException {

        T next = this.getNextFeatureOnChrom(startChrom, from);
        if (next != null) {
            return next;
        } // There is no feature left on the starting chrom at the given position.

        // Search the remaining chroms in order starting from position 1:
        List<String> chroms = this.getChromListStartingAt(startChrom);
        if (chroms.contains(startChrom)) {
            // The start chroms is searched last. The `if` statement controls for the case where
            // the startChrom is not present at all in this track.
            chroms.remove(startChrom);
            chroms.add(startChrom);
        }
        for (String chrom : chroms) {
            next =
                    this.getNextFeatureOnChrom(
                            chrom, 0); // Use 0 so if the next feature starts at the beginning of the chrom,
            // i.e. at start=1, it is not missed. See issue #50
            if (next != null) {
                return next;
            }
        }
        return null;
    }

    /**
     * Get the next feature on chrom after "from" position or null if no feature found. This function
     * should be used only by getNextFeature() which is more general since it searches all the
     * chromosomes in turn.
     *
     * @throws IOException
     * @throws InvalidGenomicCoordsException
     */
    protected T getNextFeatureOnChrom(String chrom, int from)
            throws IOException, InvalidGenomicCoordsException {

        int qend = Math.max((from - 1), 0);

        TabixBigBedIterator iter = this.getReader().query(chrom, qend, Integer.MAX_VALUE);

        // We accumulate here a number of feature to test for visibility and we return the first one
        // that passes. Testing one by one would be slow.
        List<T> buffer = new ArrayList<>();

        while (true) {
            String line = iter.next();
            if (buffer.size() < 1000 && line != null) {
                T x = this.createFeature(line);
                buffer.add(x);
                continue;
            }
            if (!buffer.isEmpty()) {
                this.removeInvisibleFeatures(buffer);
                for (int i = 0; i < buffer.size(); i++) {
                    if (buffer.get(i).getFrom() > from) {
                        return buffer.get(i);
                    }
                }
                buffer.clear();
            }
            if (line == null) {
                return null;
            }
        }
    }

    private T getPreviousFeature(String startChrom, int pos)
            throws IOException, InvalidGenomicCoordsException {

        T prev = getPreviousFeatureOnChrom(startChrom, pos);
        if (prev != null) {
            return prev;
        }
        // There is no feature left on the starting chrom at the given position.

        // Search the remaining chroms in order starting from last:
        List<String> chroms = Lists.reverse(this.getChromListStartingAt(startChrom));
        for (String chrom : chroms) {
            prev = getPreviousFeatureOnChrom(chrom, Integer.MAX_VALUE);
            if (prev != null) {
                return prev;
            }
        }
        return null;
    }

    /**
     * Experimental method Find the previous feature relative to the given chrom-pos position. Read
     * the chunk of file *before* the pos position and return the last feature found. I.e. find the
     * feature immediately to to left of pos and not overlapping pos.
     */
    private T getPreviousFeatureOnChrom(String chrom, int pos)
            throws IOException, InvalidGenomicCoordsException {

        int chunkFrom = Math.max(0, pos - 1000000);
        int chunkTo = pos - 1; // -1 because we don't include the current position
        T last = null;

        while (chunkTo > 0) {

            TabixBigBedIterator iter = this.getReader().query(chrom, chunkFrom, chunkTo);

            List<T> buffer = new ArrayList<T>();

            while (true) {
                // Find the last feature in this chunk where end coordinate is less then pos
                String line = iter.next();
                if (buffer.size() < 1000 && line != null) {
                    T x = createFeature(line);
                    buffer.add(x);
                    continue;
                }
                if (!buffer.isEmpty()) {
                    this.removeInvisibleFeatures(buffer);
                    for (T intervalFeature : buffer) {
                        if (intervalFeature.getTo() < pos) {
                            last = intervalFeature;
                        }
                    }
                    buffer.clear();
                }
                if (line == null) {
                    break;
                }
            }
            if (last != null) {
                break; // The last feature is not null and valid. Stop looking
            } else {
                // Move to previous chunk;
                chunkTo = chunkFrom;
                chunkFrom = Math.max(0, chunkFrom - 1000000);
            }
        }
        if (last == null) {
            return null;
        }
        List<T> xLast = new ArrayList<>();
        xLast.add(last);
        this.removeInvisibleFeatures(xLast);
        if (last.getTo() > pos || xLast.isEmpty()) {
            return null;
        }
        return last;
    }

    protected GenomicCoords startEndOfNextFeature(GenomicCoords currentGc, boolean getPrevious)
            throws InvalidGenomicCoordsException, IOException {
        T nextFeature;
        if (getPrevious) {
            nextFeature = getPreviousFeature(currentGc.getChrom(), currentGc.getFrom());
        } else {
            nextFeature = getNextFeature(currentGc.getChrom(), currentGc.getTo());
        }
        if (nextFeature == null) {
            return currentGc;
        }
        GenomicCoords nextGc =
                new GenomicCoords(
                        Utils.coordinatesToString(
                                nextFeature.getChrom(), nextFeature.getFrom(), nextFeature.getTo()),
                        currentGc.getUserWindowSize(),
                        currentGc.getSamSeqDict(),
                        currentGc.getFastaFile());
        return nextGc;
    }

    public GenomicCoords findNextMatch(GenomicCoords currentGc, Pattern pattern)
            throws IOException, InvalidGenomicCoordsException {

        T nextFeature =
                findNextRegexInGenome(pattern, currentGc.getChrom(), currentGc.getTo());
        if (nextFeature == null) {
            return currentGc;
        }
        GenomicCoords nextGc =
                new GenomicCoords(
                        Utils.coordinatesToString(
                                nextFeature.getChrom(),
                                nextFeature.getFrom(),
                                nextFeature.getFrom() + currentGc.getGenomicWindowSize() - 1),
                        currentGc.getUserWindowSize(),
                        currentGc.getSamSeqDict(),
                        currentGc.getFastaFile());
        return nextGc;
    }

    /**
     * Execute findAllChromRegexInGenome() and return the extreme coordinates of the matched features
     */
    protected GenomicCoords genomicCoordsAllChromMatchInGenome(
            Pattern pattern, GenomicCoords currentGc) throws IOException, InvalidGenomicCoordsException {

        List<T> matchedFeatures = this.findAllChromMatchInGenome(pattern, currentGc);

        if (matchedFeatures.isEmpty()) {
            return currentGc;
        }

        // Now get the coords of the first and last feature matched.
        String chrom = matchedFeatures.get(0).getChrom();
        int startFrom = matchedFeatures.get(0).getFrom();
        int endTo = matchedFeatures.get(matchedFeatures.size() - 1).getTo();
        GenomicCoords allMatchesGc =
                new GenomicCoords(
                        Utils.coordinatesToString(chrom, startFrom, endTo),
                        currentGc.getUserWindowSize(),
                        currentGc.getSamSeqDict(),
                        currentGc.getFastaFile());
        return allMatchesGc;
    }

    /**
     * Find all the feature matching regex. Only the feature on one chromosome are returned and this
     * chromsome is the first one to have a match. The search starts from the beginning of the current
     * chrom and if nothing is found continues to the other chroms.
     *
     * @throws InvalidGenomicCoordsException
     */
    private List<T> findAllChromMatchInGenome(Pattern pattern, GenomicCoords currentGc)
            throws IOException, InvalidGenomicCoordsException {

        // Accumulate features here
        List<T> matchedFeatures = new ArrayList<T>();

        // We start search from input chrom
        List<String> chromSearchOrder = null;
        chromSearchOrder = getChromListStartingAt(currentGc.getChrom());

        chromSearchOrder.add(currentGc.getChrom());
        for (String curChrom : chromSearchOrder) {

            TabixBigBedIterator iter = this.getReader().query(curChrom, 0, Integer.MAX_VALUE);

            while (true) {
                String line = iter.next();
                if (line == null) {
                    break;
                }
                boolean matched = pattern.matcher(line).find();
                if (matched) {
                    T x = this.createFeature(line);
                    matchedFeatures.add(x);
                    // if(this.featureIsVisible(x.getRaw())){
                    //    matchedFeatures.add(x);
                    // }
                }
            }
            // Now you have all the features matching regex on this chromosome. We need to exclude those
            // that are not visible
            this.removeInvisibleFeatures(matchedFeatures);

            if (!matchedFeatures.isEmpty()) {
                // At least one feature matching regex found on this chrom.
                // Check we are at the same position as the beginning. if so, continue to other chroms
                if (matchedFeatures.get(0).getChrom().equals(currentGc.getChrom())
                        && matchedFeatures.get(0).getFrom() == currentGc.getFrom()
                        && matchedFeatures.get(matchedFeatures.size() - 1).getTo() == currentGc.getTo()) {
                    // Discard results and keep searching other chroms.
                    matchedFeatures = new ArrayList<>();
                } else {
                    break;
                }
            }
        } // Loop chrom
        return matchedFeatures;
    }

    @Override
    public String printToScreen() throws InvalidGenomicCoordsException, InvalidColourException {
        for (T x : this.getIntervalFeatureList()) {
            if (this.getTrackFormat().equals(TrackFormat.GFF)
                    || this.getTrackFormat().equals(TrackFormat.GTF)) {
                x.setGtfAttributeForName(this.gtfAttributeForName);
            } else if (this.getTrackFormat().equals(TrackFormat.BED)
                    || this.getTrackFormat().equals(TrackFormat.BIGBED)) {
                x.setBedFieldName(Integer.valueOf(this.bedFieldForName));
            }
        }
        List<String> printable = new ArrayList<>();
        int nLines = 0;
        try {
            for (List<T> listToPrint : this.stackFeatures()) {
                nLines++;
                if (nLines > this.yMaxLines) {
                    // Limit the number of lines in output
                    break;
                }
                printable.add(this.printToScreenOneLine(listToPrint));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return StringUtils.join(printable, "\n").replaceAll("\n$", "");
    }

    protected abstract T mergeFeatures(
            T a,
            T b,
            boolean screenCoords)
            throws InvalidGenomicCoordsException, InvalidColourException;


    protected List<List<T>> stackFeatures()
            throws InvalidGenomicCoordsException, IOException, InvalidColourException {

        List<T> intervals;
        List<T> flatListOfTx = this.flatListOfPrintableFeatures();
        if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.COLLAPSED)) {
            intervals = Utils.mergeIntervalFeatures(flatListOfTx, false, this::mergeFeatures);
        } else if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.ONELINE)) {
            intervals = Utils.mergeIntervalFeatures(flatListOfTx, true, this::mergeFeatures);
        } else {
            intervals = flatListOfTx;
        }

        // Make a copy of the T list. Items will be popped out as they are
        // added to individual lines.
        List<T> flatList = new ArrayList<>(intervals);

        List<List<T>> listOfLines = new ArrayList<>();
        if (flatList.isEmpty()) {
            return listOfLines;
        }
        List<T> line = new ArrayList<>();
        line.add(flatList.get(0));
        flatList.remove(0);
        listOfLines.add(line);

        while (true) {
            ArrayList<T> trToRemove = new ArrayList<T>();
            // Find a read in input whose start is greater then end of current
            for (int i = 0; i < flatList.size(); i++) {
                T intervalFeature = flatList.get(i);
                if (intervalFeature.getScreenFrom()
                        > line.get(line.size() - 1).getScreenTo()
                        + this.getGap()) { // +2 because we want some space between adjacent reads
                    listOfLines.get(listOfLines.size() - 1).add(intervalFeature); // Append to the last line.
                    trToRemove.add(intervalFeature);
                }
            } // At the end of the loop you have put in line as many reads as you can.
            for (T intervalFeature : trToRemove) {
                flatList.remove(intervalFeature);
            }
            // Create a new line, add the first T in list
            if (!flatList.isEmpty()) {
                line = new ArrayList<T>();
                line.add(flatList.get(0));
                listOfLines.add(line);
                flatList.remove(0);
            } else {
                break;
            }
        }
        return listOfLines;
    }

    /**
     * Return a string of a single line of (typically de-stacked) reads
     *
     * @throws IOException
     * @throws InvalidGenomicCoordsException
     * @throws InvalidColourException
     */
    protected String printToScreenOneLine(List<T> listToPrint)
            throws InvalidColourException {

        List<String> printable =
                new ArrayList<
                        String>(); // Each item in this list occupies a character space in the terminal.
        // NB: Each item is String not char because it might contain the ansi formatting.
        for (int i = 0; i < this.getGc().getMapping().size(); i++) { // First create empty line
            printable.add(" ");
        }
        for (T intervalFeature : listToPrint) {
            if (intervalFeature.getScreenFrom() == -1) {
                continue; // See test canProcessIndelAtWindowBoundary for how this can happen
            }
            List<FeatureChar> text = intervalFeature.getIdeogram(false, false);

            int i = 0;
            for (int j = intervalFeature.getScreenFrom(); j <= intervalFeature.getScreenTo(); j++) {
                printable.set(j, text.get(i).format(this.isNoFormat()));
                i++;
            }
        }
        return StringUtils.join(printable, "");
    }

    /**
     * List where the original records have been grouped into transcripts. If there are no
     * transcripts, just return the input feature(s) as it is. TODO: Process here also squash, merge
     * and gap?
     *
     * @throws IOException
     * @throws InvalidGenomicCoordsException
     * @throws InvalidColourException
     * @throws InvalidCommandLineException
     */
    private List<T> flatListOfPrintableFeatures()
            throws InvalidGenomicCoordsException, IOException, InvalidColourException {

        for (T x : this.getIntervalFeatureList()) {
            x.getIdeogram(true, false);
        }
        this.changeFeatureColour(this.getColourForRegex());

        List<T> flatList = new ArrayList<T>();

        List<Double> mapToScreen = this.getGc().getMapping();

        if (this.getTrackFormat().equals(TrackFormat.GFF)
                || this.getTrackFormat().equals(TrackFormat.GTF)) {

            Map<String, List<T>> tx;
            if (this.getTrackFormat().equals(TrackFormat.GFF)) {
                tx = this.groupByGFFAttribute();
            } else if (this.getTrackFormat().equals(TrackFormat.GTF)) {
                tx = this.groupByGTFAttribute();
            } else {
                throw new RuntimeException("This is not a GTF or GFF track!");
            }
            for (String txId : tx.keySet()) {
                if (txId.equals("_na_")) {
                    // Features that are not part of transcript
                    for (T x : tx.get(txId)) {
                        x.getIdeogram(false, false);
                        flatList.add(x);
                    }
                } else {
                    flatList.add(this.collapseGFFTranscript(tx.get(txId), mapToScreen));
                }
            }

        } else {
            for (T x : this.getIntervalFeatureList()) {
                flatList.add(x);
            }
        }

        Collections.sort(flatList);
        for (T x : flatList) {
            // Add name to the ideogram of each feature.
            x.getIdeogram(false, true);
        }
        return flatList;
    }

    protected String getUnformattedTitle() {

        String sq = "";
        if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.COLLAPSED)) {
            sq = "; collapsed";
        }
        String gapped = "";
        if (this.getGap() == 0) {
            gapped = "; ungapped";
        }
        String title =
                this.getTrackTag()
                        + ";"
                        + " N: "
                        + this.intervalFeatureList.size()
                        + sq
                        + gapped
                        + this.getTitleForActiveFilters();
        return title;
    }

    @Override
    protected String getTitleForActiveFilters() {
        List<String> title = new ArrayList<String>();
        if (!this.getAwk().equals(Filter.DEFAULT_AWK.getValue())) {
            title.add("awk");
        }
        if (!this.getShowRegex().pattern().equals(Filter.DEFAULT_SHOW_REGEX.getValue())
                || !this.getHideRegex().pattern().equals(Filter.DEFAULT_HIDE_REGEX.getValue())) {
            title.add("grep");
        }
        if (!title.isEmpty()) {
            return "; filters: " + title.toString();
        } else {
            return "";
        }
    }

    @Override
    public String getTitle()
            throws InvalidColourException, InvalidGenomicCoordsException, IOException {

        if (this.isHideTitle()) {
            return "";
        }
        return this.formatTitle(this.getUnformattedTitle()) + "\n";
    }

    /**
     * Searching the current chrom starting at "from" to find the *next* feature matching the given
     * string. If not found, search the other chroms, if not found restart from the beginning of the
     * current chrom until the "from" position is reached.
     *
     * @throws InvalidGenomicCoordsException
     */
    protected T findNextRegexInGenome(Pattern pattern, String chrom, int from)
            throws IOException, InvalidGenomicCoordsException {

        int startingPoint = from - 1; // -1 because tabix.query from is 0 based (seems so at least)
        List<String> chromSearchOrder = this.getChromListStartingAt(chrom);
        chromSearchOrder.add(chrom);

        for (String curChrom : chromSearchOrder) {

            TabixBigBedIterator iter = this.getReader().query(curChrom, startingPoint, Integer.MAX_VALUE);
            List<T> buffer = new ArrayList<>();

            while (true) {
                String line = iter.next();
                if (line != null) {
                    boolean matched = pattern.matcher(line).find();
                    if (matched) {
                        T x = this.createFeature(line);
                        if (x.getFrom() > startingPoint) {
                            buffer.add(x);
                        }
                    }
                }
                if (line == null || buffer.size() > 10) {
                    this.removeInvisibleFeatures(buffer);
                    if (!buffer.isEmpty()) {
                        return buffer.get(0);
                    }
                }
                if (line == null) {
                    break;
                }
            }
            startingPoint = 0;
        }
        return null; // Not found anywhere
    }

    /**
     * Return the list of known chroms with startChrom in first position. chroms: chr1 chr2 chr3 chr4
     * chr5 chr6 chr7 startChrom: chr3 return: chr3 chr4 chr5 chr6 chr7 chr1 chr2
     */
    private List<String> getChromListStartingAt(String startChrom) {

        List<String> chroms = this.getChromosomeNames();

        List<String> chromsStartingAt = new ArrayList<String>();

        int idx = chroms.indexOf(startChrom);
        if (idx == -1) {
            // If startChrom is not present at all in the bed/gtf file.
            chromsStartingAt.addAll(chroms);
            chromsStartingAt.add(0, startChrom);
        } else {
            chromsStartingAt.addAll(chroms.subList(idx, chroms.size()));
            chromsStartingAt.addAll(chroms.subList(0, idx));

            // Sanity check
            if (chroms.size() != chromsStartingAt.size()) {
                throw new RuntimeException(
                        "Error reordering chroms. Expected "
                                + chroms.size()
                                + " chroms got "
                                + chromsStartingAt.size());
            }
            if (!(chromsStartingAt.containsAll(chroms) && chroms.containsAll(chromsStartingAt))) {
                throw new RuntimeException("Error re-ordering chromosomes");
            }
        }
        if (!(chromsStartingAt.containsAll(chroms))) {
            throw new RuntimeException("Not all known chroms have been included");
        }
        return chromsStartingAt;
    }

    /**
     * Group the features in this genomic window by GFF attribute (typically a transcripts). Features
     * that don't have the attribute make each a length=1 list.
     */
    private Map<String, List<T>> groupByGFFAttribute() {

        // * First collect the IDs of the transcripts

        // Key is transcript ID e.g. ENST00001234.
        // Values is all the IntervalFeatures captured by this ID and part of a transcript.
        // I.e. their are in txFeature set,
        Map<String, List<T>> txIds = new LinkedHashMap<String, List<T>>();

        // This key:value is for records which are not part transcripts. E.g. features like "chromosome"
        // or rRNA.
        txIds.put("_na_", new ArrayList<T>());

        // Now populate the lists of values by assigning to each key the transcript records:
        for (T x : this.getIntervalFeatureList()) {

            if (FormatGTF.getTxSuperFeatures().contains(x.getFeature().toLowerCase())) {

                // Transcript feature. E.g.
                // 7 ensembl_havana mRNA 5527151 5530709 . - .
                // ID=transcript:ENST00000331789;Parent=gene:ENSG00000075624;Name=ACTB-...
                String txId = x.getGFFValueFromKey("ID");
                if (!txIds.containsKey(txId)) {
                    txIds.put(txId, new ArrayList<T>());
                }
                txIds.get(txId).add(x);
            } else if (FormatGTF.getTxSubFeatures().contains(x.getFeature().toLowerCase())) {
                // Part of transcript, e.g:
                // 7 ensembl_havana exon 5527151 5527891 . - .
                // Parent=transcript:ENST00000331789;Name=ENSE00001902654;constitutive=0;ensembl_end_pha
                String txId = x.getGFFValueFromKey("Parent");
                if (txId == null) {
                    txId = "_na_";
                }
                if (!txIds.containsKey(txId)) {
                    txIds.put(txId, new ArrayList<T>());
                }
                txIds.get(txId).add(x);
            } else {
                // Not a transcript or part thereof. E.g.
                // 7 . biological_region 5529708 5529709 0.999 - . logic_name=eponine
                txIds.get("_na_").add(x);
            }
        }
        // We don't need to return the full Map, only the list of lists (groups) would suffice.
        // However, we need to separate the group of non-trascripts (_na_ key)
        return txIds;
    }

    /** Group the features in this genomic window by GTF attribute (typically a transcripts). */
    private Map<String, List<T>> groupByGTFAttribute() {

        // * First collect the IDs of the transcripts

        // Key is transcript ID e.g. ENST00001234.
        // Values is all the IntervalFeatures captured by this ID and part of a transcript.
        // I.e. their are in txFeature set,
        Map<String, List<T>> txIds = new LinkedHashMap<String, List<T>>();

        // This key:value is for records which are not part transcripts. E.g. features like "chromosome"
        // or rRNA.
        txIds.put("_na_", new ArrayList<T>());

        // Now populate the lists of values by assigning to each key the transcript records:
        for (T x : this.getIntervalFeatureList()) {

            if (FormatGTF.getTxSuperFeatures().contains(x.getFeature().toLowerCase())
                    || FormatGTF.getTxSubFeatures().contains(x.getFeature().toLowerCase())) {
                // Transcript feature. E.g.
                // chr7 hg19_wgEncodeGencodeBasicV19 exon       5566782 5567522 0.000000 - . gene_id
                // "ENST00000331789.5"; transcript_id "ENST00000331789.5";
                String txId = x.getGFFValueFromKey("transcript_id");
                if (!txIds.containsKey(txId)) {
                    txIds.put(txId, new ArrayList<T>());
                }
                txIds.get(txId).add(x);
            } else {
                // Not a transcript or part thereof. E.g.
                // 7 . biological_region 5529708 5529709 0.999 - . logic_name=eponine
                txIds.get("_na_").add(x);
            }
        }
        // We don't need to return the full Map, only the list of lists (groups) would suffice.
        // However, we need to separate the group of non-trascripts (_na_ key)
        return txIds;
    }

    /**
     * Collapse the list of features in a single T representing the transcript. The
     * elements of txFeatures are expected to represent the entire transcript, nothing more (e.g.
     * "chromosome"). The transcript may not be biologically complete as part of it may be outside the
     * current genomic coords.
     *
     * <p>mapToScreen: Mapping of genomic coordinates to screen coordinates. This could be obtained
     * inside this method but better to pass it from outside as it can take time to get the terminal
     * window size several times.
     *
     * @throws InvalidGenomicCoordsException
     * @throws InvalidColourException
     */
    private T collapseGFFTranscript(
            List<T> txFeatures, List<Double> mapToScreen)
            throws InvalidGenomicCoordsException, InvalidColourException {

        if (txFeatures.size() == 0) {
            System.err.println("Unexpected transcript: Length zero!");
            throw new RuntimeException();
        }

        // Collect the genomic and screen coordinates of this transcript
        int gFrom = Integer.MAX_VALUE;
        int gTo = 0;
        int screenFrom = Integer.MAX_VALUE;
        int screenTo = 0;
        for (T x : txFeatures) {

            if (x.getFrom() < gFrom) {
                gFrom = x.getFrom();
            }
            if (x.getTo() > gTo) {
                gTo = x.getTo();
            }
            if (x.getScreenFrom() < screenFrom) {
                screenFrom = x.getScreenFrom();
            }
            if (x.getScreenTo() > screenTo) {
                screenTo = x.getScreenTo();
            }
        }

        T transcript =
                new T(txFeatures.get(0).getChrom(), gFrom, gTo, TrackFormat.GFF);
        transcript.setStrand(txFeatures.get(0).getStrand());
        transcript.mapToScreen(mapToScreen);

        // Now we need to prepare the ideogram
        int txIdeogramSize = screenTo - screenFrom + 1;
        List<FeatureChar> ideogram = new ArrayList<FeatureChar>(txIdeogramSize);
        for (int i = 0; i < txIdeogramSize; i++) {
            FeatureChar c = new FeatureChar();
            c.setText('-');
            ideogram.add(c); // Default character to print. Typically this should apply to introns only.
        }

        for (String txSubType : FormatGTF.getTxSubFeatures()) {

            for (T subFeature : txFeatures) {

                if (subFeature.getFeature().toLowerCase().equals(txSubType)) {
                    // Replace the featureChars in the novel transcript with the those from the individual
                    // features
                    // cccccccc                  <- subfeature#1
                    //               eeee     <- subfeature#2
                    // ---------------------- <- novel ideogram to be replaced
                    List<FeatureChar> subFeatureIdeogram = subFeature.getIdeogram(false, false);
                    int offset = subFeature.getScreenFrom() - screenFrom;
                    for (FeatureChar x : subFeatureIdeogram) {
                        ideogram.set(offset, x);
                        offset++;
                    }
                }
            }
        }
        // Now we get the name for this transcript
        String txName = "."; // Default: No name
        outerloop:
        for (String txSuperType : FormatGTF.getTxSuperFeatures()) {
            for (T x : txFeatures) {
                if (x.getFeature().toLowerCase().equals(txSuperType)) {
                    txName = x.getName();
                }
                if (txName != null && !txName.isEmpty() && !txName.equals(".")) {
                    break outerloop; // A name found, break
                }
            }
        }
        if (txName == null || txName.isEmpty() || txName.equals(".")) {
            // If a name has not been found among the superfeatures, look at the
            // individual components (exons, CDS, etc)
            outerloop:
            for (String txSuperType : FormatGTF.getTxSubFeatures()) {
                for (T x : txFeatures) {
                    if (x.getFeature().toLowerCase().equals(txSuperType)) {
                        txName = x.getName();
                    }
                    if (txName != null && !txName.isEmpty() && !txName.equals(".")) {
                        break outerloop; // A name found, break
                    }
                }
            }
        }
        transcript.setName(txName);
        transcript.setIdeogram(ideogram, false);
        return transcript;
    }

    protected List<T> getIntervalFeatureList() {
        return intervalFeatureList;
    }

    protected void setIntervalFeatureList(List<T> intervalFeatureList) {
        this.intervalFeatureList = intervalFeatureList;
    }

    @Override
    public ArrayList<String> getChromosomeNames() {
        ArrayList<String> x = new ArrayList<String>(this.getReader().getChromosomes());
        Collections.sort(x);
        return x;
    }

    TabixBigBedReader getReader() {

        if (this.bigBedReader != null && this.tabixReader != null) {
            System.err.println("You cannot have both tabix and bigBed readers set!");
            throw new RuntimeException();
        }

        if (this.bigBedReader != null) {
            return new TabixBigBedReader(this.bigBedReader);
        } else if (this.tabixReader != null) {
            return new TabixBigBedReader(this.tabixReader);
        } else {
            System.err.println("Tabix and bigBed reader both null.");
            throw new RuntimeException();
        }
    }

    protected TabixReader getTabixReader(String tabixFile) throws IOException {
        return new TabixReader(new File(tabixFile).getAbsolutePath());
    }

    /** This setter is for TrackBookmark to work. */
    protected void setTabixReader(TabixReader tabixReader) {
        this.tabixReader = tabixReader;
    }

    protected TabixReader getTabixReader() {
        return this.tabixReader;
    }

    @Override
    protected List<String> getRecordsAsStrings() {

        List<String> featureList = new ArrayList<>();
        for (T ift : intervalFeatureList) {
            featureList.add(ift.getRaw());
        }
        return featureList;
    }

    @Override
    protected void setColourForRegex(List<Argument> xcolourForRegex) {
        if (xcolourForRegex == null) {
            this.colourForRegex = null;
            return;
        } else {
            if (this.colourForRegex == null) {
                this.colourForRegex = new ArrayList<Argument>();
            }
            for (Argument p : xcolourForRegex) {
                this.colourForRegex.add(p);
            }
        }
    }

    private List<Argument> getColourForRegex() {
        return this.colourForRegex;
    }

    @Override
    protected void changeFeatureColour(List<Argument> list)
            throws InvalidColourException, IOException {
        if (list == null) {
            return;
        }

        for (Argument arg : list) {
            String regex = arg.getKey();
            String colour = arg.getArg();
            String[] rawrecs = new String[this.getIntervalFeatureList().size()];
            for (int i = 0; i < this.getIntervalFeatureList().size(); i++) {
                rawrecs[i] = this.getIntervalFeatureList().get(i).getRaw();
            }
            boolean[] matched = Utils.matchByAwkOrRegex(rawrecs, regex);
            for (int i = 0; i < matched.length; i++) {
                boolean m = matched[i];
                if (arg.isInvert()) {
                    m = !m;
                }
                if (m) {
                    for (FeatureChar f : this.getIntervalFeatureList().get(i).getIdeogram(false, false)) {
                        f.setBgColour(colour);
                        f.setFgColour(Xterm256.getContrastColour(colour));
                    }
                }
            }
        }
    }

    @Override
    public void reload()
            throws InvalidGenomicCoordsException,
            IOException,
            ClassNotFoundException,
            InvalidRecordException,
            SQLException {
        if (!Files.isSameFile(Paths.get(this.getWorkFilename()), Paths.get(this.getFilename()))) {
            TrackIntervalFeature tr = new TrackIntervalFeature(this.getFilename(), this.getGc());
            String fname = this.getWorkFilename();
            Files.move(
                    Paths.get(tr.getWorkFilename()),
                    Paths.get(fname),
                    java.nio.file.StandardCopyOption.REPLACE_EXISTING);
            Files.move(
                    Paths.get(tr.getWorkFilename() + FileExtensions.TABIX_INDEX),
                    Paths.get(fname + FileExtensions.TABIX_INDEX),
                    java.nio.file.StandardCopyOption.REPLACE_EXISTING);
        }
        this.tabixReader = this.getTabixReader(this.getWorkFilename());
        this.update();
    }

    @Override
    public void setFeatureName(String nameFieldOrAttribute) {
        if (this.getTrackFormat().equals(TrackFormat.GFF)
                || this.getTrackFormat().equals(TrackFormat.GTF)) {
            this.gtfAttributeForName = nameFieldOrAttribute;
        } else if (this.getTrackFormat().equals(TrackFormat.BED)
                || this.getTrackFormat().equals(TrackFormat.BIGBED)) {
            if (nameFieldOrAttribute.equals("-na")) {
                this.bedFieldForName = -1;
            } else {
                try {
                    this.bedFieldForName =
                            Integer.valueOf(nameFieldOrAttribute)
                                    - 1; // User's input is 1-based, convert ot 0-based
                } catch (NumberFormatException e) {
                    System.err.println("Cannot convert " + nameFieldOrAttribute + " to integer");
                    throw e;
                }
            }
        }
    }

    protected int getScoreColIdx() {
        return scoreColIdx;
    }

    protected void setScoreColIdx(int scoreColIdx)
            throws ClassNotFoundException,
            IOException,
            InvalidGenomicCoordsException,
            InvalidRecordException,
            SQLException {
        this.scoreColIdx = scoreColIdx;
        this.update();
    }

    private VCFCodec getVCFCodec() {
        return null;
    }

    @Override
    public void close() {
        if (this.tabixReader != null) {
            this.tabixReader.close();
        }
        if (this.bigBedReader != null) {
            this.bigBedReader.close();
        }
    }
}
