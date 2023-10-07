package tracks;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Joiner;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMSequenceRecord;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackBedgraph extends TrackIntervalFeature {

    public TrackBedgraph(String filename, GenomicCoords gc) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
        super(filename, gc);
        this.setTrackFormat(TrackFormat.BEDGRAPH);
    }
    
    protected TrackBedgraph(){
        
    }
    
    /* ----------- METHODS ----------- */

    /** Get values for bedgraph
     * @param intervalFeatureList 
     * @throws InvalidRecordException 
     * @throws InvalidGenomicCoordsException 
     * */
    private void bedGraphToScores(List<IntervalFeature> intervalFeatureList) throws IOException, InvalidRecordException, InvalidGenomicCoordsException{
        
        List<ScreenWiggleLocusInfo> screenWigLocInfoList= new ArrayList<ScreenWiggleLocusInfo>();
        for(int i= 0; i < getGc().getUserWindowSize(); i++){
            screenWigLocInfoList.add(new ScreenWiggleLocusInfo());
        }

        for(IntervalFeature ift : intervalFeatureList) {
            ift.mapToScreen(this.getGc().getMapping());
            for(int i= ift.getScreenFrom(); i <= ift.getScreenTo(); i++){
                screenWigLocInfoList.get(i).increment(ift.getScore());
            }
        }
        
        List<Float> screenScores= new ArrayList<Float>();
        for(ScreenWiggleLocusInfo x : screenWigLocInfoList){
            screenScores.add((float)x.getMeanScore());
        }
        this.setScreenScores(screenScores);
        return;
    }
    
    @Override
    public void update() throws IOException, InvalidRecordException, InvalidGenomicCoordsException, ClassNotFoundException, SQLException {

        this.intervalFeatureList = this.getFeaturesInInterval(
                this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());
        
        if(this.getScoreColIdx() < 4){
            System.err.println("Invalid index for bedgraph column of data value. Expected >=4. Got " + this.getScoreColIdx());
            this.scoreColIdx= 4;
            throw new InvalidRecordException();
        }
        this.bedGraphToScores(this.intervalFeatureList);
    }

    @Override
    public String printToScreen() throws InvalidColourException{

        if(this.getyMaxLines() == 0){return "";}
        
        TextProfile textProfile= new TextProfile(this.getScreenScores(), this.getyMaxLines(), this.getYLimitMin(), this.getYLimitMax());
        
        ArrayList<String> lineStrings= new ArrayList<String>();
        for(int i= (textProfile.getProfile().size() - 1); i >= 0; i--){
            List<String> xl= textProfile.getProfile().get(i);
            lineStrings.add(StringUtils.join(xl, ""));
        }

        String printable= Joiner.on("\n").join(lineStrings);
        if(!this.isNoFormat()){
            new Xterm256();
            printable= "\033[48;5;"
            + Config.get256Color(ConfigKey.background)
            + ";38;5;"
            + Xterm256.colorNameToXterm256(this.getTitleColour())
            + "m"
            + printable;
        }
        return printable;
    }
    
    @Override
    public String getTitle() throws InvalidColourException, InvalidGenomicCoordsException, IOException{

        if(this.isHideTitle()){
            return "";
        }
        
        Float[] range = Utils.range(this.getScreenScores());
        String[] rounded= Utils.roundToSignificantDigits(range[0], range[1], 2);

        String ymin= this.getYLimitMin().isNaN() ? "auto" : this.getYLimitMin().toString();
        String ymax= this.getYLimitMax().isNaN() ? "auto" : this.getYLimitMax().toString();
        
        String xtitle= this.getTrackTag() 
                + "; ylim[" + ymin + " " + ymax + "]" 
                + "; range[" + rounded[0] + " " + rounded[1] + "]";
        
        String filters = this.getTitleForActiveFilters();
        xtitle += filters;
        
        return this.formatTitle(xtitle) + "\n";
    }
    
}
