package samTextViewer;

import com.esotericsoftware.yamlbeans.YamlReader;
import com.esotericsoftware.yamlbeans.YamlWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import jline.console.history.History;
import jline.console.history.MemoryHistory;
import org.apache.commons.io.FileUtils;

public class ASCIIGenomeHistory {

  private static final String DEFAULT_FILENAME =
      Main.DEFAULT_ASCIIGENOME_DIR.getAbsolutePath() + File.separator + "history.yaml";
  private final int MAX_FILES = 200;
  private final int MAX_CMDS = 2000;
  private final String fileName;
  List<String> files = new ArrayList<String>();
  List<String> positions = new ArrayList<String>();
  List<String> commands = new ArrayList<String>();
  Map<Integer, List<String>> sessions = new LinkedHashMap<Integer, List<String>>();
  private List<String> reference;

  // C O N S T R U C T O R S

  @SuppressWarnings("unchecked")
  protected ASCIIGenomeHistory(String yaml) throws IOException {
    if (yaml == null) {
      this.fileName = DEFAULT_FILENAME;
    } else {
      this.fileName = yaml;
    }

    if (yaml == null) {
      // Return an empty history object ready to be filled in.
      return;
    }

    try {
      moveOldHistoryFile();
      this.convertOldHistoryToYaml();
    } catch (Exception e) {
      //
    }

    try {
      YamlReader reader = new YamlReader(new FileReader(yaml));
      Map<String, Object> values = (HashMap<String, Object>) reader.read();
      List<String> cmd = (List<String>) values.get("commands");
      if (cmd != null) {
        this.commands = cmd;
      }
      List<String> pos = (List<String>) values.get("positions");
      if (pos != null) {
        this.positions = pos;
      }
      List<String> files = (List<String>) values.get("files");
      if (files != null) {
        this.files = files;
      }

      List<String> reference = (List<String>) values.get("reference");
      if (reference != null
          && reference.size() > 0
          && reference.get(0) != null
          && !reference.get(0).trim().isEmpty()) {
        this.reference = reference;
      }

      Map<Integer, List<String>> sessions =
          (LinkedHashMap<Integer, List<String>>) values.get("sessions");
      if (files != null) {
        this.sessions = sessions;
      }
      reader.close();
    } catch (Exception e) {
      System.err.println("Cannot read history file '" + yaml + "'");
    }
  }

  protected ASCIIGenomeHistory() throws IOException {
    this(DEFAULT_FILENAME);
  }

  // M E T H O D S

  /**
   * Read the asciigenome history file and put it a list as current history. Or return empty history
   * file does not exist or can't be read.
   */
  public History getCommandHistory() {
    History cmdHistory = new MemoryHistory();
    for (String x : this.getCommands()) {
      cmdHistory.add(x);
    }
    return cmdHistory;
  }

  public List<String> getFiles() {
    return files;
  }

  public void setFiles(List<String> files) {
    this.files = files;
  }

  public List<String> getPositions() {
    return positions;
  }

  public void setPositions(List<String> positions) {
    this.positions = positions;
  }

  private List<String> getCommands() {
    return commands;
  }

  public void setCommands(List<String> commands) {
    this.commands = commands;
  }

  protected void write(File outYaml) throws IOException {
    Map<String, Object> asciigenome_history = new HashMap<String, Object>();

    // List of commands
    List<String> lastCommands = new ArrayList<String>();
    int max_cmds = MAX_CMDS;
    for (String cmd : this.getCommands()) {
      if (max_cmds == 0) {
        break;
      }
      max_cmds--;
      lastCommands.add(cmd.trim());
    }
    asciigenome_history.put("commands", lastCommands);

    // List of files
    List<String> lastFiles = this.getFiles();
    lastFiles = lastFiles.subList(Math.max(0, lastFiles.size() - MAX_FILES), lastFiles.size());
    // Convert ArrayList#subList to List so the yaml file does show an odd data type.
    asciigenome_history.put("files", new ArrayList<String>(lastFiles));

    // Positions
    List<String> lastPos = this.getPositions();
    asciigenome_history.put("positions", lastPos);

    // Reference
    List<String> ref = this.getReference();
    if (ref != null && ref.size() > 0 && ref.get(0) != null && !ref.get(0).trim().isEmpty()) {
      asciigenome_history.put("reference", new ArrayList<String>(ref));
    }
    // Write yaml
    YamlWriter writer = new YamlWriter(new FileWriter(outYaml));
    writer.write(asciigenome_history);
    writer.close();
  }

  protected void write() throws IOException {
    File file = new File(DEFAULT_FILENAME);
    file.getParentFile().mkdirs();
    this.write(file);
  }

  public String getFileName() {
    return fileName;
  }

  /**
   * Convert old format history file to yaml. Once done, the old file is deleted.
   *
   * @throws IOException
   */
  private void convertOldHistoryToYaml() throws IOException {

    if (this.getFileName() != null && new File(this.getFileName()).exists()) {
      // If a yaml history file already exist don't do anything
      return;
    }
    File old = new File(System.getProperty("user.home") + File.separator + ".asciigenome_history");
    if (!old.exists() || this.getFileName() == null) {
      return;
    }
    List<String> oldHist = FileUtils.readLines(old, "UTF-8");
    List<String> pos = new ArrayList<String>();
    List<String> cmds = new ArrayList<String>();
    List<String> files = new ArrayList<String>();
    for (String line : oldHist) {
      line = line.trim();
      if (line.startsWith("## pos ##")) {
        pos.add(line.replaceAll("## pos ##", "").trim());
      } else if (line.startsWith("## cmd ##")) {
        cmds.add(line.replaceAll("## cmd ##", "").trim());
      } else if (line.startsWith("## file ##")) {
        files.add(line.replaceAll("## file ##", "").trim());
      } else {
        cmds.add(line);
      }
    }
    ASCIIGenomeHistory ag = new ASCIIGenomeHistory(null);
    ag.setCommands(cmds);
    ag.setFiles(files);
    ag.setPositions(pos);
    ag.write(new File(this.getFileName()));
    ag.write();
    old.delete();
  }

  public List<String> getReference() {
    return this.reference;
  }

  public void setReference(List<String> reference) {
    this.reference = reference;
  }

  private void moveOldHistoryFile() {
    File oldFile = new File(System.getProperty("user.home") + File.separator + ".asciigenome.yaml");
    if (oldFile.isFile()) {
      System.err.println(
          "Warning: The history file from older version of ASCIIGenome is going to be moved to "
              + DEFAULT_FILENAME);
      File newFile = new File(DEFAULT_FILENAME);
      newFile.getParentFile().mkdirs();
      oldFile.renameTo(newFile);
    }
  }
}
