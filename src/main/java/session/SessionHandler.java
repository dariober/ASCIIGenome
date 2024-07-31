package session;

import static com.fasterxml.jackson.dataformat.yaml.YAMLGenerator.Feature.MINIMIZE_QUOTES;
import static com.fasterxml.jackson.dataformat.yaml.YAMLGenerator.Feature.WRITE_DOC_START_MARKER;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SequenceWriter;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import com.fasterxml.jackson.dataformat.yaml.YAMLMapper;
import com.google.common.base.Joiner;
import exceptions.SessionException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

import samTextViewer.Main;

public class SessionHandler {
  public static final File DEFAULT_SESSION_FILE =
      new File(Main.DEFAULT_ASCIIGENOME_DIR.getAbsolutePath() + File.separator + "session.yml");
  private List<Session> sessions = new ArrayList<>();
  private File sessionFile = DEFAULT_SESSION_FILE;

  public SessionHandler() {}

  public SessionHandler(File sessionFile) throws IOException, SessionException {
    if (!DEFAULT_SESSION_FILE.exists()) {
      DEFAULT_SESSION_FILE.getParentFile().mkdirs();
      DEFAULT_SESSION_FILE.createNewFile();
    }
    this.sessionFile = sessionFile.getAbsoluteFile();
    this.sessions = read(this.sessionFile);
  }

  private static List<Session> read(File sessionYamlFile) throws IOException, SessionException {
    if (!sessionYamlFile.exists() || !sessionYamlFile.canRead()) {
      throw new SessionException(
          "File '" + sessionYamlFile.getAbsolutePath() + "' does not exist or is not readable.");
    }
    InputStream yaml = Files.newInputStream(Paths.get(sessionYamlFile.getPath()));
    List<Session> sessions =
        new YAMLMapper().readValue(yaml, new TypeReference<List<Session>>() {});
    yaml.close();
    validateSessionList(sessions);
    sortSessionsByLastRead(sessions);
    return sessions;
  }

  private static void validateSessionList(List<Session> sessions) throws SessionException {
    Set<String> names = new HashSet<>();
    for (Session s : sessions) {
      if (names.contains(s.getSessionName())) {
        throw new SessionException("Duplicate session names found");
      }
      names.add(s.getSessionName());
    }
  }

  private static void addOrReplaceSession(List<Session> sessions, Session session) {
    if (sessions.contains(session)) {
      sessions.set(sessions.indexOf(session), session);
    } else {
      sessions.add(session);
    }
    sortSessionsByLastRead(sessions);
  }

  public static void saveAs(File sessionFile, Session session, String sessionName)
      throws IOException, SessionException {
    session.setSessionName(sessionName);
    List<Session> sessions = new ArrayList<>();
    if (sessionFile.exists() && sessionFile.length() > 0) {
      sessions = read(sessionFile);
    }
    addOrReplaceSession(sessions, session);
    YAMLFactory yf = new YAMLFactory();
    yf.configure(WRITE_DOC_START_MARKER, false);
    ObjectMapper mapper = new ObjectMapper(yf);
    FileOutputStream fos = new FileOutputStream(sessionFile);
    SequenceWriter sw = mapper.writerWithDefaultPrettyPrinter().writeValues(fos);
    sw.write(sessions);
    sw.close();
  }

  private static void sortSessionsByLastRead(List<Session> sessions) {
    sessions.sort(Comparator.comparing(Session::getLastRead).reversed());
  }

  public Session get(String sessionNameOrIndex) throws SessionException {
    for (Session s : this.getSessions()) {
      if (s.getSessionName().equals(sessionNameOrIndex)) {
        return s;
      }
    }
    int i;
    try {
      i = Integer.parseInt(sessionNameOrIndex);
    } catch (NumberFormatException e) {
      throw new SessionException("Cannot find session with name '" + sessionNameOrIndex + "'");
    }
    sortSessionsByLastRead(this.getSessions());
    if (i < 1 || i > this.getSessions().size())
      throw new SessionException(
          "Session index must between 1 (for most recent session) and the number of sessions in"
              + " you sessions file ("
              + this.getSessions().size()
              + ")");
    return this.getSessions().get(i - 1);
  }

  public List<Session> getSessions() {
    return sessions;
  }

  public void setSessions(List<Session> sessions) {
    this.sessions = sessions;
  }

  public File getSessionFile() {
    return sessionFile;
  }

  public void setSessionFile(File sessionFile) {
    this.sessionFile = sessionFile;
  }

  public String print(int upto, boolean mostRecentLast) throws IOException {
    List<String> out = new ArrayList<>();
    int i = 0;
    for (Session x : this.getSessions()) {
      if (upto <= i) {
        break;
      }
      YAMLFactory yf = new YAMLFactory();
      yf.configure(WRITE_DOC_START_MARKER, false);
      yf.configure(MINIMIZE_QUOTES, true);
      ObjectMapper mapper = new ObjectMapper(yf);
      String sw = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(x);
      out.add(sw);
      i += 1;
    }
    if (mostRecentLast) {
      Collections.reverse(out);
    }
    out.add("Session file: " + this.getSessionFile());
    return Joiner.on("\n").join(out);
  }
}
