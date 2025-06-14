package session;

import static org.junit.Assert.*;

import colouring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import exceptions.InvalidTrackTypeException;
import exceptions.SessionException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import org.junit.Test;
import samTextViewer.GenomicCoords;
import samTextViewer.InteractiveInputTest;
import samTextViewer.TrackProcessor;
import tracks.TrackPileup;
import tracks.TrackReads;
import tracks.TrackSet;

public class SessionHandlerTest {

  @Test
  public void canDeleteSession() throws IOException, SessionException {
    InputStream yaml = Files.newInputStream(Paths.get("test_data/session.yaml"));
    SessionHandler sh = new SessionHandler(new File("test_data/session.yaml"));
    sh.get("newSession");
    boolean deleted = sh.deleteSession("newSession");
    assertTrue(deleted);

    deleted = sh.deleteSession("newSession");
    assertFalse(deleted);
  }

  @Test
  public void testGetSessionByNameOrIndex() throws IOException, SessionException {
    InputStream yaml = Files.newInputStream(Paths.get("test_data/session.yaml"));
    SessionHandler sh = new SessionHandler(new File("test_data/session.yaml"));
    assertEquals("no-fastafile", sh.get("1").getSessionName());
  }

  @Test
  public void testSaveSession()
      throws IOException,
          InvalidGenomicCoordsException,
          SQLException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidConfigException,
          SessionException {
    new Config(null);
    TrackProcessor proc = InteractiveInputTest.gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    proc.getGenomicCoordsHistory().current().setFastaFile("test_data/chr7.fa");
    TrackPileup tr =
        new TrackPileup("test_data/ds051.actb.bam", proc.getGenomicCoordsHistory().current());
    tr.setTitleColour("red");
    proc.getTrackSet().addTrack(tr, "ds");
    proc.getTrackSet()
        .addTrack(
            new TrackReads("test_data/ds051.actb.bam", proc.getGenomicCoordsHistory().current()),
            "ds");
    new File("tmp.yml").delete();
    SessionHandler.saveAs(new File("tmp.yml"), proc.toSession(), "newSession");
    SessionHandler.saveAs(new File("tmp.yml"), proc.toSession(), "newSession2");
    SessionHandler sh = new SessionHandler(new File("tmp.yml"));
    assertEquals(2, sh.getSessions().size());
    assertEquals("newSession2", sh.get("1").getSessionName());
    assertEquals("newSession", sh.get("2").getSessionName());
    new File("tmp.yml").delete();
  }

  @Test
  public void testCanLoadGenomicCoords()
      throws IOException, InvalidGenomicCoordsException, InvalidColourException, SessionException {
    SessionHandler sh = new SessionHandler(new File("test_data/session.yaml"));
    assertEquals(4, sh.getSessions().size());

    Session session = sh.get("newSession");
    GenomicCoords gc = session.getGenome().toGenomicCoords();
    assertEquals(5566800, (long) gc.getFrom());
    assertTrue(new String(gc.getSequenceFromFasta()).startsWith("TTATTCAACTGGT"));
    assertTrue(gc.getChromIdeogram(10, true).startsWith("1"));
    assertEquals(1, gc.getSamSeqDict().size());

    session = sh.get("no-fastafile");
    gc = session.getGenome().toGenomicCoords();
    assertEquals(80, (long) gc.getFrom());
    assertEquals(25, gc.getSamSeqDict().size());
  }

  @Test
  public void testCanLoadTracks()
      throws IOException,
          InvalidTrackTypeException,
          InvalidGenomicCoordsException,
          InvalidConfigException,
          SessionException {
    new Config(null);
    SessionHandler sh = new SessionHandler(new File("test_data/session.yaml"));
    Session session = sh.get("newSession");
    TrackSet ts = session.toTrackSet();
    assertEquals("ds#1", ts.getTrackList().get(0).getTrackTag());
    assertEquals("red", ts.getTrackList().get(0).getTitleColour());
    assertEquals("ds#2", ts.getTrackList().get(1).getTrackTag());
  }

  @Test
  public void testConstructors() {
    new SessionTrack();
    new Session();
    new SessionHandler();
    new SessionGenome();
  }
}
