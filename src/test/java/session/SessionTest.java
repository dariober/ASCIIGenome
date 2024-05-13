package session;

import coloring.Config;
import com.esotericsoftware.yamlbeans.YamlException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;

public class SessionTest {
    @Test
    public void canReadSession() throws Exception {
        new Config(null);
        Session session = new Session("test_data/session.yaml", "last", 80);
        System.err.println(session.getTrackSet().getTrackList());
    }
}
