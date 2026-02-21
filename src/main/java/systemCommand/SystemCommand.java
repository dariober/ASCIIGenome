package systemCommand;

import com.google.errorprone.annotations.MustBeClosed;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;
import samTextViewer.Utils;
import tracks.AbstractTrack;

public class SystemCommand {

  private File tempFile;

  public SystemCommand() {}

  @SuppressWarnings("MustBeClosedChecker")
  @MustBeClosed
  public Stream<String> streamLinesThroughSystemCommand(
      Stream<String> rawRecordLines, String header, String cmd) throws IOException {
    Stream<String> rawLinesWithHeader;
    if (header == null || header.trim().isEmpty()) {
      rawLinesWithHeader = rawRecordLines;
    } else {
      List<String> l = new ArrayList<>();
      l.add(header);
      rawLinesWithHeader = Stream.concat(l.stream(), rawRecordLines);
    }

    ArrayList<String> args = new ArrayList<>();
    args.addAll(List.of("bash", "-euo", "pipefail", "-c", cmd));
    return this.execSystemCommand(rawLinesWithHeader, args);
  }

  public boolean[] passAwkFilter(String[] inLines, String awkScript) throws IOException {
    boolean[] matched = new boolean[inLines.length];
    List<String> output = new ArrayList<>();
    try (Stream<String> pass = this.streamLinesThroughAwk(Stream.of(inLines), awkScript)) {
      output = pass.toList();
    } finally {
      this.deleteTempFile();
    }
    // Check input and output. If an input line is found in output at the
    // same line number where it should be, add True else False.
    int j = 0;
    for (int i = 0; i < inLines.length; i++) {
      String inLine = inLines[i];
      if (output.size() > j) {
        String outLine = output.get(j);
        if (inLine.equals(outLine)) {
          matched[i] = true;
          j++;
        } else {
          matched[i] = false;
        }
      } else {
        matched[i] = false;
      }
    }
    return matched;
  }

  @SuppressWarnings("MustBeClosedChecker")
  @MustBeClosed
  public Stream<String> streamLinesThroughAwk(Stream<String> inLines, String awkScript)
      throws IOException {
    awkScript = awkScript.trim();

    if (awkScript.isEmpty()) {
      return inLines;
    }

    File awkTmpFile = this.prepareAwkScript(awkScript);
    String cmd = "awk -v OFS='\t' -F '\t' -f " + awkTmpFile.getAbsolutePath();

    ArrayList<String> args = new ArrayList<>();
    args.addAll(List.of("bash", "-euo", "pipefail", "-c", cmd));

    return this.execSystemCommand(inLines, args);
  }

  private File prepareAwkScript(String awkScript) throws IOException {
    String script = AbstractTrack.awkFunc + "\n" + awkScript.trim();
    File tmp = Utils.createTempFile(".asciigenome.", ".awk", true);
    Files.writeString(Path.of(tmp.getAbsolutePath()), script);
    this.setTempFile(tmp);
    return tmp;
  }

  @MustBeClosed
  private Stream<String> execSystemCommand(Stream<String> input, List<String> cmd)
      throws IOException {

    ProcessBuilder pb = new ProcessBuilder(cmd);
    Process process = pb.start();

    // Thread that feeds stdin
    Thread writerThread =
        new Thread(
            () -> {
              try (BufferedWriter writer =
                  new BufferedWriter(
                      new OutputStreamWriter(process.getOutputStream(), StandardCharsets.UTF_8))) {

                input.forEach(
                    line -> {
                      try {
                        writer.write(line);
                        writer.newLine();
                      } catch (IOException e) {
                        throw new UncheckedIOException(e);
                      }
                    });

              } catch (IOException | UncheckedIOException e) {
                process.destroyForcibly();
              }
            });

    writerThread.setDaemon(true);
    writerThread.start();

    // Thread that continuously drains stderr to avoid blocking
    StringBuilder stderrBuffer = new StringBuilder();
    Thread stderrThread =
        new Thread(
            () -> {
              try (BufferedReader errReader =
                  new BufferedReader(
                      new InputStreamReader(process.getErrorStream(), StandardCharsets.UTF_8))) {

                String line;
                while ((line = errReader.readLine()) != null) {
                  stderrBuffer.append(line).append(System.lineSeparator());
                }

              } catch (IOException e) {
                process.destroyForcibly();
              }
            });

    stderrThread.setDaemon(true);
    stderrThread.start();

    BufferedReader stdoutReader =
        new BufferedReader(new InputStreamReader(process.getInputStream(), StandardCharsets.UTF_8));

    Stream<String> stdoutStream = stdoutReader.lines();

    return stdoutStream.onClose(
        () -> {
          try {
            stdoutReader.close();

            writerThread.join();
            stderrThread.join();

            int exit = process.waitFor();
            if (exit != 0) {
              throw new RuntimeException(
                  cmd + "\n" + "Process exited with code " + exit + "\nstderr:\n" + stderrBuffer);
            }

          } catch (Exception e) {
            process.destroyForcibly();
            throw new RuntimeException(e);
          }
        });
  }

  private void setTempFile(File tempFile) {
    this.tempFile = tempFile;
  }

  public void deleteTempFile() throws IOException {
    if (this.tempFile != null) {
      Files.deleteIfExists(this.tempFile.toPath());
    }
  }
}
