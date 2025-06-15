package commandHelp;

// With some patience you should replace the hardcoded command names with these.
public enum Command {
  featureColour("featureColour"),
  featureDisplayMode("featureDisplayMode"),
  print("print"),
  q("q"),
  h("h"),
  f("f"),
  b("b"),
  open_bracket("["),
  close_bracket("]"),
  ff("ff"),
  bb("bb"),
  zi("zi"),
  zo("zo"),
  extend("extend"),
  l("l"),
  r("r"),
  goTo("goto"),
  INT("INT"),
  plus("+"),
  minus("-"),
  p("p"),
  n("n"),
  next("next"),
  find("find"),
  seqRegex("seqRegex"),
  translate("translate"),
  bookmark("bookmark"),
  grep("grep"),
  nameForFeatures("nameForFeatures"),
  gap("gap"),
  readsAsPairs("readsAsPairs"),
  trackHeight("trackHeight"),
  colourTrack("colourTrack"),
  hideTitle("hideTitle"),
  editNames("editNames"),
  addHeader("addHeader"),
  ylim("ylim"),
  dataCol("dataCol"),
  setGenome("setGenome"),
  showGenome("showGenome"),
  infoTracks("infoTracks"),
  open("open"),
  sessionOpen("sessionOpen"),
  sessionSave("sessionSave"),
  sessionList("sessionList"),
  sessionDelete("sessionDelete"),
  reload("reload"),
  recentlyOpened("recentlyOpened"),
  dropTracks("dropTracks"),
  orderTracks("orderTracks"),
  posHistory("posHistory"),
  history("history"),
  rpm("rpm"),
  samtools("samtools"),
  BSseq("BSseq"),
  save("save");

  private final String commandDescription;

  private Command(String value) {
    commandDescription = value;
  }

  public String getCmdDescr() {
    return commandDescription;
  }
}
