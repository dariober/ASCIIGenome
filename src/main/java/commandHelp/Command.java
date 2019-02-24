package commandHelp;

// With some patience you should replace the hardcoded command names with these. 
public enum Command {
	featureColorForRegex("featureColorForRegex"), 
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
	bookmark("bookmark"),
	grep("grep"),
	gffNameAttr("gffNameAttr"),
	gap("gap"),
	readsAsPairs("readsAsPairs"),
	trackHeight("trackHeight"),
	colorTrack("colorTrack"),
	hideTitle("hideTitle"),
	editNames("editNames"),
	ylim("ylim"),
	dataCol("dataCol"),
	setGenome("setGenome"),
	showGenome("showGenome"),
	infoTracks("infoTracks"),
	open("open"),
	reload("reload"),
	recentlyOpened("recentlyOpened"),
	dropTracks("dropTracks"),
	orderTracks("orderTracks"),
	posHistory("posHistory"),
	history("history"),
	rpm("rpm"),
	samtools("samtools"),
	BSseq("BSseq"),
	save("save"),
	sessionSave("sessionSave");
	
    private final String commandDescription;

    private Command(String value) {
    	commandDescription = value;
    }

    public String getCmdDescr() {
        return commandDescription;
    }
}
