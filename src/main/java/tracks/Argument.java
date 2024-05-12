package tracks;

/**
 * Class to add some feature to a generic argument. Used to define whether the regex in
 * TrackSet.setFeatureColorForRegex should change colour for feature matching or not. arg: colour to
 * use invert: invert matching?
 */
class Argument {
  private String key;
  private String arg;
  private boolean invert = false;

  Argument(String key, String arg, boolean invert) {
    this.setKey(key);
    this.setArg(arg);
    this.setInvert(invert);
  }

  protected String getArg() {
    return arg;
  }

  protected void setArg(String arg) {
    this.arg = arg;
  }

  protected boolean isInvert() {
    return invert;
  }

  protected void setInvert(boolean invert) {
    this.invert = invert;
  }

  protected String getKey() {
    return key;
  }

  protected void setKey(String key) {
    this.key = key;
  }
}
