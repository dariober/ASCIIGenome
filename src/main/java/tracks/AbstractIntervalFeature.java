package tracks;

abstract class AbstractIntervalFeature {

  private String bgColour = null;

  protected String getBgColour() {
    return bgColour;
  }

  protected void setBgColour(String bgColour) {
    this.bgColour = bgColour;
  }
}
