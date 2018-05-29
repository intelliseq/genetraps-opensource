package pl.intelliseq.genetraps.api.dx;

public interface IProcessManager {
    String runCommand(String command);

    String runCommandAndGetJobId(String command);

    String runTouch(String args);

    String runMkdir(String args);

    Object runJSONDescribe(String args);

    String runFastqc(String args);

    String runUrlFetch(String inputUrl, String sampleNumber, String... tags);

    String runBwa(String left, String right, String outputFolder);
}
